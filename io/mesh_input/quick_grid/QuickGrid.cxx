/*
 * QuickGrid.cxx
 *
 *  Created on: Mar 12, 2011
 *      Author: jamitch
 */

#include "QuickGrid.h"
#include "mpi.h"
#include <tr1/memory>
#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <strstream>
#include <set>
#include <cstdlib>

namespace QUICKGRID {

using namespace std;
using std::tr1::shared_ptr;
using UTILITIES::Minus;
using UTILITIES::Dot;

void print_meta_data(const QuickGridData& gridData) {
	std::strstream s;

	s << "QUICKGRID.print_meta_data(const QuickGridData& gridData)\n";
	s << "\tdimension : " << gridData.dimension << "\n";
	s << "\tglobalNumPoints : " << gridData.globalNumPoints << "\n";
	s << "\tsizeNeighborhoodList : " << gridData.sizeNeighborhoodList << "\n";
	s << "\tnumExport : " << gridData.numExport << "\n";
	std::string unpack = gridData.unPack ? "true" : "false";
	s << "\tunPack : " << unpack << "\n";

	int *ptr = gridData.neighborhoodPtr.get();
	int *neigh = gridData.neighborhood.get();
	for(size_t n=0;n<gridData.numPoints;n++,ptr++){
		int num_neigh = *neigh; neigh++;
		s << "\tNeighborhood : " << gridData.myGlobalIDs.get()[n]
		  << "; neigh ptr : " << *ptr << "; num neigh : " << num_neigh << "\n";

		for(size_t p=0;p<num_neigh;p++,neigh++){
			if(0 == p%10 && p != 0)
				s << "\n";
			s << *neigh << ", ";
		}
		s << "\n";
	}
	std::cout << s.str();

}


bool SphericalNormFunction (const double* u, const double* v, double r) {
	double dx = v[0]-u[0];
	double dy = v[1]-u[1];
	double dz = v[2]-u[2];
	return dx*dx+dy*dy+dz*dz - r*r < 0.0;
}

bool NoOpNormFunction (const double* u, const double* v, double r) {
	return true;
}

QuickGridData allocatePdGridData(size_t numCells, size_t dimension){

	// coordinates
	Array<double> X(numCells*dimension);

	// volume
	Array<double> V(numCells);

	// Global ids for cells on this processor
	Array<int> globalIds(numCells);

	// array in indices that point to neighborhood for a given localId
	Array<int> neighborhoodPtr(numCells);

	// Flag for marking points that get exported during load balance
	Array<char> exportFlag(numCells);


	// Initialize all the above data to zero
	double *xPtr = X.get();
	double *vPtr = V.get();
	int *gIdsPtr = globalIds.get();
	int *nPtr = neighborhoodPtr.get();
	char *exportFlagPtr = exportFlag.get();
	for(size_t p=0;p<numCells;p++){

		for(size_t d=0;d<dimension;d++)
			xPtr[p*dimension+d]=0;

		vPtr[p]=0;
		gIdsPtr[p]=0;
		nPtr[p]=0;
		exportFlagPtr[p]=0;
	}

	/*
	 * Neighborhood data is consistent but essentially empty.
	 * Points, Ids, volume are allocated but need to be filled with
	 * correct values
	 */
	QuickGridData gridData;

	/*
	 * Initialize neighborhood list to a consistent state
	 * 1) Set sizeNeighborhoodList=1
	 * 2) Create a new and empty neighborhood
	 * 3) Set neighborhood pointer for each point to 0
	 */
	int sizeNeighborhoodList=1;
	Array<int> neighborhoodList(sizeNeighborhoodList);
	int *neighborhood = neighborhoodList.get();
	/*
	 * number of neighbors for every point is zero
	 */
	*neighborhood = 0;

	gridData.dimension = dimension;
	gridData.globalNumPoints = 0;
	gridData.numPoints = numCells;
	gridData.sizeNeighborhoodList = sizeNeighborhoodList;
	gridData.numExport=0;
	gridData.myGlobalIDs = globalIds.get_shared_ptr();
	gridData.myX = X.get_shared_ptr();
	gridData.cellVolume = V.get_shared_ptr();
	gridData.neighborhood = neighborhoodList.get_shared_ptr();
	gridData.neighborhoodPtr = neighborhoodPtr.get_shared_ptr();
	gridData.exportFlag = exportFlag.get_shared_ptr();
	gridData.unPack = true;

	return gridData;
}

Array<double> getDiscretization(const Spec1D& spec){
	size_t numCells = spec.getNumCells();
	Array<double> ptr(numCells);
	double x0=spec.getX0();
	double cellSize=spec.getCellSize();
	double p = x0+cellSize/2.0;
	double *x = ptr.get();
	for(size_t i=0;i<numCells;p+=cellSize,i++)
		x[i]=p;
	return ptr;
}

Array<double> getDiscretization(const Spec1D& xSpec, const Spec1D& ySpec){

	Array<double> xx = getDiscretization(xSpec);
	Array<double> yy = getDiscretization(ySpec);
	double*x=xx.get();
	double*y=yy.get();

	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t numCells = nx*ny;
	Array<double> g(2*numCells);
	double* gPtr = g.get();
	double *yPtr = y;
	for(size_t j=0;j<ny;j++,yPtr++){
		double *xPtr = x;
		for(size_t i=0;i<nx;i++,xPtr++){
			*gPtr=*xPtr; gPtr++;
			*gPtr=*yPtr; gPtr++;
		}
	}

	return g;

}

Array<double> getDiscretization(const Spec1D& xSpec, const Spec1D& ySpec, const Spec1D& zSpec){
	// Set points and cells
	// note number of points is same as number of cells
	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t nz = zSpec.getNumCells();
	size_t numCells = nx*ny*nz;

	Array<double> xx = getDiscretization(xSpec);
	Array<double> yy = getDiscretization(ySpec);
	Array<double> zz = getDiscretization(zSpec);
	double*x=xx.get();
	double*y=yy.get();
	double*z=zz.get();

	size_t dimension=3;
	Array<double> X(numCells*dimension);
	double *XPtr = X.get();
	double point[3]={0.0,0.0,0.0};
	for(size_t k=0;k<nz;k++){
		point[2]=z[k];
		for(size_t j=0;j<ny;j++){
			point[1]=y[j];
			for(size_t i=0;i<nx;i++){
				point[0]=x[i];
				for(size_t p=0;p<3;p++,XPtr++)
					*XPtr = point[p];
			}
		}
	}

	return X;
}

Array<double> getDiscretization(const SpecRing2D& spec) {

	// Compute set of radii for each ring
	Spec1D rSpec(spec.getNumRings(), spec.getrI(), spec.getRingThickness());
	Array<double> rPtr = getDiscretization(rSpec);
	Vector3D c = spec.getCenter();

	size_t numRays = spec.getNumRays();
	size_t numRings = spec.getNumRings();

	Spec1D thetaSpec(numRays, 0.0, 2.0*M_PI);
	Array<double> thPtr = getDiscretization(thetaSpec);


	// Total number of cells
	size_t numCells = spec.getNumCells();

	// At each point store (x,y,z=0) + (c[0]+c[1]+c[2])
	Array<double> gPtr(3*numCells);
	double *g = gPtr.get();

	// Outer loop on rays
	double *theta = thPtr.get();

	for(size_t ny=0;ny<numRays;ny++,theta++){
		// Loop over rings
		double *r = rPtr.get();
		for(size_t nx=0;nx<numRings;nx++,r++){
			*g = (*r)*cos(*theta) + c[0]; g++;
			*g = (*r)*sin(*theta) + c[1]; g++;
			*g = c[2];                    g++;
		}
	}
	return gPtr;
}



Array<double>  getDiscretization(const SpecRing2D& spec, const Spec1D& axisSpec){
	Array<double> ptsPtr = getDiscretization(spec);
	Array<double> zPtr = getDiscretization(axisSpec);

	size_t nz = axisSpec.getNumCells();
	size_t numCellsRing = spec.getNumCells();
	size_t numCells = numCellsRing*nz;
	Array<double> gPtr(3*numCells);
	double *g = gPtr.get();

	// Loop over VTK points and set
	double *z = zPtr.get();
	for(size_t n=0;n<nz;n++,z++){
		double *pts = ptsPtr.get();
		for(size_t i=0;i<numCellsRing;i++){
			// this does x and y and z
			*g = *pts; g++, pts++;
			*g = *pts; g++, pts++;
			*g = *z;   g++; pts++;

		}
	}

	return gPtr;
}

Horizon Spec1D::getCellHorizon(double h) const {
	int n = getCellNeighborhoodSize(h);
	return Horizon(n,this->numCells);
}

RingHorizon Spec1D::getRingCellHorizon(double h,double ringRadius) const {
	int n = getCellNeighborhoodSize(h,ringRadius);
	return RingHorizon(n,this->numCells);
}

/*
 * Input:  double horizon
 * Output: integer horizon -- This is the one-sided size of the
 * neighborhood expressed in the number of cells which define
 * the discretization
 */
size_t Spec1D::getCellNeighborhoodSize(double horizon, double ringRadius) const {
	double dx = getCellSize()*ringRadius;
	double dCx = horizon/dx;
	size_t nCx = static_cast<size_t>(dCx);
//	std::cout << "PdQPointSet1d::getCellNeighborhoodSize" << std::endl;
//	std::cout << "\tincoming horizon = " << h << "; ratio dCx = h/dx = " << dCx << "; nCx = " << nCx << std::endl;
	if(abs(dCx - nCx) >= .5 )
		nCx += 1;
//	std::cout << "Final nCx = " << nCx << std::endl;

	// this handles the pathology when the horizon > 1.5 * xLength
	if(nCx > numCells)
     nCx = numCells;

	return nCx;
}


SpecRing2D::SpecRing2D(Vector3D center, double innerRadius, double outerRadius, int numRings)
: c(center), rI(innerRadius), r0(outerRadius), numRings(numRings), numRays(0), numCells(0) {

	// Compute set of radii for each ring
	double ringThickness = r0-rI;
	double dr=ringThickness/numRings;
	// Use an average radius
	double R = (rI+r0)/2;
	// Try to make cells with equal length sides -- compute angular increment
	double dTheta = dr/R;
	numRays = (size_t)(2 * M_PI / dTheta) + 1;
	numCells = numRays * this->numRings;

}

double SpecRing2D::getRingThickness() const { return abs(r0-rI); }








QuickGridMeshGenerationIterator::QuickGridMeshGenerationIterator(size_t numProcs, NormFunctionPointer norm)
: iteratorProc(0),
  numProcs(numProcs),
  neighborHoodNorm(norm)
{}

std::pair<Cell3D,QuickGridData> QuickGridMeshGenerationIterator::beginIterateProcs(QuickGridData& pdGridDataProc0) {
	iteratorProc = 0;
	Cell3D cellLocator(0,0,0);
	std::pair<Cell3D,QuickGridData> returnVal = computePdGridData(iteratorProc, cellLocator, pdGridDataProc0, neighborHoodNorm);
	iteratorProc++;
	return returnVal;
}

std::pair<Cell3D,QuickGridData>  QuickGridMeshGenerationIterator::nextProc(Cell3D cellLocator, QuickGridData& pdGridDataProcN) {
	std::pair<Cell3D,QuickGridData> returnVal = computePdGridData(iteratorProc, cellLocator, pdGridDataProcN, neighborHoodNorm);
	iteratorProc++;
	return returnVal;
}

std::vector<size_t> QuickGridMeshGenerationIterator::getNumCellsPerProcessor(size_t globalNumCells, size_t numProcs) {
	// compute cellsPerProc
	std::vector<size_t> cellsPerProc;
	int numCellsPerProc = globalNumCells/numProcs;
	int numCellsLastProc = numCellsPerProc + globalNumCells % numProcs;
	cellsPerProc  = vector<size_t>(numProcs,numCellsPerProc);
	cellsPerProc[numProcs-1] = numCellsLastProc;
	return cellsPerProc;
}

TensorProduct3DMeshGenerator::TensorProduct3DMeshGenerator
(
		size_t numProc,
		double radiusHorizon,
		const Spec1D& xSpec,
		const Spec1D& ySpec,
		const Spec1D& zSpec,
		NormFunctionPointer norm
)
:
QuickGridMeshGenerationIterator(numProc,norm),
horizonRadius(radiusHorizon), globalNumberOfCells(0), specs(3,xSpec), H(3,xSpec.getCellHorizon(radiusHorizon))
{

	globalNumberOfCells = xSpec.getNumCells()*ySpec.getNumCells()*zSpec.getNumCells();

	// Need to correct for the initial value set at default constructor of vector
	// specs[0] = xSpec;
	specs[1] = ySpec; H[1] = ySpec.getCellHorizon(radiusHorizon);
	specs[2] = zSpec; H[2] = zSpec.getCellHorizon(radiusHorizon);

	cellsPerProc = QuickGridMeshGenerationIterator::getNumCellsPerProcessor(globalNumberOfCells,numProc);

	// This is for computing the tensor product space
	xx = getDiscretization(xSpec);
	yy = getDiscretization(ySpec);
	zz = getDiscretization(zSpec);


}

QuickGridData TensorProduct3DMeshGenerator::allocatePdGridData() const {

	// Number of cells on this processor -- number of cells on last
	//   processor is always <= numCellsPerProc
	size_t numCellsAllocate = cellsPerProc[getNumProcs()-1];

	size_t dimension = 3;
	return QUICKGRID::allocatePdGridData(numCellsAllocate, dimension);

}

std::pair<Cell3D,QuickGridData> TensorProduct3DMeshGenerator::computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm) const {

//	std::cout << "CellsPerProcessor3D::computePdGridData proc = " << proc << std::endl;
	Spec1D xSpec = specs[0];
	Spec1D ySpec = specs[1];
	Spec1D zSpec = specs[2];
	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t nz = zSpec.getNumCells();

	// Each cell has the same volume
	double cellVolume = xSpec.getCellSize()*ySpec.getCellSize()*zSpec.getCellSize();

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	Horizon hY = H[1];
	Horizon hZ = H[2];


	// Number of cells on this processor
	size_t myNumCells = cellsPerProc[proc];

	// Coordinates used for tensor product space
	const double*x=xx.get();
	const double*y=yy.get();
	const double*z=zz.get();

	// Discretization on this processor
	pdGridData.dimension = 3;
	pdGridData.globalNumPoints = nx*ny*nz;
	pdGridData.numPoints = myNumCells;
	int *gidsPtr = pdGridData.myGlobalIDs.get();
	double *XPtr = pdGridData.myX.get();
	double *volPtr = pdGridData.cellVolume.get();
	// allocate neighborhood for incoming data since each one of these has a different length
	int sizeNeighborhoodList = getSizeNeighborList(proc, cellLocator);
	pdGridData.sizeNeighborhoodList = sizeNeighborhoodList;
	Array<int> neighborhoodList(sizeNeighborhoodList);
	pdGridData.neighborhood = neighborhoodList.get_shared_ptr();
	int *neighborPtr = neighborhoodList.get();
	int updateSizeNeighborhoodList = 0;

	// Compute discretization on processor
	size_t kStart = cellLocator.k;
	size_t jStart = cellLocator.j;
	size_t iStart = cellLocator.i;

	size_t cell = 0;
	double point[3]={0.0,0.0,0.0};
	double pointNeighbor[3] = {0.0,0.0,0.0};
	for(size_t k=kStart;k<nz && cell<myNumCells;k++){
		// kStart = 0; this one doesn't matter
		cellLocator.k = k;
		point[2]=z[k];
		size_t zStart = hZ.start(k);
		size_t nCz = hZ.numCells(k);

		for(size_t j=jStart;j<ny && cell<myNumCells;j++){
			jStart=0; // begin at jStart only the first time
			cellLocator.j = j;
			point[1]=y[j];
			size_t yStart = hY.start(j);
			size_t nCy = hY.numCells(j);

			for(size_t i=iStart;i<nx && cell<myNumCells;i++,cell++,cellLocator.i++){
				iStart=0; // begin at iStart only the first time
				// global id
				*gidsPtr = i + j * nx + k * nx * ny; gidsPtr++;

				point[0]=x[i];
				// copy point data

				for(size_t p=0;p<3;p++,XPtr++)
					*XPtr = point[p];

				// copy cell volume
				*volPtr = cellVolume; volPtr++;

				size_t xStart = hX.start(i);
				size_t nCx = hX.numCells(i);
				// number of cells in this (i,j,k) neighborhood
//				*neighborPtr = computeNumNeighbors(i,j,k); neighborPtr++;
				// Compute neighborhood
				int *numNeigh = neighborPtr; neighborPtr++;
				int countNeighbors = 0;

				for(size_t kk=zStart;kk<nCz+zStart;kk++){
					for(size_t jj=yStart;jj<nCy+yStart;jj++){
						for(size_t ii=xStart;ii<nCx+xStart;ii++){
							// skip this cell since its its own neighbor
							if(ii == i && jj == j && kk == k) continue;

							pointNeighbor[0] = x[ii];
							pointNeighbor[1] = y[jj];
							pointNeighbor[2] = z[kk];
							if(norm(pointNeighbor,point,horizonRadius)){
								int globalId = ii + jj * nx + kk * nx * ny;
								*neighborPtr = globalId; neighborPtr++;
								countNeighbors++;
							}
						}
					}
				}

				*numNeigh = countNeighbors;
				updateSizeNeighborhoodList += countNeighbors;

			}

			if(cellLocator.i == nx){
				cellLocator.i = 0;
				if(cell==myNumCells){
					// current j holds current last row which needs to be incremented
					cellLocator.j += 1;
					if(cellLocator.j == ny){
						cellLocator.j = 0;
						// current k holds current last row which needs to be incremented
						cellLocator.k += 1;
					}
				}
			}

		}
	}

	// post process and compute neighborhoodPtr
	pdGridData.sizeNeighborhoodList = updateSizeNeighborhoodList+myNumCells;
	int *neighPtr = pdGridData.neighborhoodPtr.get();
	int *neighborListPtr = neighborhoodList.get();
	int sum=0;
	for(int n=0;n<myNumCells;n++){
		neighPtr[n]=sum;
		int numCells = *neighborListPtr; neighborListPtr += (numCells+1);
		sum += (numCells+1);
	}

	return std::pair<Cell3D,QuickGridData>(cellLocator,pdGridData);
}

size_t TensorProduct3DMeshGenerator::getSizeNeighborList(size_t proc, Cell3D cellLocator) const {


	Spec1D xSpec = specs[0];
	Spec1D ySpec = specs[1];
	Spec1D zSpec = specs[2];
	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t nz = zSpec.getNumCells();

	/*
	 * Compacted list of neighbors
	 * 1) for i = 0,..(numCells-1)
	 * 2)   numCellNeighborhood(i)
	 * 3)   neighorhood(i) -- does not include "i"
	 */
	// for each cell store numCells +  neighborhood
	// int size =  numCells    +                     sum(numCells(i),i);
	size_t numCells = cellsPerProc[proc];
	size_t size=numCells;

	// Perform sum over all neighborhoods
	size_t kStart = cellLocator.k;
	size_t jStart = cellLocator.j;
	size_t iStart = cellLocator.i;
	size_t cell = 0;

	for(size_t k=kStart;k<nz;k++){
		// kStart = 0; this one doesn't matter

		for(size_t j=jStart;j<ny;j++){
			jStart=0; // begin at jStart only the first time

			for(size_t i=iStart;i<nx && cell<numCells;i++,cell++){
				// begin at iStart only the first time
				iStart = 0;
				size += computeNumNeighbors(i,j,k);
			}
		}
	}

	return size;
}

size_t TensorProduct3DMeshGenerator::computeNumNeighbors(size_t i, size_t j, size_t k) const {

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	Horizon hY = H[1];
	Horizon hZ = H[2];

	size_t nCx = hX.numCells(i);
	size_t nCy = hY.numCells(j);
	size_t nCz = hZ.numCells(k);

	/*
	 * Number of cells in this (i,j,k) neighborhood
	 * Note that '1' is subtracted to remove "this cell"
	 */

	return nCx * nCy * nCz - 1;
}



TensorProductCylinderMeshGenerator::TensorProductCylinderMeshGenerator
(
		size_t nProcs,
		double radiusHorizon,
		const SpecRing2D& rSpec,
		const Spec1D& axisSpec,
		NormFunctionPointer norm
) :
QuickGridMeshGenerationIterator(nProcs,norm),
horizonRadius(radiusHorizon),
globalNumberOfCells(0),
ringSpec(rSpec),
specs(3,axisSpec),
ringHorizon(0,0),
H(3,axisSpec.getCellHorizon(radiusHorizon))
{
	// Number of cells in cylinder
	globalNumberOfCells = rSpec.getNumCells()*axisSpec.getNumCells();
	// Set up tensor product specs
	// radial direction
	specs[0] = ringSpec.getRaySpec();
	// theta direction
	specs[1] = ringSpec.getRingSpec();
	// cylinder axis -- z direction
	// already set on initialization
	//specs[2] = axisSpec;

	// need to reset ringHorizon;
	/*
	 * Uses average radius to compute cell size;
	 * If careful attention is given to this everything will be fine.
	 * But for special cases -- such as when size of elements along the
	 * different axes are not approximately the same this may produce
	 * unexpected results although they would be consistent.
	 */
	ringHorizon = specs[1].getRingCellHorizon(radiusHorizon,(rSpec.getrI()+rSpec.getr0())/2.0);

	// Now set up the horizon objects -- note that the ring horizon is special and needs special treatment
	H[0] = specs[0].getCellHorizon(radiusHorizon); // r
	H[1] = specs[1].getCellHorizon(radiusHorizon); // theta
	// H[2] = specs[2].getCellHorizon(radiusHorizon); // z -- this one is already set on initialization

	// compute cellsPerProc
	cellsPerProc = QuickGridMeshGenerationIterator::getNumCellsPerProcessor(globalNumberOfCells,nProcs);

	rPtr = getDiscretization(specs[0]);
	thetaPtr = getDiscretization(specs[1]);
	zPtr = getDiscretization(specs[2]);


}

std::pair<Cell3D,QuickGridData> TensorProductCylinderMeshGenerator::computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm) const {

	/**
	 * This routine attempts to make as close of an analogy to CellsPerProcessor3D::computePdGridData(...)
	 * as makes sense:
	 *    x <--> r
	 *    y <--> theta
	 *    z <--> z
	 */
//	std::cout << "CellsPerProcessorCylinder::computePdGridData proc = " << proc << std::endl;
	Spec1D xSpec = specs[0]; // r
	Spec1D ySpec = specs[1]; // theta
	Spec1D zSpec = specs[2]; // z
	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int nz = zSpec.getNumCells();

	// Each cell has the same volume arc size in radians --
	//   Cells have different volumes bases upon "r"
	double dr = xSpec.getCellSize();
	double dz = zSpec.getCellSize();
	double cellRads = ySpec.getCellSize();
	double cellVolume;

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	Horizon hY = H[1];
	Horizon hZ = H[2];

	// Number of cells on this processor
	int myNumCells = cellsPerProc[proc];

	// Coordinates used for tensor product space
	const double*x=rPtr.get();
	const double*y=thetaPtr.get();
	const double*z=zPtr.get();

	// Discretization on this processor
	pdGridData.dimension = 3;
	pdGridData.globalNumPoints = nx*ny*nz;
	pdGridData.numPoints = myNumCells;
	int *gidsPtr = pdGridData.myGlobalIDs.get();
	double *XPtr = pdGridData.myX.get();
	double *volPtr = pdGridData.cellVolume.get();
	// allocate neighborhood for incoming data since each one of these has a different length
	int sizeNeighborhoodList = getSizeNeighborList(proc, cellLocator);
	pdGridData.sizeNeighborhoodList = sizeNeighborhoodList;
	Array<int> neighborhoodList(sizeNeighborhoodList);
	pdGridData.neighborhood = neighborhoodList.get_shared_ptr();
	int *neighborPtr = neighborhoodList.get();
	int updateSizeNeighborhoodList = 0;

	// Compute discretization on processor
	int kStart = cellLocator.k;
	int jStart = cellLocator.j;
	int iStart = cellLocator.i;

	size_t cell = 0;
	double point[3]={0.0,0.0,0.0};
	double pointNeighbor[3] = {0.0,0.0,0.0};
	for(size_t k=kStart;k<nz && cell<myNumCells;k++){
		// kStart = 0; this one doesn't matter
		cellLocator.k = k;
		point[2]=z[k];
		size_t zStart = hZ.start(k);
		size_t nCz = hZ.numCells(k);

		for(size_t j=jStart;j<ny && cell<myNumCells;j++){
			jStart=0; // begin at jStart only the first time
			cellLocator.j = j;
			double theta = y[j];
			double c = cos(theta);
			double s = sin(theta);
			for(size_t i=iStart;i<nx && cell<myNumCells;i++,cell++,cellLocator.i++){
				iStart=0; // begin at iStart only the first time
				// global id
				*gidsPtr = i + j * nx + k * nx * ny; gidsPtr++;

				double r = x[i];
				point[0] = r * c;
				point[1] = r * s;
				// copy point data
				for(size_t p=0;p<3;p++,XPtr++)
					*XPtr = point[p];

				// copy cell volume
				cellVolume = r*dr*cellRads*dz;
				*volPtr = cellVolume; volPtr++;

				size_t xStart = hX.start(i);
				size_t nCx = hX.numCells(i);
				// number of cells in this (i,j,k) neighborhood
//				*neighborPtr = computeNumNeighbors(i,j,k); neighborPtr++;
				// Compute neighborhood
				int *numNeigh = neighborPtr; neighborPtr++;
				size_t countNeighbors = 0;
				for(size_t kk=zStart;kk<nCz+zStart;kk++){
					RingHorizon::RingHorizonIterator yHorizonIter = ringHorizon.horizonIterator(j);

					while(yHorizonIter.hasNextCell()){
						size_t jj = yHorizonIter.nextCell();
						for(size_t ii=xStart;ii<nCx+xStart;ii++){
							// skip this cell since its its own neighbor
							if(ii == i && jj == j && kk == k) continue;

							// Neighborhood norm calculation
							double r = x[ii];
							double theta = y[jj];
							double c = cos(theta);
							double s = sin(theta);
							pointNeighbor[0] = r * c;
							pointNeighbor[1] = r * s;
							pointNeighbor[2] = z[kk];

							if(norm(pointNeighbor,point,horizonRadius)){
								int globalId = ii + jj * nx + kk * nx * ny;
								*neighborPtr = globalId; neighborPtr++;
								countNeighbors++;
							}
						}
					}
				}
				*numNeigh = countNeighbors;
				updateSizeNeighborhoodList += countNeighbors;
			}

			if(cellLocator.i == nx){
				cellLocator.i = 0;
				if(cell==myNumCells){
					// current j holds current last row which needs to be incremented
					cellLocator.j += 1;
					if(cellLocator.j == ny){
						cellLocator.j = 0;
						// current k holds current last row which needs to be incremented
						cellLocator.k += 1;
					}
				}
			}

		}
	}
//
	// post process and compute neighborhoodPtr
	pdGridData.sizeNeighborhoodList = updateSizeNeighborhoodList+myNumCells;
	int *neighPtr = pdGridData.neighborhoodPtr.get();
	int *neighborListPtr = neighborhoodList.get();
	int sum=0;
	for(int n=0;n<myNumCells;n++){
		neighPtr[n]=sum;
		int numCells = *neighborListPtr; neighborListPtr += (numCells+1);
		sum += (numCells+1);
	}

	return std::pair<Cell3D,QuickGridData>(cellLocator,pdGridData);
}

/**
 * This function is distinct in that special care must be taken to
 * properly account for a cylinder
 */
size_t TensorProductCylinderMeshGenerator::computeNumNeighbors(size_t i, size_t j, size_t k) const {

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	// Horizon hY = H[1]; // cylinder uses ringHorizon
	RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(j);
	Horizon hZ = H[2];

	size_t nCx = hX.numCells(i);
	size_t nCy = hIter.numCells();
	size_t nCz = hZ.numCells(k);
	/*
	 * Number of cells in this (i,j,k) neighborhood
	 * Note that '1' is subtracted to remove "this cell"
	 */

	return nCx * nCy * nCz - 1;
}

size_t TensorProductCylinderMeshGenerator::getSizeNeighborList(size_t proc, Cell3D cellLocator) const {


	Spec1D xSpec = specs[0];
	Spec1D ySpec = specs[1];
	Spec1D zSpec = specs[2];
	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t nz = zSpec.getNumCells();

	/*
	 * Compacted list of neighbors
	 * 1) for i = 0,..(numCells-1)
	 * 2)   numCellNeighborhood(i)
	 * 3)   neighorhood(i) -- does not include "i"
	 */
	// for each cell store numCells +  neighborhood
	// int size =  numCells    +                     sum(numCells(i),i);
	size_t numCells = cellsPerProc[proc];
	size_t size=numCells;

	// Perform sum over all neighborhoods
	size_t kStart = cellLocator.k;
	size_t jStart = cellLocator.j;
	size_t iStart = cellLocator.i;
	size_t cell = 0;

	for(size_t k=kStart;k<nz;k++){
		// kStart = 0; this one doesn't matter

		for(size_t j=jStart;j<ny;j++){
			jStart=0; // begin at jStart only the first time

			for(size_t i=iStart;i<nx && cell<numCells;i++,cell++){
				// begin at iStart only the first time
				iStart = 0;
				size += computeNumNeighbors(i,j,k);
			}
		}
	}

	return size;
}

QuickGridData TensorProductCylinderMeshGenerator::allocatePdGridData() const {

	// Number of cells on this processor -- number of cells on last
	//   processor is always <= numCellsPerProc
	size_t numCellsAllocate = cellsPerProc[getNumProcs()-1];

	size_t dimension = 3;
	return QUICKGRID::allocatePdGridData(numCellsAllocate, dimension);

}

TensorProductSolidCylinder::TensorProductSolidCylinder
(
		size_t nProcs,
		double radius,
		size_t numRings,
		double zStart,
		double cylinderLength
) :
QuickGridMeshGenerationIterator(nProcs),
coreRadius(0),
outerRadius(radius),
coreCellVolume(0),
globalNumberOfCells(0),
specs(3,Spec1D(1,0,1))
{

	/*
	 * Compute inner radius of inner most ring
	 */
	double R = radius;
	size_t nR = numRings;
	double rI = R/nR/sqrt(M_PI);
	coreRadius = rI;

	/*
	 * Create radial spec
	 */
	Spec1D raySpec(nR,rI,R-rI);

	/*
	 * Compute dTheta -- ring secgment size -- attempts to create cells that have
	 * approximately the same area
	 */
	double dTheta = 2.0 * sqrt(M_PI) / ( sqrt(M_PI) * nR + 1);

	/*
	 * Create ring spec
	 */
	size_t numRays = (int)(2 * M_PI / dTheta) + 1;
	Spec1D ringSpec(numRays, 0.0, 2.0*M_PI);

	/*
	 * Now get approximate cell size along length of cylinder
	 */
	double cellSize = (R-rI)/nR;

	/*
	 * Now create discretization along axis of cylinder
	 */
	size_t numCellsAxis = (size_t)(cylinderLength/cellSize)+1;
	Spec1D axisSpec(numCellsAxis,zStart,cylinderLength);

	/*
	 * Compute "coreCellVolume"
	 */
	coreCellVolume = M_PI * coreRadius * coreRadius * axisSpec.getCellSize();

	/*
	 * Add in 1 cell for the "core"
	 */
	globalNumberOfCells = (numRings * numRays + 1) * numCellsAxis;

//	cout << "TensorProductSolidCylinder; numRings, numRays, numCellsAxis, globalNumberOfCells = "
//			<< numRings << ", "
//			<< numRays << ", "
//			<< numCellsAxis << ", "
//			<< globalNumberOfCells
//			<<  endl;

	// compute cellsPerProc
	cellsPerProc = QuickGridMeshGenerationIterator::getNumCellsPerProcessor(globalNumberOfCells,nProcs);

//	cout << "cellsPerProc.size() = " << cellsPerProc.size() << endl;
//	for(std::size_t i=0;i<cellsPerProc.size();i++)
//		cout << "proc, cellsPerProc = " << i << ", " << cellsPerProc[i] << endl;

	/*
	 * Save specs
	 */
	specs[0] = raySpec;
	specs[1] = ringSpec;
	specs[2] = axisSpec;

	rPtr = getDiscretization(specs[0]);
	thetaPtr = getDiscretization(specs[1]);
	zPtr = getDiscretization(specs[2]);


}

QuickGridData TensorProductSolidCylinder::allocatePdGridData() const {

	// Number of cells on this processor -- number of cells on last
	//   processor is always <= numCellsPerProc
	size_t numCellsAllocate = cellsPerProc[getNumProcs()-1];

	size_t dimension = 3;
	return QUICKGRID::allocatePdGridData(numCellsAllocate, dimension);

}

std::pair<Cell3D,QuickGridData> TensorProductSolidCylinder::computePdGridData(size_t proc, Cell3D cellLocator, QuickGridData& pdGridData, NormFunctionPointer norm) const {

	/**
	 * This routine attempts to make as close of an analogy to CellsPerProcessor3D::computePdGridData(...)
	 * as makes sense:
	 *    x <--> r
	 *    y <--> theta
	 *    z <--> z
	 */
//	std::cout << "CellsPerProcessorCylinder::computePdGridData proc = " << proc << std::endl;
	Spec1D xSpec = specs[0]; // r
	Spec1D ySpec = specs[1]; // theta
	Spec1D zSpec = specs[2]; // z
	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t nz = zSpec.getNumCells();

	// Each cell has the same volume arc size in radians --
	//   Cells have different volumes bases upon "r"
	double dr = xSpec.getCellSize();
	double dz = zSpec.getCellSize();
	double cellRads = ySpec.getCellSize();
	double cellVolume;

	// Number of cells on this processor
	size_t myNumCells = cellsPerProc[proc];


	// Coordinates used for tensor product space
	const double*x=rPtr.get();
	const double*y=thetaPtr.get();
	const double*z=zPtr.get();

	// Discretization on this processor
	pdGridData.dimension = 3;
	size_t numRings = nx;
	size_t numRays = ny;
	/*
	 * Note that the number of cells is the tensor product of nx, ny, and nz plus the number of core cells with is nz
	 */
	pdGridData.globalNumPoints = numRings*numRays*nz + nz;
	pdGridData.numPoints = myNumCells;
	int *gidsPtr = pdGridData.myGlobalIDs.get();
	double *XPtr = pdGridData.myX.get();
	double *volPtr = pdGridData.cellVolume.get();

//	cout << "\tCellsPerProcessorCylinder::computePdGridData; proc, myNumCells = " << proc << ", " << myNumCells  << endl;
//	cout << "CellsPerProcessorCylinder::computePdGridData globalNumPoints = " << pdGridData.globalNumPoints  << endl;
	// Compute discretization on processor
	size_t kStart = cellLocator.k;
	size_t jStart = cellLocator.j;
	size_t iStart = cellLocator.i;

	size_t cell = 0;
	double point[3]={0.0,0.0,0.0};
	for(size_t k=kStart;k<nz && cell<myNumCells;k++){
		// kStart = 0; this one doesn't matter
		cellLocator.k = k;
		point[2]=z[k];

		/*
		 * We only add the core cell at the start of every x-y plane of points otherwise NOT
		 */
		if(0==jStart && 0==iStart){
			/*
			 * CORE CELL
			 */
			/*
			 * Place core cell at start of cells in x-y disk
			 * NOTE that number of cells in disk = nx * ny + 1 core cell
			 */
			*gidsPtr = k * ( nx * ny + 1); gidsPtr++;
			/*
			 * increment cell counter
			 */
			cell++;

			*volPtr = coreCellVolume; volPtr++;
			point[0] = 0;
			point[1] = 0;
			// copy point data
			for(int p=0;p<3;p++,XPtr++)
				*XPtr = point[p];

//			cout << "\tproc, core cell gid  = " << proc <<  ", " << k * ( nx * ny + 1) << endl;
		}


		for(size_t j=jStart;j<ny && cell<myNumCells;j++){
			jStart=0; // begin at jStart only the first time
			cellLocator.j = j;
			double theta = y[j];
			double c = cos(theta);
			double s = sin(theta);
			for(size_t i=iStart;i<nx && cell<myNumCells;i++,cell++,cellLocator.i++){
				iStart=0; // begin at iStart only the first time
				/*
				 * global id: note that global ids start with on offset of "1" to include core cell at start
				 */
				*gidsPtr = 1 + i + j * nx + k * ( nx * ny + 1); gidsPtr++;

//				cout << "\tproc, cell gid  = " << proc <<  ", " << (1 + i + j * nx + k * ( nx * ny + 1)) << endl;
				double r = x[i];
				point[0] = r * c;
				point[1] = r * s;
				// copy point data
				for(int p=0;p<3;p++,XPtr++)
					*XPtr = point[p];

				// copy cell volume
				cellVolume = r*dr*cellRads*dz;
				*volPtr = cellVolume; volPtr++;

			}

			if(cellLocator.i == nx){
				cellLocator.i = 0;
				if(cell==myNumCells){
					// current j holds current last row which needs to be incremented
					cellLocator.j += 1;
					if(cellLocator.j == ny){
						cellLocator.j = 0;
						// current k holds current last row which needs to be incremented
						cellLocator.k += 1;
					}
				}
			}

		}
	}

	/*
	 * Allocate 1 entry for number of neighbors -- note that neighborPtr points to 0 for every point (PdQuickGrid::allocatePdGridData)
	 * and therefore the length of the neighborhood list for every point is zero
	 */
	pdGridData.sizeNeighborhoodList=1;
	Array<int> neighborhoodList(pdGridData.sizeNeighborhoodList);
	pdGridData.neighborhood = neighborhoodList.get_shared_ptr();
	int *neighborhood = neighborhoodList.get();
	/*
	 * number of neighbors for every point is zero
	 */
	*neighborhood = 0;
	return std::pair<Cell3D,QuickGridData>(cellLocator,pdGridData);
}


/**
 * This function produces an unbalanced discretization although for some geometries it
 * may not be too bad
 */
QuickGridData getDiscretization(size_t rank, QuickGridMeshGenerationIterator &cellIter)
{
	MPI_Status status;
	int ack = 0;
	int ackTag = 0;
	int numPointsTag=1;
	int idsTag = 2;
	int coordinatesTag = 3;
	int globalNumPointsTag=4;
	int sizeNeighborhoodListTag=5;
	int neighborhoodTag=6;
	int volumeTag=7;
	int neighborhoodPtrTag=8;
	QuickGridData gridData;
	int dimension = cellIter.getDimension();

	if(0 == rank){

		QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
		QuickGridData  pdGridDataProcN = cellIter.allocatePdGridData();
		std::pair<Cell3D,QuickGridData> p0Data = cellIter.beginIterateProcs(pdGridDataProc0);
		gridData = p0Data.second;
		Cell3D nextCellLocator = p0Data.first;


		while(cellIter.hasNextProc()){
			int proc = cellIter.proc();
			std::pair<Cell3D,QuickGridData> data = cellIter.nextProc(nextCellLocator,pdGridDataProcN);
			QuickGridData gridData = data.second;
			nextCellLocator = data.first;


			// Need to send this data to proc
			int globalNumPoints = gridData.globalNumPoints;
			int numPoints = gridData.numPoints;
			int sizeNeighborhoodList = gridData.sizeNeighborhoodList;
			shared_ptr<int> gIds = gridData.myGlobalIDs;
			shared_ptr<double> X = gridData.myX;
			shared_ptr<double> V = gridData.cellVolume;
			shared_ptr<int> neighborhood = gridData.neighborhood;
			shared_ptr<int> neighborhoodPtr = gridData.neighborhoodPtr;

			MPI_Send(&numPoints, 1, MPI_INT, proc, numPointsTag, MPI_COMM_WORLD);
			MPI_Recv(&ack, 1, MPI_INT, proc, ackTag, MPI_COMM_WORLD, &status);
			MPI_Send(&globalNumPoints, 1, MPI_INT, proc, globalNumPointsTag, MPI_COMM_WORLD);
			MPI_Send(&sizeNeighborhoodList, 1, MPI_INT, proc, sizeNeighborhoodListTag, MPI_COMM_WORLD);
			MPI_Send(gIds.get(), numPoints, MPI_INT, proc, idsTag, MPI_COMM_WORLD);
			MPI_Send(X.get(), dimension*numPoints, MPI_DOUBLE, proc, coordinatesTag, MPI_COMM_WORLD);
			MPI_Send(V.get(),numPoints, MPI_DOUBLE, proc, volumeTag, MPI_COMM_WORLD);
			MPI_Send(neighborhood.get(), sizeNeighborhoodList, MPI_INT, proc, neighborhoodTag, MPI_COMM_WORLD);
			MPI_Send(neighborhoodPtr.get(), numPoints, MPI_INT, proc, neighborhoodPtrTag, MPI_COMM_WORLD);

		}
	    /* signal all procs it is OK to go on */
	    ack = 0;
	    for(int proc=1;proc<cellIter.getNumProcs();proc++){
	      MPI_Send(&ack, 1, MPI_INT, proc, ackTag, MPI_COMM_WORLD);
	    }

	}
	else {
		// Receive data from processor 0
		// Create this procs 'GridData'
		int numPoints=0;
		int globalNumPoints = 0;
		int sizeNeighborhoodList = 0;
		MPI_Recv(&numPoints, 1, MPI_INT, 0, numPointsTag, MPI_COMM_WORLD, &status);
		ack = 0;
		if (numPoints > 0){
			QuickGridData 	gData = QUICKGRID::allocatePdGridData(numPoints,dimension);
			std::tr1::shared_ptr<double> g=gData.myX;
			std::tr1::shared_ptr<double> cellVolume=gData.cellVolume;
			std::tr1::shared_ptr<int> gIds=gData.myGlobalIDs;
			std::tr1::shared_ptr<int> neighborhoodPtr=gData.neighborhoodPtr;
			MPI_Send(&ack, 1, MPI_INT, 0, ackTag, MPI_COMM_WORLD);
			MPI_Recv(&globalNumPoints, 1, MPI_INT, 0, globalNumPointsTag, MPI_COMM_WORLD, &status);
			MPI_Recv(&sizeNeighborhoodList, 1, MPI_INT, 0, sizeNeighborhoodListTag, MPI_COMM_WORLD, &status);
			Array<int> neighborhood(sizeNeighborhoodList);
			MPI_Recv(gIds.get(), numPoints, MPI_INT, 0, idsTag, MPI_COMM_WORLD, &status);
			MPI_Recv(g.get(), dimension*numPoints, MPI_DOUBLE, 0, coordinatesTag, MPI_COMM_WORLD, &status);
			MPI_Recv(cellVolume.get(), numPoints, MPI_DOUBLE, 0,volumeTag, MPI_COMM_WORLD, &status);
			MPI_Recv(neighborhood.get(), sizeNeighborhoodList, MPI_INT, 0, neighborhoodTag, MPI_COMM_WORLD, &status);
			MPI_Recv(neighborhoodPtr.get(), numPoints, MPI_INT, 0, neighborhoodPtrTag, MPI_COMM_WORLD, &status);

			gData.dimension = dimension;
			gData.globalNumPoints = globalNumPoints;
			gData.numPoints = numPoints;
			gData.sizeNeighborhoodList = sizeNeighborhoodList;
			gData.neighborhood=neighborhood.get_shared_ptr();

			gridData = gData;

		}
		else if (numPoints == 0){
			MPI_Send(&ack, 1, MPI_INT, 0, ackTag, MPI_COMM_WORLD);
		}
		else{
			MPI_Finalize();
			std::exit(1);
		}

	    MPI_Recv(&ack, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
	    if (ack < 0){
	      MPI_Finalize();
	      std::exit(1);
	    }

	}
	return gridData;
}

}
