/*! \file PdQuickGrid.cxx */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "PdQuickGrid.h"
#include <tr1/memory>
#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
#include <set>
#include "Epetra_Comm.h"

namespace PdQuickGrid {

using namespace std;
using std::tr1::shared_ptr;


/*
 * Declaration for private functions
 */
std::pair<int, std::tr1::shared_ptr<int> > getSharedGlobalIds(const PdGridData& pdGridData);
const Epetra_BlockMap getOverlap(int dimension, int numShared, int*shared, int numOwned, const int* owned, const Epetra_Comm& comm);

bool SphericalNormFunction (const double* u, const double* v, double r) {
	double dx = v[0]-u[0];
	double dy = v[1]-u[1];
	double dz = v[2]-u[2];
	return dx*dx+dy*dy+dz*dz - r*r < 0.0;
}

bool NoOpNormFunction (const double* u, const double* v, double r) {
	return true;
}

std::tr1::shared_ptr<int> getLocalOwnedIds(const PdGridData& gridData, const Epetra_BlockMap& overlapMap){
	shared_ptr<int> localIds(new int[gridData.numPoints],PdQuickGrid::Deleter<int>());
	int *lIds = localIds.get();
	int *end = localIds.get()+gridData.numPoints;
	int *gIds = gridData.myGlobalIDs.get();
	for(; lIds != end;lIds++, gIds++)
		*lIds = overlapMap.LID(*gIds);
	return localIds;
}

std::tr1::shared_ptr<int> getLocalNeighborList(const PdGridData& gridData, const Epetra_BlockMap& overlapMap){
	shared_ptr<int> localNeighborList(new int[gridData.sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	int *localNeig = localNeighborList.get();
	int *neighPtr = gridData.neighborhoodPtr.get();
	int *neigh = gridData.neighborhood.get();
	for(int p=0;p<gridData.numPoints;p++){
		int ptr = neighPtr[p];
		int numNeigh = neigh[ptr];
		localNeig[ptr]=numNeigh;
		for(int n=1;n<=numNeigh;n++){
			int gid = neigh[ptr+n];
			int localId = overlapMap.LID(gid);
			localNeig[ptr+n] = localId;
		}
	}
	return localNeighborList;
}

const Epetra_BlockMap getOwnedMap(const Epetra_Comm& comm,const PdGridData& gridData, int ndf){
	int numShared=0;
	int *sharedPtr=NULL;
	int numOwned = gridData.numPoints;
	const int *ownedPtr = gridData.myGlobalIDs.get();
	return getOverlap(ndf, numShared,sharedPtr,numOwned,ownedPtr,comm);
}

/*
 * This function is private
 */
const Epetra_BlockMap getOverlap(int ndf, int numShared, int*shared, int numOwned,const  int* owned, const Epetra_Comm& comm){

	int numPoints = numShared+numOwned;
	shared_ptr<int> ids(new int[numPoints],PdQuickGrid::Deleter<int>());
	int *ptr = ids.get();

	for(int j=0;j<numOwned;j++,ptr++)
		*ptr=owned[j];

	for(int j=0;j<numShared;j++,ptr++)
		*ptr=shared[j];

	return Epetra_BlockMap(-1,numPoints, ids.get(),ndf, 0,comm);

}

const Epetra_BlockMap getOverlapMap(const Epetra_Comm& comm, const PdGridData& gridData, int ndf){
	std::pair<int, std::tr1::shared_ptr<int> > sharedPair = PdQuickGrid::getSharedGlobalIds(gridData);
	std::tr1::shared_ptr<int> sharedPtr = sharedPair.second;
	int numShared = sharedPair.first;
	int *shared = sharedPtr.get();
	int *owned = gridData.myGlobalIDs.get();
	int numOwned = gridData.numPoints;
	return getOverlap(ndf,numShared,shared,numOwned,owned,comm);
}

std::pair<int, std::tr1::shared_ptr<int> > getSharedGlobalIds(const PdGridData& gridData){
	std::set<int> ownedIds(gridData.myGlobalIDs.get(),gridData.myGlobalIDs.get()+gridData.numPoints);
	std::set<int> shared;
	int *neighPtr = gridData.neighborhoodPtr.get();
	int *neigh = gridData.neighborhood.get();
	std::set<int>::const_iterator ownedIdsEnd = ownedIds.end();
	for(int p=0;p<gridData.numPoints;p++){
		int ptr = neighPtr[p];
		int numNeigh = neigh[ptr];
		for(int n=1;n<=numNeigh;n++){
			int id = neigh[ptr+n];
			/*
			 * look for id in owned points
			 */
			if(ownedIdsEnd == ownedIds.find(id)){
				/*
				 * add this point to shared
				 */
				shared.insert(id);
			}
		}
	}

	// Copy set into shared ptr
	shared_ptr<int> sharedGlobalIds(new int[shared.size()],PdQuickGrid::Deleter<int>());
	int *sharedPtr = sharedGlobalIds.get();
	set<int>::iterator it;
	for ( it=shared.begin() ; it != shared.end(); it++, sharedPtr++ )
		*sharedPtr = *it;

	return std::pair<int, std::tr1::shared_ptr<int> >(shared.size(),sharedGlobalIds);
}

std::tr1::shared_ptr<double> getDiscretization(const PdQRing2d& spec) {

	// Compute set of radii for each ring
	PdQPointSet1d rSpec(spec.getNumRings(), spec.getrI(), spec.getRingThickness());
	shared_ptr<double> rPtr = PdQuickGrid::getDiscretization(rSpec);
	valarray<double> c = spec.getCenter();

	int numRays = spec.getNumRays();
	int numRings = spec.getNumRings();

	PdQPointSet1d thetaSpec(numRays, 0.0, 2.0*M_PI);
	shared_ptr<double> thPtr = PdQuickGrid::getDiscretization(thetaSpec);


	// Total number of cells
	int numCells = spec.getNumCells();

	// At each point store (x,y,z=0) + (c[0]+c[1]+c[2])
	shared_ptr<double> gPtr(new double[3*numCells],PdQuickGrid::Deleter<double>());
	double *g = gPtr.get();

	// Outer loop on rays
	double *theta = thPtr.get();

	for(int ny=0;ny<numRays;ny++,theta++){
		// Loop over rings
		double *r = rPtr.get();
		for(int nx=0;nx<numRings;nx++,r++){
			*g = (*r)*cos(*theta) + c[0]; g++;
			*g = (*r)*sin(*theta) + c[1]; g++;
			*g = c[2];                    g++;
		}
	}
	return gPtr;
}


std::tr1::shared_ptr<double>  getDiscretization(const PdQRing2d& spec, const PdQPointSet1d& axisSpec){
	std::tr1::shared_ptr<double> ptsPtr = getDiscretization(spec);
	std::tr1::shared_ptr<double> zPtr = getDiscretization(axisSpec);

	int nz = axisSpec.getNumCells();
	int numCellsRing = spec.getNumCells();
	int numCells = numCellsRing*nz;
	shared_ptr<double> gPtr(new double[3*numCells],PdQuickGrid::Deleter<double>());
	double *g = gPtr.get();

	// Loop over VTK points and set
	double *z = zPtr.get();
	for(int n=0;n<nz;n++,z++){
		double *pts = ptsPtr.get();
		for(int i=0;i<numCellsRing;i++){
			// this does x and y and z
			*g = *pts; g++, pts++;
			*g = *pts; g++, pts++;
			*g = *z;   g++; pts++;

		}
	}

	return gPtr;
}

shared_ptr<double> getDiscretization(const PdQPointSet1d& spec){
	int numCells = spec.getNumCells();
	shared_ptr<double> ptr(new double[numCells],PdQuickGrid::Deleter<double>());
	double x0=spec.getX0();
	double cellSize=spec.getCellSize();
	double p = x0+cellSize/2.0;
	double *x = ptr.get();
	for(int i=0;i<numCells;p+=cellSize,i++)
		x[i]=p;
	return ptr;
}

shared_ptr<double> getDiscretization(const PdQPointSet1d& xSpec, const PdQPointSet1d& ySpec){
	std::tr1::shared_ptr<double> xx = getDiscretization(xSpec);
	std::tr1::shared_ptr<double> yy = getDiscretization(ySpec);
	double*x=xx.get();
	double*y=yy.get();

	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int numCells = nx*ny;
	std::tr1::shared_ptr<double> g(new double[2*numCells],PdQuickGrid::Deleter<double>());
	double* gPtr = g.get();
	double *yPtr = y;
	for(int j=0;j<ny;j++,yPtr++){
		double *xPtr = x;
		for(int i=0;i<nx;i++,xPtr++){
			*gPtr=*xPtr; gPtr++;
			*gPtr=*yPtr; gPtr++;
		}
	}

	return g;
}

Horizon PdQPointSet1d::getCellHorizon(double h) const {
	int n = getCellNeighborhoodSize(h);
	return Horizon(n,this->numCells);
}

RingHorizon PdQPointSet1d::getRingCellHorizon(double h,double ringRadius) const {
	int n = getCellNeighborhoodSize(h,ringRadius);
	return RingHorizon(n,this->numCells);
}

/*
 * Input:  double horizon
 * Output: integer horizon -- This is the one-sided size of the
 * neighborhood expressed in the number of cells which define
 * the discretization
 */
int PdQPointSet1d::getCellNeighborhoodSize(double horizon, double ringRadius) const {
	double dx = getCellSize()*ringRadius;
	double dCx = horizon/dx;
	int nCx = static_cast<int>(dCx);
//	std::cout << "PdQPointSet1d::getCellNeighborhoodSize" << std::endl;
//	std::cout << "\tincoming horizon = " << h << "; ratio dCx = h/dx = " << dCx << "; nCx = " << nCx << std::endl;
	if(abs(dCx - nCx) >= .5 )
		nCx += 1;
//	std::cout << "Final nCx = " << nCx << std::endl;
	return nCx;
}


PdQRing2d::PdQRing2d(valarray<double> center, double innerRadius, double outerRadius, int numRings)
: c(center), rI(innerRadius), r0(outerRadius), numRings(numRings), numRays(0), numCells(0) {

	// Compute set of radii for each ring
	double ringThickness = r0-rI;
	double dr=ringThickness/numRings;
	// Use an average radius
	double R = (rI+r0)/2;
	// Try to make cells with equal length sides -- compute angular increment
	double dTheta = dr/R;
	numRays = (int)(2 * M_PI / dTheta) + 1;
	numCells = numRays * this->numRings;

}

double PdQRing2d::getRingThickness() const { return abs(r0-rI); }

PdGridData allocatePdGridData(int numCells, int dimension){

	// coordinates
	shared_ptr<double> X(new double[numCells*dimension],PdQuickGrid::Deleter<double>());

	// volume
	shared_ptr<double> V(new double[numCells],PdQuickGrid::Deleter<double>());

	// Global ids for cells on this processor
	shared_ptr<int> globalIds(new int[numCells],PdQuickGrid::Deleter<int>());

	// array in indices that point to neighborhood for a given localId
	shared_ptr<int> neighborhoodPtr(new int[numCells],PdQuickGrid::Deleter<int>());

	// Flag for marking points that get exported during load balance
	shared_ptr<char> exportFlag(new char[numCells],PdQuickGrid::Deleter<char>());


	// Initialize all the above data to zero
	double *xPtr = X.get();
	double *vPtr = V.get();
	int *gIdsPtr = globalIds.get();
	int *nPtr = neighborhoodPtr.get();
	char *exportFlagPtr = exportFlag.get();
	for(int p=0;p<numCells;p++){

		for(int d=0;d<dimension;d++)
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
	PdGridData gridData;

	/*
	 * Initialize neighborhood list to a consistent state
	 * 1) Set sizeNeighborhoodList=1
	 * 2) Create a new and empty neighborhood
	 * 3) Set neighborhood pointer for each point to 0
	 */
	int sizeNeighborhoodList=1;
	shared_ptr<int> neighborhoodList(new int[sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
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
	gridData.myGlobalIDs = globalIds;
	gridData.myX = X;
	gridData.cellVolume = V;
	gridData.neighborhood = neighborhoodList;
	gridData.neighborhoodPtr = neighborhoodPtr;
	gridData.exportFlag = exportFlag;
	gridData.unPack = true;

	return gridData;
}

shared_ptr<double> getDiscretization(const PdQPointSet1d& xSpec, const PdQPointSet1d& ySpec, const PdQPointSet1d& zSpec){
	// Set points and cells
	// note number of points is same as number of cells
	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int nz = zSpec.getNumCells();
	int numCells = nx*ny*nz;

	std::tr1::shared_ptr<double> xx = getDiscretization(xSpec);
	std::tr1::shared_ptr<double> yy = getDiscretization(ySpec);
	std::tr1::shared_ptr<double> zz = getDiscretization(zSpec);
	double*x=xx.get();
	double*y=yy.get();
	double*z=zz.get();

	int dimension=3;
	shared_ptr<double> X(new double[numCells*dimension],PdQuickGrid::Deleter<double>());
	double *XPtr = X.get();
	double point[3]={0.0,0.0,0.0};
	for(int k=0;k<nz;k++){
		point[2]=z[k];
		for(int j=0;j<ny;j++){
			point[1]=y[j];
			for(int i=0;i<nx;i++){
				point[0]=x[i];
				for(int p=0;p<3;p++,XPtr++)
					*XPtr = point[p];
			}
		}
	}

	return X;
}

PdQuickGridMeshGenerationIterator::PdQuickGridMeshGenerationIterator(int numProcs, NormFunctionPointer norm)
: iteratorProc(0),
  numProcs(numProcs),
  neighborHoodNorm(norm)
{}

std::pair<Cell3D,PdGridData> PdQuickGridMeshGenerationIterator::beginIterateProcs(PdGridData& pdGridDataProc0) {
	iteratorProc = 0;
	Cell3D cellLocator(0,0,0);
	std::pair<Cell3D,PdGridData> returnVal = computePdGridData(iteratorProc, cellLocator, pdGridDataProc0, neighborHoodNorm);
	iteratorProc++;
	return returnVal;
}

std::pair<Cell3D,PdGridData>  PdQuickGridMeshGenerationIterator::nextProc(Cell3D cellLocator, PdGridData& pdGridDataProcN) {
	std::pair<Cell3D,PdGridData> returnVal = computePdGridData(iteratorProc, cellLocator, pdGridDataProcN, neighborHoodNorm);
	iteratorProc++;
	return returnVal;
}

std::vector<int> PdQuickGridMeshGenerationIterator::getNumCellsPerProcessor(int globalNumCells, int numProcs) {
	// compute cellsPerProc
	std::vector<int> cellsPerProc;
	int numCellsPerProc = globalNumCells/numProcs;
	int numCellsLastProc = numCellsPerProc + globalNumCells % numProcs;
	cellsPerProc  = vector<int>(numProcs,numCellsPerProc);
	cellsPerProc[numProcs-1] = numCellsLastProc;
	return cellsPerProc;
}

TensorProduct3DMeshGenerator::TensorProduct3DMeshGenerator
(
		int numProc,
		double radiusHorizon,
		const PdQPointSet1d& xSpec,
		const PdQPointSet1d& ySpec,
		const PdQPointSet1d& zSpec,
		NormFunctionPointer norm
)
:
PdQuickGridMeshGenerationIterator(numProc,norm),
horizonRadius(radiusHorizon), globalNumberOfCells(0), specs(3,xSpec), H(3,xSpec.getCellHorizon(radiusHorizon))
{

	globalNumberOfCells = xSpec.getNumCells()*ySpec.getNumCells()*zSpec.getNumCells();

	// Need to correct for the initial value set at default constructor of vector
	// specs[0] = xSpec;
	specs[1] = ySpec; H[1] = ySpec.getCellHorizon(radiusHorizon);
	specs[2] = zSpec; H[2] = zSpec.getCellHorizon(radiusHorizon);

	cellsPerProc = PdQuickGridMeshGenerationIterator::getNumCellsPerProcessor(globalNumberOfCells,numProc);

	// This is for computing the tensor product space
	xx = getDiscretization(xSpec);
	yy = getDiscretization(ySpec);
	zz = getDiscretization(zSpec);


}

PdGridData TensorProduct3DMeshGenerator::allocatePdGridData() const {

	// Number of cells on this processor -- number of cells on last
	//   processor is always <= numCellsPerProc
	int numCellsAllocate = cellsPerProc[getNumProcs()-1];

	int dimension = 3;
	return PdQuickGrid::allocatePdGridData(numCellsAllocate, dimension);

}

std::pair<Cell3D,PdGridData> TensorProduct3DMeshGenerator::computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm) const {

//	std::cout << "CellsPerProcessor3D::computePdGridData proc = " << proc << std::endl;
	PdQPointSet1d xSpec = specs[0];
	PdQPointSet1d ySpec = specs[1];
	PdQPointSet1d zSpec = specs[2];
	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int nz = zSpec.getNumCells();

	// Each cell has the same volume
	double cellVolume = xSpec.getCellSize()*ySpec.getCellSize()*zSpec.getCellSize();

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	Horizon hY = H[1];
	Horizon hZ = H[2];


	// Number of cells on this processor
	int myNumCells = cellsPerProc[proc];

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
	shared_ptr<int> neighborhoodList(new int[sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	pdGridData.neighborhood = neighborhoodList;
	int *neighborPtr = neighborhoodList.get();
	int updateSizeNeighborhoodList = 0;

	// Compute discretization on processor
	int kStart = cellLocator.k;
	int jStart = cellLocator.j;
	int iStart = cellLocator.i;

	int cell = 0;
	double point[3]={0.0,0.0,0.0};
	double pointNeighbor[3] = {0.0,0.0,0.0};
	for(int k=kStart;k<nz && cell<myNumCells;k++){
		// kStart = 0; this one doesn't matter
		cellLocator.k = k;
		point[2]=z[k];
		int zStart = hZ.start(k);
		int nCz = hZ.numCells(k);

		for(int j=jStart;j<ny && cell<myNumCells;j++){
			jStart=0; // begin at jStart only the first time
			cellLocator.j = j;
			point[1]=y[j];
			int yStart = hY.start(j);
			int nCy = hY.numCells(j);

			for(int i=iStart;i<nx && cell<myNumCells;i++,cell++,cellLocator.i++){
				iStart=0; // begin at iStart only the first time
				// global id
				*gidsPtr = i + j * nx + k * nx * ny; gidsPtr++;

				point[0]=x[i];
				// copy point data

				for(int p=0;p<3;p++,XPtr++)
					*XPtr = point[p];

				// copy cell volume
				*volPtr = cellVolume; volPtr++;

				int xStart = hX.start(i);
				int nCx = hX.numCells(i);
				// number of cells in this (i,j,k) neighborhood
//				*neighborPtr = computeNumNeighbors(i,j,k); neighborPtr++;
				// Compute neighborhood
				int *numNeigh = neighborPtr; neighborPtr++;
				int countNeighbors = 0;

				for(int kk=zStart;kk<nCz+zStart;kk++){
					for(int jj=yStart;jj<nCy+yStart;jj++){
						for(int ii=xStart;ii<nCx+xStart;ii++){
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

	return std::pair<Cell3D,PdGridData>(cellLocator,pdGridData);
}

int TensorProduct3DMeshGenerator::getSizeNeighborList(int proc, Cell3D cellLocator) const {


	PdQPointSet1d xSpec = specs[0];
	PdQPointSet1d ySpec = specs[1];
	PdQPointSet1d zSpec = specs[2];
	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int nz = zSpec.getNumCells();

	/*
	 * Compacted list of neighbors
	 * 1) for i = 0,..(numCells-1)
	 * 2)   numCellNeighborhood(i)
	 * 3)   neighorhood(i) -- does not include "i"
	 */
	// for each cell store numCells +  neighborhood
	// int size =  numCells    +                     sum(numCells(i),i);
	int numCells = cellsPerProc[proc];
	int size=numCells;

	// Perform sum over all neighborhoods
	int kStart = cellLocator.k;
	int jStart = cellLocator.j;
	int iStart = cellLocator.i;
	int cell = 0;

	for(int k=kStart;k<nz;k++){
		// kStart = 0; this one doesn't matter

		for(int j=jStart;j<ny;j++){
			jStart=0; // begin at jStart only the first time

			for(int i=iStart;i<nx && cell<numCells;i++,cell++){
				// begin at iStart only the first time
				iStart = 0;
				size += computeNumNeighbors(i,j,k);
			}
		}
	}

	return size;
}

int TensorProduct3DMeshGenerator::computeNumNeighbors(int i, int j, int k) const {

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	Horizon hY = H[1];
	Horizon hZ = H[2];

	int nCx = hX.numCells(i);
	int nCy = hY.numCells(j);
	int nCz = hZ.numCells(k);

	/*
	 * Number of cells in this (i,j,k) neighborhood
	 * Note that '1' is subtracted to remove "this cell"
	 */

	return nCx * nCy * nCz - 1;
}

TensorProductCylinderMeshGenerator::TensorProductCylinderMeshGenerator
(
		int nProcs,
		double radiusHorizon,
		const PdQRing2d& rSpec,
		const PdQPointSet1d& axisSpec,
		NormFunctionPointer norm
) :
PdQuickGridMeshGenerationIterator(nProcs,norm),
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
	cellsPerProc = PdQuickGridMeshGenerationIterator::getNumCellsPerProcessor(globalNumberOfCells,nProcs);

	rPtr = getDiscretization(specs[0]);
	thetaPtr = getDiscretization(specs[1]);
	zPtr = getDiscretization(specs[2]);


}

std::pair<Cell3D,PdGridData> TensorProductCylinderMeshGenerator::computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm) const {

	/**
	 * This routine attempts to make as close of an analogy to CellsPerProcessor3D::computePdGridData(...)
	 * as makes sense:
	 *    x <--> r
	 *    y <--> theta
	 *    z <--> z
	 */
//	std::cout << "CellsPerProcessorCylinder::computePdGridData proc = " << proc << std::endl;
	PdQPointSet1d xSpec = specs[0]; // r
	PdQPointSet1d ySpec = specs[1]; // theta
	PdQPointSet1d zSpec = specs[2]; // z
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
	shared_ptr<int> neighborhoodList(new int[sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	pdGridData.neighborhood = neighborhoodList;
	int *neighborPtr = neighborhoodList.get();
	int updateSizeNeighborhoodList = 0;

	// Compute discretization on processor
	int kStart = cellLocator.k;
	int jStart = cellLocator.j;
	int iStart = cellLocator.i;

	int cell = 0;
	double point[3]={0.0,0.0,0.0};
	double pointNeighbor[3] = {0.0,0.0,0.0};
	for(int k=kStart;k<nz && cell<myNumCells;k++){
		// kStart = 0; this one doesn't matter
		cellLocator.k = k;
		point[2]=z[k];
		int zStart = hZ.start(k);
		int nCz = hZ.numCells(k);

		for(int j=jStart;j<ny && cell<myNumCells;j++){
			jStart=0; // begin at jStart only the first time
			cellLocator.j = j;
			double theta = y[j];
			double c = cos(theta);
			double s = sin(theta);
			for(int i=iStart;i<nx && cell<myNumCells;i++,cell++,cellLocator.i++){
				iStart=0; // begin at iStart only the first time
				// global id
				*gidsPtr = i + j * nx + k * nx * ny; gidsPtr++;

				double r = x[i];
				point[0] = r * c;
				point[1] = r * s;
				// copy point data
				for(int p=0;p<3;p++,XPtr++)
					*XPtr = point[p];

				// copy cell volume
				cellVolume = r*dr*cellRads*dz;
				*volPtr = cellVolume; volPtr++;

				int xStart = hX.start(i);
				int nCx = hX.numCells(i);
				// number of cells in this (i,j,k) neighborhood
//				*neighborPtr = computeNumNeighbors(i,j,k); neighborPtr++;
				// Compute neighborhood
				int *numNeigh = neighborPtr; neighborPtr++;
				int countNeighbors = 0;
				for(int kk=zStart;kk<nCz+zStart;kk++){
					RingHorizon::RingHorizonIterator yHorizonIter = ringHorizon.horizonIterator(j);

					while(yHorizonIter.hasNextCell()){
						int jj = yHorizonIter.nextCell();
						for(int ii=xStart;ii<nCx+xStart;ii++){
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

	return std::pair<Cell3D,PdGridData>(cellLocator,pdGridData);
}

/**
 * This function is distinct in that special care must be taken to
 * properly account for a cylinder
 */
int TensorProductCylinderMeshGenerator::computeNumNeighbors(int i, int j, int k) const {

	// Horizon in each of the coordinate directions
	Horizon hX = H[0];
	// Horizon hY = H[1]; // cylinder uses ringHorizon
	RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(j);
	Horizon hZ = H[2];

	int nCx = hX.numCells(i);
	int nCy = hIter.numCells();
	int nCz = hZ.numCells(k);
	/*
	 * Number of cells in this (i,j,k) neighborhood
	 * Note that '1' is subtracted to remove "this cell"
	 */

	return nCx * nCy * nCz - 1;
}

int TensorProductCylinderMeshGenerator::getSizeNeighborList(int proc, Cell3D cellLocator) const {


	PdQPointSet1d xSpec = specs[0];
	PdQPointSet1d ySpec = specs[1];
	PdQPointSet1d zSpec = specs[2];
	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int nz = zSpec.getNumCells();

	/*
	 * Compacted list of neighbors
	 * 1) for i = 0,..(numCells-1)
	 * 2)   numCellNeighborhood(i)
	 * 3)   neighorhood(i) -- does not include "i"
	 */
	// for each cell store numCells +  neighborhood
	// int size =  numCells    +                     sum(numCells(i),i);
	int numCells = cellsPerProc[proc];
	int size=numCells;

	// Perform sum over all neighborhoods
	int kStart = cellLocator.k;
	int jStart = cellLocator.j;
	int iStart = cellLocator.i;
	int cell = 0;

	for(int k=kStart;k<nz;k++){
		// kStart = 0; this one doesn't matter

		for(int j=jStart;j<ny;j++){
			jStart=0; // begin at jStart only the first time

			for(int i=iStart;i<nx && cell<numCells;i++,cell++){
				// begin at iStart only the first time
				iStart = 0;
				size += computeNumNeighbors(i,j,k);
			}
		}
	}

	return size;
}

PdGridData TensorProductCylinderMeshGenerator::allocatePdGridData() const {

	// Number of cells on this processor -- number of cells on last
	//   processor is always <= numCellsPerProc
	int numCellsAllocate = cellsPerProc[getNumProcs()-1];

	int dimension = 3;
	return PdQuickGrid::allocatePdGridData(numCellsAllocate, dimension);

}

TensorProductSolidCylinder::TensorProductSolidCylinder
(
		int nProcs,
		double radius,
		int numRings,
		double zStart,
		double cylinderLength
) :
PdQuickGridMeshGenerationIterator(nProcs),
coreRadius(0),
outerRadius(radius),
coreCellVolume(0),
globalNumberOfCells(0),
specs(3,PdQPointSet1d(1,0,1))
{

	/*
	 * Compute inner radius of inner most ring
	 */
	double R = radius;
	int nR = numRings;
	double rI = R/nR/sqrt(M_PI);
	coreRadius = rI;

	/*
	 * Create radial spec
	 */
	PdQuickGrid::PdQPointSet1d raySpec(nR,rI,R-rI);

	/*
	 * Compute dTheta -- ring secgment size -- attempts to create cells that have
	 * approximately the same area
	 */
	double dTheta = 2.0 * sqrt(M_PI) / ( sqrt(M_PI) * nR + 1);

	/*
	 * Create ring spec
	 */
	int numRays = (int)(2 * M_PI / dTheta) + 1;
	PdQuickGrid::PdQPointSet1d ringSpec(numRays, 0.0, 2.0*M_PI);

	/*
	 * Now get approximate cell size along length of cylinder
	 */
	double cellSize = (R-rI)/nR;

	/*
	 * Now create discretization along axis of cylinder
	 */
	int numCellsAxis = (int)(cylinderLength/cellSize)+1;
	PdQuickGrid::PdQPointSet1d axisSpec(numCellsAxis,zStart,cylinderLength);

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
	cellsPerProc = PdQuickGridMeshGenerationIterator::getNumCellsPerProcessor(globalNumberOfCells,nProcs);

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

PdGridData TensorProductSolidCylinder::allocatePdGridData() const {

	// Number of cells on this processor -- number of cells on last
	//   processor is always <= numCellsPerProc
	int numCellsAllocate = cellsPerProc[getNumProcs()-1];

	int dimension = 3;
	return PdQuickGrid::allocatePdGridData(numCellsAllocate, dimension);

}

std::pair<Cell3D,PdGridData> TensorProductSolidCylinder::computePdGridData(int proc, Cell3D cellLocator, PdGridData& pdGridData, NormFunctionPointer norm) const {

	/**
	 * This routine attempts to make as close of an analogy to CellsPerProcessor3D::computePdGridData(...)
	 * as makes sense:
	 *    x <--> r
	 *    y <--> theta
	 *    z <--> z
	 */
//	std::cout << "CellsPerProcessorCylinder::computePdGridData proc = " << proc << std::endl;
	PdQPointSet1d xSpec = specs[0]; // r
	PdQPointSet1d ySpec = specs[1]; // theta
	PdQPointSet1d zSpec = specs[2]; // z
	int nx = xSpec.getNumCells();
	int ny = ySpec.getNumCells();
	int nz = zSpec.getNumCells();

	// Each cell has the same volume arc size in radians --
	//   Cells have different volumes bases upon "r"
	double dr = xSpec.getCellSize();
	double dz = zSpec.getCellSize();
	double cellRads = ySpec.getCellSize();
	double cellVolume;

	// Number of cells on this processor
	int myNumCells = cellsPerProc[proc];


	// Coordinates used for tensor product space
	const double*x=rPtr.get();
	const double*y=thetaPtr.get();
	const double*z=zPtr.get();

	// Discretization on this processor
	pdGridData.dimension = 3;
	int numRings = nx;
	int numRays = ny;
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
	int kStart = cellLocator.k;
	int jStart = cellLocator.j;
	int iStart = cellLocator.i;

	int cell = 0;
	double point[3]={0.0,0.0,0.0};
	for(int k=kStart;k<nz && cell<myNumCells;k++){
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


		for(int j=jStart;j<ny && cell<myNumCells;j++){
			jStart=0; // begin at jStart only the first time
			cellLocator.j = j;
			double theta = y[j];
			double c = cos(theta);
			double s = sin(theta);
			for(int i=iStart;i<nx && cell<myNumCells;i++,cell++,cellLocator.i++){
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
	shared_ptr<int> neighborhoodList(new int[pdGridData.sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	pdGridData.neighborhood = neighborhoodList;
	int *neighborhood = neighborhoodList.get();
	/*
	 * number of neighbors for every point is zero
	 */
	*neighborhood = 0;
	return std::pair<Cell3D,PdGridData>(cellLocator,pdGridData);
}

}
