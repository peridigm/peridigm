/*
 * utPdQ2CylinderRingHorizon.cxx
 *
 *  Created on: Jan 6, 2010
 *      Author: jamitch
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "PdVTK.h"
#include "vtkXMLStructuredGridWriter.h"
#include "vtkXMLStructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "PdQuickGrid.h"
#include <tr1/memory>
#include <valarray>
#include <iostream>
#include <cmath>


using namespace PdQuickGrid;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

int numNeighbors(int kZ);
void ringHorizon()
{

	/*
	 * Construct ring spec
	 */
	// set ring at point = (0,0,0)
	// Note that ring is only translated to the location
	std::valarray<double> center(0.0,3);
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 0.020;
	double outerRadius = 0.025;
	double a = .050;
	double cylinderLength = 2.0*a;
	int numRings = 3;
	PdQuickGrid::PdQRing2d ring2dSpec(center,innerRadius,outerRadius,numRings);
	//	vtkSmartPointer<vtkUnstructuredGrid> g4 = PdVTK::getGrid(ring2dSpec);
	//	PdVTK::WriteUnstructured("ring2d.vtu",g4);
	/*
	 * This produces 85 cells in a ring for the Q2 cylinder geometry
	 */
	PdQPointSet1d thetaSpec(ring2dSpec.getNumRays(), 0.0, 2.0*M_PI);
	/*
	 * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1
	 *
	 */
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	int numCellsAxis = (int)(cylinderLength/cellSize) + 1;
	PdQPointSet1d axisSpec(numCellsAxis,0.0,cylinderLength);

	double SCALE=2.51;
	double horizon = SCALE*cellSize;
	RingHorizon ringHorizon = thetaSpec.getRingCellHorizon(horizon,innerRadius);

	/*
	 * Testing
	 */
	/*
	 * Check cells in ring horizon
	 */
	{
		int cell=80;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(77 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(78 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(79 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(80 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(81 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(82 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=81;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(78 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(79 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(80 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(81 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(82 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=82;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(79 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(80 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(81 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(82 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=83;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(80 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(81 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(82 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=84;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(81 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(82 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(2 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=0;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(82 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(2 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(3 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}
	{
		int cell=1;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(83 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(2 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(3 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(4 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=2;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(84 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(2 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(3 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(4 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(5 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=3;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(0 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(2 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(3 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(4 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(5 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(6 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}

	{
		int cell=4;
		RingHorizon::RingHorizonIterator hIter = ringHorizon.horizonIterator(cell);
		BOOST_CHECK(7 == hIter.numCells());
		BOOST_CHECK(hIter.hasNextCell());
		// loop over cells in horizon
		BOOST_CHECK(1 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(2 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(3 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(4 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(5 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(6 == hIter.nextCell());
		BOOST_CHECK(hIter.hasNextCell());
		BOOST_CHECK(7 == hIter.nextCell());
		BOOST_CHECK(!hIter.hasNextCell());
	}


}

void Q2CylinderNeighborhoodSizes()
{

	/*
	 * Construct ring spec
	 */
	// set ring at point = (0,0,0)
	// Note that ring is only translated to the location
	std::valarray<double> center(0.0,3);
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 0.020;
	double outerRadius = 0.025;
	double a = .050;
	double cylinderLength = 2.0*a;
	int numRings = 3;
	PdQuickGrid::PdQRing2d ring2dSpec(center,innerRadius,outerRadius,numRings);
	int numRays = ring2dSpec.getNumRays();
	//	vtkSmartPointer<vtkUnstructuredGrid> g4 = PdVTK::getGrid(ring2dSpec);
	//	PdVTK::WriteUnstructured("ring2d.vtu",g4);
	/*
	 * This produces 85 cells in a ring for the Q2 cylinder geometry
	 */
	PdQPointSet1d thetaSpec(ring2dSpec.getNumRays(), 0.0, 2.0*M_PI);
	/*
	 * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1
	 */
	double SCALE=2.51;
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	int numCellsAxis = (int)(cylinderLength/cellSize)+1;
	PdQPointSet1d axisSpec(numCellsAxis,0.0,cylinderLength);
	double horizon = SCALE*cellSize;

	// Create decomposition iterator
	int numProcs = 1;
	TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec);

	/*
	 * Testing
	 */
	{
		BOOST_CHECK(15300==cellPerProcIter.getNumGlobalCells());
		BOOST_CHECK(85==numRays);
		BOOST_CHECK(60==numCellsAxis);
		int i=0,j=0,k=0;
		int proc=0;
		Cell3D cellLocator(i,j,k);
		/*
		 * The number below can be hand calculated
		 * k=0 -->83*255
		 * k=1 -->(83+21)*255
		 * k=2 -->(83+21+21)*255
		 * k>=3 -->(83+21+21+21)*255
		 * ...
		 * k=57 -->(83+21+21)*255
		 * k=58 -->(83+21)*255
		 * k=59 -->83*255
		 *
		 * sum all of the above, then add numPoints for storing length of list at each point
		 */
		BOOST_CHECK((255*8508+15300)==cellPerProcIter.getSizeNeighborList(proc,cellLocator));
	}
	/*
	 * On one processor, its possible to calculate the length of the neighborhood list by hand
	 * size = 83
	 */
	/*
	 * Check neighborhoods
	 */
	{
		/*
		 * Every point for k=0 should have 83 neighbors
		 */
		int k=0;
		for(int j=0;j<numRays;j++){
			for(int i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				BOOST_CHECK(83==cellPerProcIter.computeNumNeighbors(i,j,k));
				BOOST_CHECK(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k=1 should have 83+21 neighbors
		 */
		int k=1;
		for(int j=0;j<numRays;j++){
			for(int i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				BOOST_CHECK((83+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
				BOOST_CHECK(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k=2 should have 83+21+21 neighbors
		 */
		int k=2;
		for(int j=0;j<numRays;j++){
			for(int i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				BOOST_CHECK((83+21+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
				BOOST_CHECK(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k=3 should have 83+21+21+21 neighbors
		 */
		int k=3;
		for(int j=0;j<numRays;j++){
			for(int i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				BOOST_CHECK((83+21+21+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
				BOOST_CHECK(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k>=3 and k<?? should have 83+21+21+21 neighbors
		 */
		int k=4;
		for(int j=0;j<numRays;j++){
			for(int i=0;i<numRings;i++){
				Cell3D cellLocator(i,j,k);
				BOOST_CHECK((83+21+21+21)==cellPerProcIter.computeNumNeighbors(i,j,k));
				BOOST_CHECK(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
			}
		}
	}

	{
		/*
		 * Every point for k<=60 up to three planes from end
		 */
		for(int k=59;k>55;k--){
			for(int j=0;j<numRays;j++){
				for(int i=0;i<numRings;i++){
					Cell3D cellLocator(i,j,k);
					BOOST_CHECK((83+(59-k)*21)==cellPerProcIter.computeNumNeighbors(i,j,k));
					BOOST_CHECK(numNeighbors(k)==cellPerProcIter.computeNumNeighbors(i,j,k));
				}
			}
		}
	}

}

void Q2CylinderNeighborhoods()
{

	/*
	 * Construct ring spec
	 */
	// set ring at point = (0,0,0)
	// Note that ring is only translated to the location
	std::valarray<double> center(0.0,3);
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 0.020;
	double outerRadius = 0.025;
	double a = .050;
	double cylinderLength = 2.0*a;
	int numRings = 3;
	PdQuickGrid::PdQRing2d ring2dSpec(center,innerRadius,outerRadius,numRings);
	//	vtkSmartPointer<vtkUnstructuredGrid> g4 = PdVTK::getGrid(ring2dSpec);
	//	PdVTK::WriteUnstructured("ring2d.vtu",g4);
	/*
	 * This produces 85 cells in a ring for the Q2 cylinder geometry
	 */
	PdQPointSet1d thetaSpec(ring2dSpec.getNumRays(), 0.0, 2.0*M_PI);
	/*
	 * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1
	 */
	double SCALE=2.51;
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	int numCellsAxis = (int)(cylinderLength/cellSize)+1;
	PdQPointSet1d axisSpec(numCellsAxis,0.0,cylinderLength);
	double dz = axisSpec.getCellSize();
	double cellRads = thetaSpec.getCellSize();
	double dr = cellSize;
	double horizon = SCALE*cellSize;

	std::tr1::shared_ptr<double> rPtr = getDiscretization(ring2dSpec.getRaySpec());
	std::tr1::shared_ptr<double> thetaPtr = getDiscretization(ring2dSpec.getRingSpec());
	std::tr1::shared_ptr<double> zPtr = getDiscretization(axisSpec);

	// Create decomposition iterator
	int numProcs = 1;
	TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec);
	PdGridData gridData;
	/*
	 * Testing
	 */
	double cylinderVolume = M_PI*(outerRadius*outerRadius-innerRadius*innerRadius)*cylinderLength;
	double sumCellVol=0.0;
	{
		int i=0,j=0,k=0,proc=0;
		Cell3D cellLocator(i,j,k);
		PdGridData pdGridDataProc0 = cellPerProcIter.allocatePdGridData();
		std::pair<Cell3D,PdGridData> p0Data = cellPerProcIter.computePdGridData(proc,cellLocator,pdGridDataProc0);
		gridData = p0Data.second;
		Cell3D nextCellLocator = p0Data.first;

		BOOST_CHECK(3==gridData.dimension);
		BOOST_CHECK(15300==gridData.globalNumPoints);
		BOOST_CHECK(15300==gridData.numPoints);
		BOOST_CHECK((255*8508+15300)==gridData.sizeNeighborhoodList);
		BOOST_CHECK(0==gridData.numExport);
		int *gIds = gridData.myGlobalIDs.get();
		for(int p=0;p<gridData.numPoints;p++,gIds++)
			BOOST_CHECK(p==*gIds);


		/*
		 * Generate a random number of points and assert their coordinates and neighborhoods
		 */
		double *X = gridData.myX.get();
		double *r = rPtr.get();
		double *theta = thetaPtr.get();
		double *z = zPtr.get();
		int nx = ring2dSpec.getNumRings();
		int ny = thetaSpec.getNumCells();
		int *neighborhoodPtr = gridData.neighborhoodPtr.get();
		int *neighborhood = gridData.neighborhood.get();
		double *vol = gridData.cellVolume.get();

		for(int k=0;k<axisSpec.getNumCells();k++){
			int ranZid = k; // rand() % axisSpec.getNumCells();
			double ranZ = z[ranZid];
			for(int j=0;j<thetaSpec.getNumCells();j++){
				int ranThetaid = j; // rand() % thetaSpec.getNumCells();
				double ranTheta = theta[ranThetaid];
				for(int i=0;i<ring2dSpec.getNumRings();i++){
					int ranRid = i; // rand() % ring2dSpec.getNumRings();
					double ranR = r[ranRid];
					double x = ranR*cos(ranTheta);
					double y = ranR*sin(ranTheta);
					double z = ranZ;
					int gId =  ranRid + ranThetaid * nx + ranZid * nx * ny;
					BOOST_CHECK(x==X[3*gId]);
					BOOST_CHECK(y==X[3*gId+1]);
					BOOST_CHECK(z==X[3*gId+2]);
					int ptr = neighborhoodPtr[gId];
					BOOST_CHECK(numNeighbors(ranZid)==neighborhood[ptr]);

					/*
					 * Volume
					 */
					double v = ranR*dr*cellRads*dz;
					BOOST_CHECK(v==vol[gId]);
					sumCellVol+= vol[gId];
				}
			}
		}
		const double tolerance = 1.0e-10;
		BOOST_CHECK_CLOSE(sumCellVol,cylinderVolume,tolerance);
	//	std::cout << "cylinderVolume = " << cylinderVolume << "; sumCellVol = " << sumCellVol << std::endl;
	}

}


/*
 * For the Q2 cylinder, and mesh discretization used in this file, this
 * function returns the number of neighbors for a point in slab kZ along the axis.
 */
int numNeighbors(int kZ){
	BOOST_CHECK(kZ<60);
	if(kZ<3)
		return 83+kZ*21;
	else if(kZ>=3 && kZ <=55)
		return 83+3*21;
	else
		return 83+(59-kZ)*21;
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPdQ2CylinderRingHorizon" );
	proc->add(BOOST_TEST_CASE( &ringHorizon ));
	proc->add(BOOST_TEST_CASE( &Q2CylinderNeighborhoodSizes ));
	proc->add(BOOST_TEST_CASE( &Q2CylinderNeighborhoods ));
	framework::master_test_suite().add( proc );

	return success;

}

bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(
		int argc,
		char* argv[]
)
{

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}

