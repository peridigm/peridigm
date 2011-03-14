/*
 * ut_QuickGridHorizon.cxx
 *
 *  Created on: Mar 12, 2011
 *      Author: wow
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "../QuickGrid.h"
//#include "PdQuickGridParallel.h"
#include <iostream>


using namespace QUICKGRID;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

void cellNeighborhoodHorizon_11Cells()
{

	size_t numCells = 11;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*spec.getCellSize();
	Horizon h = spec.getCellHorizon(horizon);

	// 1d, 2d, and 3d meshes are formed by tensor product
	// A global id for the final mesh is computed by
	// giving a tuple (i,j,k) (for 3d) where
	// i: ith cell along the x-axis
	// j: jth cell along the y-axis
	// k: kth cell along the z-axis

	// Once constructed (with a distance), the horizon interface
	// then provides information about the neighborhood around the
	// ith cell

	size_t i=0;
	// Start is the index of cells to the left of "i" that are in the neighborhood
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(4 == h.numCells(i));

	i=1;
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(5 == h.numCells(i));

	i=2;
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(6 == h.numCells(i));

	i=3;
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(7 == h.numCells(i));

	i=4;
	BOOST_CHECK(1 == h.start(i));
	BOOST_CHECK(7 == h.numCells(i));

	i=5;
	BOOST_CHECK(2 == h.start(i));
	BOOST_CHECK(7 == h.numCells(i));

	i=6;
	BOOST_CHECK(3 == h.start(i));
	BOOST_CHECK(7 == h.numCells(i));

	i=7;
	BOOST_CHECK(4 == h.start(i));
	BOOST_CHECK(7 == h.numCells(i));

	i=8;
	BOOST_CHECK(5 == h.start(i));
	BOOST_CHECK(6 == h.numCells(i));

	i=9;
	BOOST_CHECK(6 == h.start(i));
	BOOST_CHECK(5 == h.numCells(i));

	i=10;
	BOOST_CHECK(7 == h.start(i));
	BOOST_CHECK(4 == h.numCells(i));

}

void cellNeighborhoodHorizon_3Cells()
{
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*spec.getCellSize();
	Horizon h = spec.getCellHorizon(horizon);

	// 1d, 2d, and 3d meshes are formed by tensor product
	// A global id for the final mesh is computed by
	// giving a tuple (i,j,k) (for 3d) where
	// i: ith cell along the x-axis
	// j: jth cell along the y-axis
	// k: kth cell along the z-axis

	// Once constructed (with a distance), the horizon interface
	// then provides information about the neighborhood around the
	// ith cell

	size_t i=0;
	// Start is the index of cells to the left of "i" that are in the neighborhood
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(3 == h.numCells(i));

	i=1;
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(3 == h.numCells(i));

	i=2;
	BOOST_CHECK(0 == h.start(i));
	BOOST_CHECK(3 == h.numCells(i));

}

void CellsPerProcessor3D_SerialTest_NumProcs_1()
{
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*spec.getCellSize();
	size_t numProc=1;
	double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();

	TensorProduct3DMeshGenerator cellIter(numProc,horizon,spec,spec,spec);
	QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	std::pair<Cell3D,QuickGridData> pdGridData = cellIter.beginIterateProcs(pdGridDataProc0);
	QuickGridData gridData = pdGridData.second;
	Cell3D nextCellLocator = pdGridData.first;

	size_t proc = cellIter.proc();
	// already moved to next proc
	BOOST_CHECK(1 == proc);

	BOOST_CHECK(27 == gridData.globalNumPoints);
	BOOST_CHECK(3 == gridData.dimension);
	int numPoints = gridData.numPoints;
	BOOST_CHECK(27 == numPoints);

	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	for(int id=0;id<27;id++,gIdsPtr++)
		BOOST_CHECK( *gIdsPtr == id );

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = 27 + 27*(27-1);
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	for(int id=0;id<27;id++){
		int numNeigh = *nPtr; nPtr++;
		BOOST_CHECK( 26 == numNeigh);
		for(int i=0;i<27;i++){
			if(i != id){
				BOOST_CHECK( i == *nPtr ); nPtr++;
			}
		}
	}

	// assert cell volumes
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+numPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}

}

void CellsPerProcessor3D_smallNeighborhoodSerialTest_NumProcs_1()
{
	// use this spec along the x and y axes
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);

	// this creates a different mesh along the z-axis
	size_t numCellsZ=2;
	Spec1D zSpec(numCellsZ,xStart,xLength);

	// This scale factor pushes the horizon just over the line so that 1 cell is included
	//  in the half neighborhood
	double SCALE=1.0;
	double horizon = SCALE*spec.getCellSize();
	size_t numProc=1;
	double cellVolume = spec.getCellSize() * spec.getCellSize() * zSpec.getCellSize();

	TensorProduct3DMeshGenerator cellIter(numProc,horizon,spec,spec,zSpec);
	QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	std::pair<Cell3D,QuickGridData> pdGridData = cellIter.beginIterateProcs(pdGridDataProc0);
	QuickGridData gridData = pdGridData.second;
	Cell3D nextCellLocator = pdGridData.first;

	size_t proc = cellIter.proc();
	// already moved to next proc
	BOOST_CHECK(1 == proc);

	BOOST_CHECK(18 == gridData.globalNumPoints);
	BOOST_CHECK(3 == gridData.dimension);
	size_t numPoints = gridData.numPoints;
	BOOST_CHECK(18 == numPoints);

	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	for(int id=0;id<18;id++,gIdsPtr++)
		BOOST_CHECK( *gIdsPtr == id );

	// assert length of neighborlist
	// sizeNeighborList = numPoints = sum(numNeighbors)
	// These are the correct answers that should be found by "CellsPerProcessor3D"
	int numNeighbors[]={7,11,7,11,17,11,7,11,7,
			            7,11,7,11,17,11,7,11,7};
	int neighborhoodAnswers[] = {
			                    1, 3, 4, 9, 10, 12, 13,
			                    0, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14,
			                    1, 4, 5, 10, 11, 13, 14,
			                    0, 1, 4, 6, 7, 9, 10, 12, 13, 15, 16,
			                    0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
			                    1, 2, 4, 7, 8, 10, 11, 13, 14, 16, 17,
			                    3, 4, 7, 12, 13, 15, 16,
			                    3, 4, 5, 6, 8, 12, 13, 14, 15, 16, 17,
			                    4, 5, 7, 13, 14, 16, 17,
			                    0, 1, 3, 4, 10, 12, 13,
			                    0, 1, 2, 3, 4, 5, 9, 11, 12, 13, 14,
			                    1, 2, 4, 5, 10, 13, 14,
			                    0, 1, 3, 4, 6, 7, 9, 10, 13, 15, 16,
			                    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17,
			                    1, 2, 4, 5, 7, 8, 10, 11, 13, 16, 17,
			                    3, 4, 6, 7, 12, 13, 16,
			                    3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 17,
			                    4, 5, 7, 8, 13, 14, 16
	};

	int sizeNeighborList = numPoints;
	for(int i=0;i<numPoints;i++)
		sizeNeighborList+=numNeighbors[i];

	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();
	int *neighPtr = gridData.neighborhoodPtr.get();

	// This iterates through all cell neighborhoods
	int p=0;
	int sum=0;
	for(int id=0;id<numPoints;id++){
		int numNeigh = *nPtr; nPtr++;
		// this asserts the pointers into the neighborhood
		BOOST_CHECK( neighPtr[id] == sum );
		sum += (1+numNeighbors[id]);
		// asserts number of neighbors
		BOOST_CHECK( numNeighbors[id] == numNeigh );
//		std::cout << "id = " << id << std::endl;
//		std::cout << "\t";
		// asserts neighborhood
		for(int i=0;i<numNeigh;i++){
//			std::cout << ", " << *nPtr;
			BOOST_CHECK( neighborhoodAnswers[p++] == *nPtr );
			nPtr++;
		}
//		std::cout<< std::endl;
	}

	// assert cell volumes
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+numPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}

}

void CellsPerProcessor3D_SerialTest_NumProcs_3()
{
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*spec.getCellSize();
	size_t numProcs=3;

	TensorProduct3DMeshGenerator cellIter(numProcs,horizon,spec,spec,spec);
	QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	QuickGridData  pdGridDataProcN = cellIter.allocatePdGridData();
	std::pair<Cell3D,QuickGridData> p0Data = cellIter.beginIterateProcs(pdGridDataProc0);

	QuickGridData gridData = p0Data.second;
	Cell3D nextCellLocator = p0Data.first;
	// proc 0
	BOOST_CHECK(27 == gridData.globalNumPoints);
	BOOST_CHECK(3 == gridData.dimension);
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(9 == myNumPoints);

	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 0;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++)
		BOOST_CHECK( *gIdsPtr == id );

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );

	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	// each cell in neighborhood has the entire list of cells in mesh except itself
	for(int id=0;id<myNumPoints;id++){
		int numNeigh = *nPtr; nPtr++;
		BOOST_CHECK( 26 == numNeigh);
		for(int i=0;i<27;i++){
			if(i != id){
				BOOST_CHECK( i == *nPtr ); nPtr++;
			}
		}
	}

	// assert cell volumes
	double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}

	// already moved to next proc
	size_t proc = 1;
	start = 9;
	while(cellIter.hasNextProc()){

		BOOST_CHECK(proc == cellIter.proc());
		std::pair<Cell3D,QuickGridData> data = cellIter.nextProc(nextCellLocator,pdGridDataProcN);

		QuickGridData gridData = data.second;
		nextCellLocator = data.first;

		BOOST_CHECK(27 == gridData.globalNumPoints);
		BOOST_CHECK(3 == gridData.dimension);
		int myNumPoints = gridData.numPoints;
		BOOST_CHECK(9 == myNumPoints);

		// assert global ids for this processor
		shared_ptr<int> gIds = gridData.myGlobalIDs;
		int *gIdsPtr = gIds.get();
		for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
			BOOST_CHECK( *gIdsPtr == id );
		}

		// there are 9 nodes per processor
		start += 9;

		// assert length of neighborlist
		// sizeNeighborList = myNumCells + myNumCells*numNeighbors
		int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
		BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );

		// assert neighbor lists; in this case each point has all points
		shared_ptr<int> neighborList = gridData.neighborhood;
		int *nPtr = neighborList.get();

		// This iterates through all cell neighborhoods
		// each cell in neighborhood has the entire list of cells in mesh except itself
		gIdsPtr = gIds.get();
		for(int id=0;id<myNumPoints;id++, gIdsPtr++){
			int numNeigh = *nPtr; nPtr++;
			BOOST_CHECK(26 == numNeigh);

			for(int i=0;i<27;i++){
				if(*gIdsPtr != i) {
					BOOST_CHECK(i == *nPtr);
					nPtr++;
				}
			}

		}

		// assert cell volumes
		double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();
		const double tolerance = 1.0e-15;
		double *v = gridData.cellVolume.get();
		double *end = v+myNumPoints;
		for(; v != end ; v++){
			BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
		}

		proc++;

	}

}

void CellsPerProcessor3D_SerialTest_NumProcs_4()
{
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*spec.getCellSize();
	size_t numProcs=4;

	TensorProduct3DMeshGenerator cellIter(numProcs,horizon,spec,spec,spec);
	QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	QuickGridData  pdGridDataProcN = cellIter.allocatePdGridData();
	std::pair<Cell3D,QuickGridData> p0Data = cellIter.beginIterateProcs(pdGridDataProc0);

	QuickGridData gridData = p0Data.second;
	Cell3D nextCellLocator = p0Data.first;
	// proc 0
	BOOST_CHECK(27 == gridData.globalNumPoints);
	BOOST_CHECK(3 == gridData.dimension);
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(6 == myNumPoints);

	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 0;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++)
		BOOST_CHECK( *gIdsPtr == id );

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );

	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	// each cell in neighborhood has the entire list of cells in mesh except itself
	for(int id=0;id<myNumPoints;id++){
		int numNeigh = *nPtr; nPtr++;
		BOOST_CHECK( 26 == numNeigh);
		for(int i=0;i<27;i++){
			if(i != id){
				BOOST_CHECK( i == *nPtr ); nPtr++;
			}
		}
	}

	// assert cell volumes
	double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}


	// already moved to next proc
	size_t proc = 1;
	start = 6;
	while(cellIter.hasNextProc()){

		BOOST_CHECK(proc == cellIter.proc());
		std::pair<Cell3D,QuickGridData> data = cellIter.nextProc(nextCellLocator,pdGridDataProcN);

		QuickGridData gridData = data.second;
		nextCellLocator = data.first;

		BOOST_CHECK(27 == gridData.globalNumPoints);
		BOOST_CHECK(3 == gridData.dimension);

		int myNumPoints = gridData.numPoints;
		int answerMyNumPoints=6;
		if(proc==numProcs-1)
			answerMyNumPoints = 9;

		BOOST_CHECK(answerMyNumPoints == myNumPoints);
		// assert global ids for this processor
		shared_ptr<int> gIds = gridData.myGlobalIDs;
		int *gIdsPtr = gIds.get();
		for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
			BOOST_CHECK( *gIdsPtr == id );
		}

		// there are 6 nodes per processor
		start += 6;

		// assert length of neighborlist
		// sizeNeighborList = myNumCells + myNumCells*numNeighbors
		int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
		BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );

		// assert neighbor lists; in this case each point has all points
		shared_ptr<int> neighborList = gridData.neighborhood;
		int *nPtr = neighborList.get();

		// This iterates through all cell neighborhoods
		// each cell in neighborhood has the entire list of cells in mesh except itself
		gIdsPtr = gIds.get();
		for(int id=0;id<myNumPoints;id++, gIdsPtr++){
			int numNeigh = *nPtr; nPtr++;
			BOOST_CHECK(26 == numNeigh);

			for(int i=0;i<27;i++){
				if(*gIdsPtr != i) {
					BOOST_CHECK(i == *nPtr);
					nPtr++;
				}
			}

		}

		// assert cell volumes
		double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();
		const double tolerance = 1.0e-15;
		double *v = gridData.cellVolume.get();
		double *end = v+myNumPoints;
		for(; v != end ; v++){
			BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
		}

		proc++;

	}
}


void CellsPerProcessor3D_SerialTest_NumProcs_5 ()
{
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	Spec1D spec(numCells,xStart,xLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*spec.getCellSize();
	size_t numProcs=5;

	TensorProduct3DMeshGenerator cellIter(numProcs,horizon,spec,spec,spec);
	QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	QuickGridData  pdGridDataProcN = cellIter.allocatePdGridData();
	std::pair<Cell3D,QuickGridData> p0Data = cellIter.beginIterateProcs(pdGridDataProc0);

	QuickGridData gridData = p0Data.second;
	Cell3D nextCellLocator = p0Data.first;
	// proc 0
	BOOST_CHECK(27 == gridData.globalNumPoints);
	BOOST_CHECK(3 == gridData.dimension);

	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(5 == myNumPoints);

	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 0;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++)
		BOOST_CHECK( *gIdsPtr == id );

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );

	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	// each cell in neighborhood has the entire list of cells in mesh except itself
	for(int id=0;id<myNumPoints;id++){
		int numNeigh = *nPtr; nPtr++;
		BOOST_CHECK( 26 == numNeigh);
		for(int i=0;i<27;i++){
			if(i != id){
				BOOST_CHECK( i == *nPtr ); nPtr++;
			}
		}
	}

	// assert cell volumes
	double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}

	// already moved to next proc
	size_t proc = 1;
	start = 5;
	while(cellIter.hasNextProc()){

		BOOST_CHECK(proc == cellIter.proc());
		std::pair<Cell3D,QuickGridData> data = cellIter.nextProc(nextCellLocator,pdGridDataProcN);

		QuickGridData gridData = data.second;
		nextCellLocator = data.first;

		BOOST_CHECK(27 == gridData.globalNumPoints);
		BOOST_CHECK(3 == gridData.dimension);

		int myNumPoints = gridData.numPoints;
		int answerMyNumPoints=5;
		if(proc==numProcs-1)
			answerMyNumPoints = 7;

		BOOST_CHECK(answerMyNumPoints == myNumPoints);
		// assert global ids for this processor
		shared_ptr<int> gIds = gridData.myGlobalIDs;
		int *gIdsPtr = gIds.get();
		for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
			BOOST_CHECK( *gIdsPtr == id );
		}

		// there are 5 points per processor
		start += 5;

		// assert length of neighborlist
		// sizeNeighborList = myNumCells + myNumCells*numNeighbors
		int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
		BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );

		// assert neighbor lists; in this case each point has all points
		shared_ptr<int> neighborList = gridData.neighborhood;
		int *nPtr = neighborList.get();

		// This iterates through all cell neighborhoods
		// each cell in neighborhood has the entire list of cells in mesh except itself
		gIdsPtr = gIds.get();
		for(int id=0;id<myNumPoints;id++, gIdsPtr++){
			int numNeigh = *nPtr; nPtr++;
			BOOST_CHECK(26 == numNeigh);

			for(int i=0;i<27;i++){
				if(*gIdsPtr != i) {
					BOOST_CHECK(i == *nPtr);
					nPtr++;
				}
			}

		}

		// assert cell volumes
		double cellVolume = spec.getCellSize() * spec.getCellSize() * spec.getCellSize();
		const double tolerance = 1.0e-15;
		double *v = gridData.cellVolume.get();
		double *end = v+myNumPoints;
		for(; v != end ; v++){
			BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
		}

		proc++;

	}
}

void CellsPerProcessor3D_Large()
{
	size_t numCellsX = 10;
	double xStart = 0.0;
	double xLength=1.0;
	Spec1D specX(numCellsX,xStart,xLength);

	size_t numCellsY = 10;
	double yStart = 0.0;
	double yLength=10.0;
	Spec1D specY(numCellsY,yStart,yLength);

	size_t numCellsZ = 1000;
	double zStart = 0.0;
	double zLength=100.0;
	Spec1D specZ(numCellsZ,zStart,zLength);
	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=2.51;
	double horizon = SCALE*specX.getCellSize();
	size_t numProcs=5;

	TensorProduct3DMeshGenerator cellIter(numProcs,horizon,specX,specY,specZ);
	QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	QuickGridData  pdGridDataProcN = cellIter.allocatePdGridData();
	std::pair<Cell3D,QuickGridData> p0Data = cellIter.beginIterateProcs(pdGridDataProc0);

	QuickGridData gridData = p0Data.second;
	Cell3D nextCellLocator = p0Data.first;
	// proc 0
	// already moved to next proc
	size_t proc = 1;

	while(cellIter.hasNextProc()){

		BOOST_CHECK(proc == cellIter.proc());
		std::pair<Cell3D,QuickGridData> data = cellIter.nextProc(nextCellLocator,pdGridDataProcN);

		QuickGridData gridData = data.second;
		nextCellLocator = data.first;

		proc++;
	}
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "ut_QuickGridHorizon" );
	proc->add(BOOST_TEST_CASE( &cellNeighborhoodHorizon_11Cells ));
	proc->add(BOOST_TEST_CASE( &cellNeighborhoodHorizon_3Cells ));
	proc->add(BOOST_TEST_CASE( &CellsPerProcessor3D_SerialTest_NumProcs_1 ));
	proc->add(BOOST_TEST_CASE( &CellsPerProcessor3D_SerialTest_NumProcs_3 ));
	proc->add(BOOST_TEST_CASE( &CellsPerProcessor3D_SerialTest_NumProcs_4 ));
	proc->add(BOOST_TEST_CASE( &CellsPerProcessor3D_SerialTest_NumProcs_5 ));
	proc->add(BOOST_TEST_CASE( &CellsPerProcessor3D_Large ));
	proc->add(BOOST_TEST_CASE( &CellsPerProcessor3D_smallNeighborhoodSerialTest_NumProcs_1 ));
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

