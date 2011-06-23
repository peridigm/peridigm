/*
 * utZoltanQuery_pointSizeInBytes.cxx
 *
 *  Created on: Dec 4, 2009
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include <iostream>


using std::tr1::shared_ptr;
using namespace boost::unit_test;

QUICKGRID::QuickGridData setUp(){
	/*
	 * This setup is identical to the setup in utPdQuickGridHorizon.cxx
	 *     function call "CellsPerProcessor3D_smallNeighborhoodSerialTest_NumProcs_1"
	 * That function tests the setup rigorously
	 *
	 */

	// use this spec along the x and y axes
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	QUICKGRID::Spec1D spec(numCells,xStart,xLength);

	// this creates a different mesh along the z-axis
	size_t numCellsZ=2;
	QUICKGRID::Spec1D zSpec(numCellsZ,xStart,xLength);

	// This scale factor pushes the horizon just over the line so that 1 cell is included
	//  in the half neighborhood
	double SCALE=1.0;
	double horizon = SCALE*spec.getCellSize();
	int numProc=1;

	QUICKGRID::TensorProduct3DMeshGenerator cellIter(numProc,horizon,spec,spec,zSpec);
	QUICKGRID::QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	std::pair<QUICKGRID::Cell3D,QUICKGRID::QuickGridData> pdGridData = cellIter.beginIterateProcs(pdGridDataProc0);
	QUICKGRID::QuickGridData gridData = pdGridData.second;

	return gridData;


}

void zoltanQuery_pointSizeInBytes_smallNeighborhood()
{
	/* This test is about exercising the zoltan call back
	 * function:
	 * void zoltanQuery_pointSizeInBytes
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numPoints,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int *sizes,
		int *ierr
)
	 *
	 */
	QUICKGRID::QuickGridData gridData = setUp();
	int numGids = 1;
	int numLids = 1;
	int numPoints = 9;

	/*
	 * The above setup has 2 slabs of 9 nodes each
	 * Use the call function as if we wanted to ship out 1 of the slabs
	 *
	 */
	ZOLTAN_ID_TYPE exportLocalIds[] = {9, 10, 11, 12, 13, 14, 15, 16, 17};
	// Global ids are not used in the size function
	ZOLTAN_ID_TYPE *exportGlobalIds = NULL;
	int sizes[] = {0,0,0,0,0,0,0,0,0};
	int ierr[] = {0};
	PDNEIGH::zoltanQuery_pointSizeInBytes(&gridData,numGids,numLids,numPoints,exportGlobalIds,exportLocalIds,sizes,ierr);

	// assert sizes
	// these are the known number of neighbors for each point being exported
	int numNeighbors[]={ 7,11,7,11,17,11,7,11,7};
	               // coordinates + volume + numNeigh + neighbors
	int sizeAnswers[] = {0,0,0,0,0,0,0,0,0};
	int ptrAnswers[] = {0,0,0,0,0,0,0,0,0};
	int dimension=3;
	int ptrStart = 4*7 + 4 + 4*11 + 4 + 17*1 + 1; // = 98
	ptrAnswers[0] = ptrStart;
	for(int p=1;p<numPoints;p++){
		ptrAnswers[p] = ptrAnswers[p-1] + numNeighbors[p-1] + 1;
	}

	for(int p=0;p<numPoints;p++){
		sizeAnswers[p]=(dimension+1)*sizeof(double)+(1+numNeighbors[p])*sizeof(int);
	}
	int *neighPtr = gridData.neighborhoodPtr.get();
	for(int p=0;p<9;p++){
		BOOST_CHECK( ptrAnswers[p] == neighPtr[exportLocalIds[p]] );
		BOOST_CHECK( sizeAnswers[p] == sizes[p] );
	}
	std::cout << std::endl;
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utZoltanQuery_pointSizeInBytes" );
	proc->add(BOOST_TEST_CASE( &zoltanQuery_pointSizeInBytes_smallNeighborhood ));
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

