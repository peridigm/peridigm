/*
 * utZoltanQuery_packPointsMultiFunction.cxx
 *
 *  Created on: Dec 16, 2009
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "Array.h"
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
	int numCellsZ=2;
	QUICKGRID::Spec1D zSpec(numCellsZ,xStart,xLength);

	// This scale factor pushes the horizon just over the line so that 1 cell is included
	//  in the half neighborhood
	double SCALE=1.0;
	double horizon = SCALE*spec.getCellSize();
	size_t numProc=1;

	QUICKGRID::TensorProduct3DMeshGenerator cellIter(numProc,horizon,spec,spec,zSpec);
	QUICKGRID::QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	std::pair<QUICKGRID::Cell3D,QUICKGRID::QuickGridData> pdGridData = cellIter.beginIterateProcs(pdGridDataProc0);
	QUICKGRID::QuickGridData gridData = pdGridData.second;

	return gridData;


}

void utZoltanQuery_packPointsMultiFunction()
{
	/* This test is about exercising the zoltan call back
	 * function:
void zoltanQuery_packPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numExport,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int *dest,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
)
	 *
	 */
	QUICKGRID::QuickGridData gridData = setUp();
	int numGids = 1;
	int numLids = 1;
	int numExport = 9;
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

	// the correctness of the following call has already been asserted in another unit test "utZoltanQuery_pointSizeInBytes"
	PDNEIGH::zoltanQuery_pointSizeInBytes(&gridData,numGids,numLids,numExport,exportGlobalIds,exportLocalIds,sizes,ierr);

	// sum over sizes to get number of bytes in buffer
	int totalBytes=0;
	for(int n=0;n<numExport;n++)
		totalBytes+=sizes[n];

	// Allocate buffer
	UTILITIES::Array<char> buff(totalBytes);
	char *buffPtr = buff.get();

	// create address array idx
	UTILITIES::Array<int> idx(numExport);
	int *idxPtr = idx.get();
	for(int n=0, p=0;n<numExport;n++){
		idxPtr[n] = p;
		p+= sizes[n];
	}

	int *dest = NULL;
	int zoltanErr = 0;
	/*
	 * Note that destination processors are not currently used
	 * In part, this is due to the  way we handle the load balancing step results within "Peridigm"
	 */
	// Now we wish to assert the packing up of bytes
	PDNEIGH::zoltanQuery_packPointsMultiFunction((void*)(&gridData),numGids,numLids,numExport,exportGlobalIds,exportLocalIds,dest,sizes,idxPtr,buffPtr,&zoltanErr);
	double *X = gridData.myX.get();
	double *V = gridData.cellVolume.get();
	int *neighPtr = gridData.neighborhoodPtr.get();
	int *neigh = gridData.neighborhood.get();

	// assert exportFlags on grid data
	char *exportFlag = gridData.exportFlag.get();
	// local ids 0 thru 8 should not have flags set
	for(int n=0;n<9;n++)
		BOOST_CHECK(0 == exportFlag[n]);

	// local ids 9 thru 17 should be for export
	for(int n=9;n<18;n++)
		BOOST_CHECK(1 == exportFlag[n]);

	int numNeighbors[]={ 7,11,7,11,17,11,7,11,7};
	// Now look at the buffer
	char *tmp = buffPtr;
	for(int n=0;n<numExport;n++){
		tmp=&buffPtr[idxPtr[n]];
		int id = exportLocalIds[n];
//		std::cout << "id = " << id << "; numBytes = " << sizes[n] << "; idx = " << idxPtr[n] << std::endl;

		// extract coordinates from buffer
		const double tolerance = 1.0e-15;
		int dimension = 3;
		int numBytes = dimension * sizeof(double);
		double x[] = {0.0,0.0,0.0};
		memcpy((void*)x,(void*)tmp,numBytes);
		for(int d=0;d<dimension;d++){
			BOOST_CHECK_CLOSE(x[d],X[dimension*id+d],tolerance);
		}

		// extract volume
		tmp+=numBytes;
		numBytes = sizeof(double);
		double v = 0;
		memcpy((void*)(&v),(void*)tmp,numBytes);
		BOOST_CHECK_CLOSE(v,V[id],tolerance);

		// extract neighborhood
		tmp+=numBytes;
		numBytes = sizeof(int);
		int *num; num = (int*)tmp;
		int numNeigh = *num;
		BOOST_CHECK(numNeigh == numNeighbors[n]);
		tmp+=numBytes;
		numBytes = numNeigh*sizeof(int);
		UTILITIES::Array<int> list(numNeigh);
		int *listPtr = list.get();
		memcpy((void*)(listPtr),(void*)tmp,numBytes);
		for(int j=0;j<numNeigh;j++){
			BOOST_CHECK(listPtr[j] == neigh[neighPtr[id]+1+j]);
		}

	}

	std::cout << std::endl;
}

void utZoltanQuery_unPackPointsMultiFunction_1()
{
	/* This test is about exercising the zoltan call back
	 * function:
void zoltanQuery_unPackPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numImport,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
)
	 *
	 */
	QUICKGRID::QuickGridData gridData = setUp();
	int numGids = 1;
	int numLids = 1;
	int numExport = 9;
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

	// the correctness of the following call has already been asserted in another unit test "utZoltanQuery_pointSizeInBytes"
	PDNEIGH::zoltanQuery_pointSizeInBytes(&gridData,numGids,numLids,numExport,exportGlobalIds,exportLocalIds,sizes,ierr);

	// sum over sizes to get number of bytes in buffer
	int totalBytes=0;
	for(int n=0;n<numExport;n++)
		totalBytes+=sizes[n];

	// Allocate buffer
	UTILITIES::Array<char> buff(totalBytes);
	char *buffPtr = buff.get();

	// create address array idx
	UTILITIES::Array<int> idx(numExport);
	int *idxPtr = idx.get();
	for(int n=0, p=0;n<numExport;n++){
		idxPtr[n] = p;
		p+= sizes[n];
	}

	int *dest = NULL;
	int zoltanErr = 0;
	/*
	 * Note that destination processors are not currently used
	 * In part, this is due to the  way we handle the load balancing step results within "Peridigm"
	 */
	// Correctness of the following function was asserted in the test above
	PDNEIGH::zoltanQuery_packPointsMultiFunction((void*)(&gridData),numGids,numLids,numExport,exportGlobalIds,exportLocalIds,dest,sizes,idxPtr,buffPtr,&zoltanErr);

	// Now, we need to call the unpacking function -- unpacking transforms "gridData" to import nodes
	// Lets now say that we are not importing any nodes; in this case we should get back global ids 0 thru 8
	int numImport=0;
	PDNEIGH::zoltanQuery_unPackPointsMultiFunction((void*)(&gridData),numGids,numImport,NULL,NULL,NULL,NULL,&zoltanErr);

	// By importing no points, this should have reduced the original gridData from 18 points to 9 points since they were marked for export in the pack routine
	// Assert entire data structure
	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(18 == gridData.globalNumPoints);
	BOOST_CHECK(9 == gridData.numPoints);
	BOOST_CHECK(98 == gridData.sizeNeighborhoodList);
	BOOST_CHECK(0 == gridData.numExport);
	int *ids = gridData.myGlobalIDs.get();
	for(size_t i=0;i<gridData.numPoints;i++)
		BOOST_CHECK((int)i == ids[i]);

	// Need to get a new grid data to compare it with the one coming back from unPack
	QUICKGRID::QuickGridData newGridData = setUp();
	double *newX = newGridData.myX.get(); // this has 18 points
	double *x = gridData.myX.get(); // this has 9 points -- these 9 points should be the same as the first 9 points in 'newX'
	double *newV = newGridData.cellVolume.get();
	double *v = gridData.cellVolume.get();

	int dimension=3;
	const double tolerance = 1.0e-15;
	// Note that we are only asserting the first 9 points (these
	for(size_t i=0;i<gridData.numPoints;i++){
		// Volume
		BOOST_CHECK_CLOSE(v[i],newV[i],tolerance);
		// coordinates
		for(int j=0;j<dimension;j++)
			BOOST_CHECK_CLOSE(x[i*dimension+j],newX[i*dimension+j],tolerance);
	}

	std::cout << std::endl;
}

void utZoltanQuery_unPackPointsMultiFunction_2()
{
	/* This test is about exercising the zoltan call back
	 * function:
void zoltanQuery_unPackPointsMultiFunction
(
		void *pdGridData,
		int numGids,
		int numImport,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		int *sizes,
		int *idx,
		char *buf,
		int *ierr
)
	 *
	 */
	QUICKGRID::QuickGridData gridData = setUp();
	int numGids = 1;
	int numLids = 1;
	int numExport = 9;
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

	// the correctness of the following call has already been asserted in another unit test "utZoltanQuery_pointSizeInBytes"
	PDNEIGH::zoltanQuery_pointSizeInBytes(&gridData,numGids,numLids,numExport,exportGlobalIds,exportLocalIds,sizes,ierr);

	// sum over sizes to get number of bytes in buffer
	int totalBytes=0;
	for(int n=0;n<numExport;n++)
		totalBytes+=sizes[n];

	// Allocate buffer
	UTILITIES::Array<char> buff(totalBytes);
	char *buffPtr = buff.get();

	// create address array idx
	UTILITIES::Array<int> idx(numExport);
	int *idxPtr = idx.get();
	for(int n=0, p=0;n<numExport;n++){
		idxPtr[n] = p;
		p+= sizes[n];
	}

	int *dest = NULL;
	int zoltanErr = 0;

	/*
	 * Note that destination processors are not currently used
	 * In part, this is due to the  way we handle the load balancing step results within "Peridigm"
	 */
	// Correctness of the following two function were asserted in the tests above
	PDNEIGH::zoltanQuery_packPointsMultiFunction((void*)(&gridData),numGids,numLids,numExport,exportGlobalIds,exportLocalIds,dest,sizes,idxPtr,buffPtr,&zoltanErr);
	int numImport=0;
	PDNEIGH::zoltanQuery_unPackPointsMultiFunction((void*)(&gridData),numGids,numImport,NULL,NULL,NULL,NULL,&zoltanErr);

	// At this point, grid data only the first 9 points; This was asserted in the above tests but for fun assert metadata here too!
	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(18 == gridData.globalNumPoints);
	BOOST_CHECK(9 == gridData.numPoints);
	BOOST_CHECK(98 == gridData.sizeNeighborhoodList);
	BOOST_CHECK(0 == gridData.numExport);
	int *ids = gridData.myGlobalIDs.get();
	for(size_t i=0;i<gridData.numPoints;i++)
		BOOST_CHECK((int)i == ids[i]);

	// At this point, grid data only the first 9 points; The following function should restore it to its original state
	numImport=9;
	ZOLTAN_ID_PTR importGlobalIds = exportLocalIds;
	PDNEIGH::zoltanQuery_unPackPointsMultiFunction((void*)(&gridData),numGids,numImport,importGlobalIds,sizes,idxPtr,buffPtr,&zoltanErr);

	// Now create an original PdGridData and compare with "gridData" which should be identical
	QUICKGRID::QuickGridData originalGridData = setUp();
	BOOST_CHECK(originalGridData.dimension == gridData.dimension);
	BOOST_CHECK(originalGridData.globalNumPoints == gridData.globalNumPoints);
	BOOST_CHECK(originalGridData.numPoints == gridData.numPoints);
	BOOST_CHECK(originalGridData.sizeNeighborhoodList == gridData.sizeNeighborhoodList);
	BOOST_CHECK(originalGridData.numExport == gridData.numExport);
	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(18 == gridData.globalNumPoints);
	BOOST_CHECK(18 == gridData.numPoints);
	BOOST_CHECK(196 == gridData.sizeNeighborhoodList);
	BOOST_CHECK(0 == gridData.numExport);

	// assert all quantities
	double *oX = originalGridData.myX.get();
	double *x = gridData.myX.get();
	double *oV = originalGridData.cellVolume.get();
	double *v = gridData.cellVolume.get();
	int *neighPtr = gridData.neighborhoodPtr.get();
	int *oNeighPtr = originalGridData.neighborhoodPtr.get();


	int numPoints = gridData.numPoints;
	const double tolerance = 1.0e-15;
	for(int p=0;p<numPoints;p++){
		// coordinates
		for(int d=0;d<3;d++){
			BOOST_CHECK_CLOSE(oX[p*3+d],x[p*3+d],tolerance);
		}
		// volume
		BOOST_CHECK_CLOSE(oV[p],v[p],tolerance);

		// neighborhood ptr
		BOOST_CHECK(neighPtr[p]==oNeighPtr[p]);
	}

	// assert neighborhood list
	int *neighList = gridData.neighborhood.get();
	int *oNeighList = originalGridData.neighborhood.get();
	for(int n=0;n<gridData.sizeNeighborhoodList;n++){
		BOOST_CHECK(neighList[n]==oNeighList[n]);
	}

	std::cout << std::endl;
}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utZoltanQuery_packAndUnPackPointsMultiFunction" );
	proc->add(BOOST_TEST_CASE( &utZoltanQuery_packPointsMultiFunction ));
	proc->add(BOOST_TEST_CASE( &utZoltanQuery_unPackPointsMultiFunction_1 ));
	proc->add(BOOST_TEST_CASE( &utZoltanQuery_unPackPointsMultiFunction_2 ));
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
