/*
 * utPdQuickGridMPI_np2_3x3x2.cxx
 *
 *  Created on: Dec 4, 2009
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdutMpiFixture.h"
#include <iostream>

using namespace PdQuickGrid;
using namespace Pdut;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;

static int myRank;
static int numProcs;
const int nx = 3;
const int ny = 3;
const int nz = 2;
const double xStart = 1.0;
const double xLength = 1.0;
const double yStart = 1.0;
const double yLength = 1.0;
const double zStart = 1.0;
const double zLength = 1.0;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
static PdGridData gridData;

void p0()
{
	BOOST_CHECK(0 == myRank);
	/*
	 * problem dimension is 3
	 */
	BOOST_CHECK(3 == gridData.dimension);

	/*
	 * Total number of cells in test
	 */
	BOOST_CHECK(nx*ny*nz == gridData.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(9 == myNumPoints);

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 0;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++)
		BOOST_CHECK( *gIdsPtr == id );

	// assert length of neighborlist
	// sizeNeighborList = numPoints = sum(numNeighbors)
	// These are the correct answers that should be found by "CellsPerProcessor3D"
	int numNeighbors[]={7,11,7,11,17,11,7,11,7};
	int neighborhoodAnswers[] = {
            1, 3, 4, 9, 10, 12, 13,
            0, 2, 3, 4, 5, 9, 10, 11, 12, 13, 14,
            1, 4, 5, 10, 11, 13, 14,
            0, 1, 4, 6, 7, 9, 10, 12, 13, 15, 16,
            0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
            1, 2, 4, 7, 8, 10, 11, 13, 14, 16, 17,
            3, 4, 7, 12, 13, 15, 16,
            3, 4, 5, 6, 8, 12, 13, 14, 15, 16, 17,
            4, 5, 7, 13, 14, 16, 17
	};

	int sizeNeighborList = myNumPoints;
	for(int i=0;i<myNumPoints;i++)
		sizeNeighborList+=numNeighbors[i];

	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();
	int *neighPtr = gridData.neighborhoodPtr.get();

	// This iterates through all cell neighborhoods
	int p=0;
	int sum=0;
	for(int id=0;id<myNumPoints;id++){
		int numNeigh = *nPtr; nPtr++;
		// this asserts the pointers into the neighborhood
		BOOST_CHECK( neighPtr[id] == sum);
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
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}

	std::cout << std::endl;


}

void p1()
{
	BOOST_CHECK(1 == myRank);
	/*
	 * problem dimension is 3
	 */
	BOOST_CHECK(3 == gridData.dimension);

	/*
	 * Total number of cells in test
	 */
	BOOST_CHECK(nx*ny*nz == gridData.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(9 == myNumPoints);

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 9;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
		BOOST_CHECK( *gIdsPtr == id );
	}

	// assert length of neighborlist
	// sizeNeighborList = numPoints = sum(numNeighbors)
	// These are the correct answers that should be found by "CellsPerProcessor3D"
	int numNeighbors[]={ 7,11,7,11,17,11,7,11,7};
	int neighborhoodAnswers[] = {
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

	int sizeNeighborList = myNumPoints;
	for(int i=0;i<myNumPoints;i++)
		sizeNeighborList+=numNeighbors[i];

	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = gridData.neighborhood;
	int *nPtr = neighborList.get();
	int *neighPtr = gridData.neighborhoodPtr.get();

	// This iterates through all cell neighborhoods
	int p=0;
	int sum=0; // starting point for ptr
	for(int id=0;id<myNumPoints;id++){
		int numNeigh = *nPtr; nPtr++;
		// this asserts the pointers into the neighborhood
		BOOST_CHECK( neighPtr[id] == sum);
		sum += (1+ numNeighbors[id]);
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
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		BOOST_CHECK_CLOSE(*v,cellVolume,tolerance);
	}

	std::cout << std::endl;

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	if(0 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "utPdQuickGridMPI_np2_3x3x2p0" );
		proc->add(BOOST_TEST_CASE( &p0 ));
		framework::master_test_suite().add( proc );
		return success;
	}
	if(1 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "utPdQuickGridMPI_np2_3x3x2p1" );
		proc->add(BOOST_TEST_CASE( &p1 ));
		framework::master_test_suite().add( proc );
		return success;
	}
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

	// Initialize MPI and timer
	PdutMpiFixture myMpi = PdutMpiFixture(argc,argv);

	// These are static (file scope) variables
	myRank = myMpi.rank;
	numProcs = myMpi.numProcs;
	/**
	 * This test only make sense for numProcs == 2
	 */
	if(2 != numProcs){
		std::cerr << "Unit test runtime ERROR: utPdQuickGridMPI_np2_3x3x2 only makes sense on 2 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 2 ./utPdQuickGridMPI_np2_3x3x2" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=1.0;
	double horizon = SCALE*xSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	gridData =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

//	gridData = PdQuickGrid::getDiscretization(myRank, numProcs, horizon, xSpec, ySpec, zSpec);

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
