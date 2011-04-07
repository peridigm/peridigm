/*! \file utPdSmallMeshCylinder_np4.cxx */

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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
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
const std::valarray<double> center(0.0,3);
const double ringThickness = 2.0;
const double innerRadius = ringThickness*(7/(2.0*M_PI)-1)*.5;
const double outerRadius = innerRadius+ringThickness;
const int numRings = 2;
const double cellSize=ringThickness/numRings;
const PdQuickGrid::PdQRing2d ring2dSpec(center,innerRadius,outerRadius,numRings);
const double cylinderLength = 2.0*cellSize;
const int numCellsAxis = 2;
const PdQPointSet1d axisSpec(numCellsAxis,0.0,cylinderLength);

// This is a hack to get the correct number of cells in a ring
const double SCALE=1.51;
const double horizon = SCALE*cellSize;
static PdGridData gridData;

const int neighborAnswers[] = {
		12,13,14,15,1,2,3,4,5,28,29,30,31,16,17,18,19,20,21,
		12,13,14,15,0,2,3,4,5,28,29,30,31,16,17,18,19,20,21,
		      14,15,0,1,3,4,5, 6, 7,30,31,16,17,18,19,20,21,22,23,
		      14,15,0,1,2,4,5, 6, 7,30,31,16,17,18,19,20,21,22,23,
		            0,1,2,3,5, 6, 7,8,9, 16,17,18,19,20,21,22,23,24,25,
		            0,1,2,3,4, 6, 7,8,9, 16,17,18,19,20,21,22,23,24,25,
		            2,3,4,5,7,8,9,10,11,18,19,20,21,22,23,24,25,26,27,
		            2,3,4,5,6,8,9,10,11,18,19,20,21,22,23,24,25,26,27,
		            4,5,6,7,9,10,11,12,13,20,21,22,23,24,25,26,27,28,29,
		            4,5,6,7,8,10,11,12,13,20,21,22,23,24,25,26,27,28,29,
		            6,7,8,9,11,12,13,14,15,22,23,24,25,26,27,28,29,30,31,
		            6,7,8,9,10,12,13,14,15,22,23,24,25,26,27,28,29,30,31,
		            8,9,10,11,13,14,15,0,1,24,25,26,27,28,29,30,31,16,17,
		            8,9,10,11,12,14,15,0,1,24,25,26,27,28,29,30,31,16,17,
		            10,11,12,13,15,0,1,2,3,26,27,28,29,30,31,16,17,18,19,
		            10,11,12,13,14,0,1,2,3,26,27,28,29,30,31,16,17,18,19,
		            // start new z-plane
					12,13,14,15,0,1,2,3,4,5,28,29,30,31,17,18,19,20,21,
					12,13,14,15,0,1,2,3,4,5,28,29,30,31,16,18,19,20,21,
					      14,15,0,1,2,3,4,5, 6, 7,30,31,16,17,19,20,21,22,23,
					      14,15,0,1,2,3,4,5, 6, 7,30,31,16,17,18,20,21,22,23,
					            0,1,2,3,4,5, 6, 7,8,9, 16,17,18,19,21,22,23,24,25,
					            0,1,2,3,4,5, 6, 7,8,9, 16,17,18,19,20,22,23,24,25,
					            2,3,4,5,6,7,8,9,10,11,18,19,20,21,23,24,25,26,27,
					            2,3,4,5,6,7,8,9,10,11,18,19,20,21,22,24,25,26,27,
					            4,5,6,7,8,9,10,11,12,13,20,21,22,23,25,26,27,28,29,
					            4,5,6,7,8,9,10,11,12,13,20,21,22,23,24,26,27,28,29,
					            6,7,8,9,10,11,12,13,14,15,22,23,24,25,27,28,29,30,31,
					            6,7,8,9,10,11,12,13,14,15,22,23,24,25,26,28,29,30,31,
					            8,9,10,11,12,13,14,15,0,1,24,25,26,27,29,30,31,16,17,
					            8,9,10,11,12,13,14,15,0,1,24,25,26,27,28,30,31,16,17,
					            10,11,12,13,14,15,0,1,2,3,26,27,28,29,31,16,17,18,19,
					            10,11,12,13,14,15,0,1,2,3,26,27,28,29,30,16,17,18,19

};

void p0(){

	/*
	 * This test is for proc 0 only
	 */
	BOOST_CHECK(0 == myRank);

	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(32 == gridData.globalNumPoints);
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(8 == myNumPoints);

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*19;
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	BOOST_CHECK(0 == gridData.numExport);


	double *X = gridData.myX.get();
	int *neighborhoodPtr = gridData.neighborhoodPtr.get();
	int *neighborhood = gridData.neighborhood.get();
	double *vol = gridData.cellVolume.get();
	std::tr1::shared_ptr<double> meshPtr = PdQuickGrid::getDiscretization(ring2dSpec,axisSpec);
	double dr = ring2dSpec.getRaySpec().getCellSize();
	double dz = axisSpec.getCellSize();
	double cellRads = ring2dSpec.getRingSpec().getCellSize();

	double *xx = meshPtr.get();
	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int cell = 0;
	int start=0;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++,cell++){
		BOOST_CHECK( *gIdsPtr == id );

		BOOST_CHECK(xx[id*3]==X[3*cell]);
		BOOST_CHECK(xx[id*3+1]==X[3*cell+1]);
		BOOST_CHECK(xx[id*3+2]==X[3*cell+2]);
		int ptr = neighborhoodPtr[cell];
		BOOST_CHECK(19==neighborhood[ptr]);
		for(int p=0;p<19;p++){
			BOOST_CHECK(neighborAnswers[p+id*19]==neighborhood[ptr+1+p]);
		}
		double r = sqrt(xx[id*3]*xx[id*3]+xx[id*3+1]*xx[id*3+1]);
		/*
		 * Volume
		 */
		const double tolerance = 1.0e-13;
		double v = r*dr*cellRads*dz;
		BOOST_CHECK_CLOSE(v,vol[cell],tolerance);
	}

}

void p1(){

	/*
	 * This test is for proc 1 only
	 */
	BOOST_CHECK(1 == myRank);

	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(32 == gridData.globalNumPoints);
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(8 == myNumPoints);

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*19;
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	BOOST_CHECK(0 == gridData.numExport);


	double *X = gridData.myX.get();
	int *neighborhoodPtr = gridData.neighborhoodPtr.get();
	int *neighborhood = gridData.neighborhood.get();
	double *vol = gridData.cellVolume.get();
	std::tr1::shared_ptr<double> meshPtr = PdQuickGrid::getDiscretization(ring2dSpec,axisSpec);
	double dr = ring2dSpec.getRaySpec().getCellSize();
	double dz = axisSpec.getCellSize();
	double cellRads = ring2dSpec.getRingSpec().getCellSize();

	double *xx = meshPtr.get();
	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int cell = 0;
	int start=8;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++,cell++){
		BOOST_CHECK( *gIdsPtr == id );

		BOOST_CHECK(xx[id*3]==X[3*cell]);
		BOOST_CHECK(xx[id*3+1]==X[3*cell+1]);
		BOOST_CHECK(xx[id*3+2]==X[3*cell+2]);
		int ptr = neighborhoodPtr[cell];
		BOOST_CHECK(19==neighborhood[ptr]);
		for(int p=0;p<19;p++){
			BOOST_CHECK(neighborAnswers[p+id*19]==neighborhood[ptr+1+p]);
		}
		double r = sqrt(xx[id*3]*xx[id*3]+xx[id*3+1]*xx[id*3+1]);
		/*
		 * Volume
		 */
		const double tolerance = 1.0e-13;
		double v = r*dr*cellRads*dz;
		BOOST_CHECK_CLOSE(v,vol[cell],tolerance);
	}

}

void p2(){

	/*
	 * This test is for proc 2 only
	 */
	BOOST_CHECK(2 == myRank);

	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(32 == gridData.globalNumPoints);
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(8 == myNumPoints);

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*19;
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	BOOST_CHECK(0 == gridData.numExport);


	double *X = gridData.myX.get();
	int *neighborhoodPtr = gridData.neighborhoodPtr.get();
	int *neighborhood = gridData.neighborhood.get();
	double *vol = gridData.cellVolume.get();
	std::tr1::shared_ptr<double> meshPtr = PdQuickGrid::getDiscretization(ring2dSpec,axisSpec);
	double dr = ring2dSpec.getRaySpec().getCellSize();
	double dz = axisSpec.getCellSize();
	double cellRads = ring2dSpec.getRingSpec().getCellSize();

	double *xx = meshPtr.get();
	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int cell = 0;
	int start=16;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++,cell++){
		BOOST_CHECK( *gIdsPtr == id );

		BOOST_CHECK(xx[id*3]==X[3*cell]);
		BOOST_CHECK(xx[id*3+1]==X[3*cell+1]);
		BOOST_CHECK(xx[id*3+2]==X[3*cell+2]);
		int ptr = neighborhoodPtr[cell];
		BOOST_CHECK(19==neighborhood[ptr]);
		for(int p=0;p<19;p++){
			BOOST_CHECK(neighborAnswers[p+id*19]==neighborhood[ptr+1+p]);
		}
		double r = sqrt(xx[id*3]*xx[id*3]+xx[id*3+1]*xx[id*3+1]);
		/*
		 * Volume
		 */
		const double tolerance = 1.0e-13;
		double v = r*dr*cellRads*dz;
		BOOST_CHECK_CLOSE(v,vol[cell],tolerance);
	}

}

void p3(){

	/*
	 * This test is for proc 3 only
	 */
	BOOST_CHECK(3 == myRank);

	BOOST_CHECK(3 == gridData.dimension);
	BOOST_CHECK(32 == gridData.globalNumPoints);
	int myNumPoints = gridData.numPoints;
	BOOST_CHECK(8 == myNumPoints);

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*19;
	BOOST_CHECK( sizeNeighborList == gridData.sizeNeighborhoodList );
	BOOST_CHECK(0 == gridData.numExport);


	double *X = gridData.myX.get();
	int *neighborhoodPtr = gridData.neighborhoodPtr.get();
	int *neighborhood = gridData.neighborhood.get();
	double *vol = gridData.cellVolume.get();
	std::tr1::shared_ptr<double> meshPtr = PdQuickGrid::getDiscretization(ring2dSpec,axisSpec);
	double dr = ring2dSpec.getRaySpec().getCellSize();
	double dz = axisSpec.getCellSize();
	double cellRads = ring2dSpec.getRingSpec().getCellSize();

	double *xx = meshPtr.get();
	// Assert global ids for this processor
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int cell = 0;
	int start=24;
	for(int id=start;id<gridData.numPoints+start;id++,gIdsPtr++,cell++){
		BOOST_CHECK( *gIdsPtr == id );

		BOOST_CHECK(xx[id*3]==X[3*cell]);
		BOOST_CHECK(xx[id*3+1]==X[3*cell+1]);
		BOOST_CHECK(xx[id*3+2]==X[3*cell+2]);
		int ptr = neighborhoodPtr[cell];
		BOOST_CHECK(19==neighborhood[ptr]);
		for(int p=0;p<19;p++){
			BOOST_CHECK(neighborAnswers[p+id*19]==neighborhood[ptr+1+p]);
		}
		double r = sqrt(xx[id*3]*xx[id*3]+xx[id*3+1]*xx[id*3+1]);
		/*
		 * Volume
		 */
		const double tolerance = 1.0e-13;
		double v = r*dr*cellRads*dz;
		BOOST_CHECK_CLOSE(v,vol[cell],tolerance);
	}

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	if(0 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "utPdSmallMeshCylinder_np4p0" );
		proc->add(BOOST_TEST_CASE( &p0 ));
		framework::master_test_suite().add( proc );
		return success;
	}
	if(1 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "utPdSmallMeshCylinder_np4p1" );
		proc->add(BOOST_TEST_CASE( &p1 ));
		framework::master_test_suite().add( proc );
		return success;
	}
	if(2 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "utPdSmallMeshCylinder_np4p2" );
		proc->add(BOOST_TEST_CASE( &p2 ));
		framework::master_test_suite().add( proc );
		return success;
	}

	if(3 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "utPdSmallMeshCylinder_np4p3" );
		proc->add(BOOST_TEST_CASE( &p3 ));
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
	 * This test only make sense for numProcs == 4
	 */
	if(4 != numProcs){
		std::cerr << "Unit test runtime ERROR: utPdSmallMeshCylinder_np4 only makes sense on 4 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 4 ./utPdSmallMeshCylinder_np4" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	TensorProductCylinderMeshGenerator cellIter(numProcs, horizon,ring2dSpec, axisSpec);
	gridData = getDiscretization(myRank, cellIter);

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
