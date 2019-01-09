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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "../QuickGrid.h"
#include "PdZoltan.h"
#include "PdutMpiFixture.h"
#include <iostream>

using std::shared_ptr;
using namespace Pdut;
using std::cout;

static size_t myRank;
static size_t numProcs;
const int nx = 3;
const int ny = 3;
const int nz = 2;
const double xStart = 1.0;
const double xLength = 1.0;
const double yStart = 1.0;
const double yLength = 1.0;
const double zStart = 1.0;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;
static QUICKGRID::QuickGridData gridData;

TEUCHOS_UNIT_TEST( QuickGridMPI_np2_3x3x2, p0Test) {

        if (myRank == 0) {

	TEST_ASSERT(0 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == gridData.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == gridData.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = gridData.numPoints;
	TEST_ASSERT(9 == myNumPoints);

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 0;
	for(size_t id=start;id<gridData.numPoints+start;id++,gIdsPtr++)
		TEST_ASSERT( *gIdsPtr == (int)id );

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

	TEST_ASSERT( sizeNeighborList == gridData.sizeNeighborhoodList );
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
		TEST_ASSERT( neighPtr[id] == sum);
		sum += (1+numNeighbors[id]);
		// asserts number of neighbors
		TEST_ASSERT( numNeighbors[id] == numNeigh );
//		std::cout << "id = " << id << std::endl;
//		std::cout << "\t";
		// asserts neighborhood
		for(int i=0;i<numNeigh;i++){
//			std::cout << ", " << *nPtr;
			TEST_ASSERT( neighborhoodAnswers[p++] == *nPtr );
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
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

	std::cout << std::endl;

    }


}


TEUCHOS_UNIT_TEST( QuickGridMPI_np2_3x3x2, p1Test) {


        if (myRank == 1){

	TEST_ASSERT(1 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == gridData.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == gridData.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = gridData.numPoints;
	TEST_ASSERT(9 == myNumPoints);

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 9;
	for(size_t id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
		TEST_ASSERT( *gIdsPtr == (int)id );
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

	TEST_ASSERT( sizeNeighborList == gridData.sizeNeighborhoodList );
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
		TEST_ASSERT( neighPtr[id] == sum);
		sum += (1+ numNeighbors[id]);
		// asserts number of neighbors
		TEST_ASSERT( numNeighbors[id] == numNeigh );
//		std::cout << "id = " << id << std::endl;
//		std::cout << "\t";
		// asserts neighborhood
		for(int i=0;i<numNeigh;i++){
//			std::cout << ", " << *nPtr;
			TEST_ASSERT( neighborhoodAnswers[p++] == *nPtr );
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
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

	std::cout << std::endl;

      }

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
		std::cerr << "Unit test runtime ERROR: ut_QuickGridMPI_np2_3x3x2 only makes sense on 2 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 2 ./ut_QuickGridMPI_np2_3x3x2" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// This scale factor pushes the horizon just over the line so that 3 cells are included
	//  in the half neighborhood
	double SCALE=1.0;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	gridData =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);


	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
