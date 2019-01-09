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
const size_t nx = 3;
const size_t ny = 3;
const size_t nz = 3;
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

QUICKGRID::QuickGridData getGrid() {

	// This scale factor pushes the horizon just over the line so that 2 cells are included
	//  in the half neighborhood
	double SCALE=1.5;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}


TEUCHOS_UNIT_TEST( QuickGridMPI_np3, p0Test) {

	if (myRank == 0){

        QUICKGRID::QuickGridData decomp = getGrid();

        
	TEST_ASSERT(0 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(9 == myNumPoints);

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 0;
	std::cout << "proc = 0 global ids = ";
	for(size_t id=start;id<decomp.numPoints+start;id++,gIdsPtr++){
//		TEST_ASSERT( *gIdsPtr == id );
		std::cout << *gIdsPtr << ", ";
	}
	std::cout << std::endl;

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
	TEST_ASSERT( sizeNeighborList == decomp.sizeNeighborhoodList );

	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = decomp.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	// each cell in neighborhood has the entire list of cells in mesh except itself
	gIdsPtr = gIds.get();
	for(int id=0;id<myNumPoints;id++, gIdsPtr++){
		int numNeigh = *nPtr; nPtr++;
		TEST_ASSERT(26 == numNeigh);

		for(int i=0;i<27;i++){
			if(*gIdsPtr != i) {
				TEST_ASSERT(i == *nPtr);
				nPtr++;
			}
		}

	}

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

    }


}

TEUCHOS_UNIT_TEST( QuickGridMPI_np3, p1Test) 
{
	 if (myRank == 1){

        QUICKGRID::QuickGridData decomp = getGrid();

       
	TEST_ASSERT(1 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
        TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(9 == myNumPoints);

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 9;
	std::cout << "proc = 1 global ids = ";
	for(size_t id=start;id<decomp.numPoints+start;id++,gIdsPtr++){
//		TEST_ASSERT( *gIdsPtr == id );
		std::cout << *gIdsPtr << ", ";
	}
	std::cout << std::endl;

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
	TEST_ASSERT( sizeNeighborList == decomp.sizeNeighborhoodList );

	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = decomp.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	// each cell in neighborhood has the entire list of cells in mesh except itself
	gIdsPtr = gIds.get();
	for(int id=0;id<myNumPoints;id++, gIdsPtr++){
		int numNeigh = *nPtr; nPtr++;
		TEST_ASSERT(26 == numNeigh);

		for(int i=0;i<27;i++){
			if(*gIdsPtr != i) {
				TEST_ASSERT(i == *nPtr);
				nPtr++;
			}
		}

	}

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

    }

}


TEUCHOS_UNIT_TEST( QuickGridMPI_np3, p2Test) 
{
	if (myRank == 2){

        QUICKGRID::QuickGridData decomp = getGrid();
	TEST_ASSERT(2 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(9 == myNumPoints);

	// assert length of neighborlist
	// sizeNeighborList = myNumCells + myNumCells*numNeighbors
	int sizeNeighborList = myNumPoints + myNumPoints*(27-1);
	TEST_ASSERT( sizeNeighborList == decomp.sizeNeighborhoodList );

	/*
	 * Assert global ids on this processor
	 */
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int start = 18;
	std::cout << "proc = 2 global ids = ";
	for(size_t id=start;id<decomp.numPoints+start;id++,gIdsPtr++){
//		TEST_ASSERT( *gIdsPtr == id );
		std::cout << *gIdsPtr << ", ";
	}
	std::cout << std::endl;

	// assert neighbor lists; in this case each point has all points
	shared_ptr<int> neighborList = decomp.neighborhood;
	int *nPtr = neighborList.get();

	// This iterates through all cell neighborhoods
	// each cell in neighborhood has the entire list of cells in mesh except itself
	gIdsPtr = gIds.get();
	for(int id=0;id<myNumPoints;id++, gIdsPtr++){
		int numNeigh = *nPtr; nPtr++;
		TEST_ASSERT(26 == numNeigh);

		for(int i=0;i<27;i++){
			if(*gIdsPtr != i) {
				TEST_ASSERT(i == *nPtr);
				nPtr++;
			}
		}

	}

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

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
	 * This test only make sense for numProcs == 3
	 */
	if(3 != numProcs){
		std::cerr << "Unit test runtime ERROR: ut_QuickGridMPI_np3 only makes sense on 3 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 3 ./ut_QuickGridMPI_np3" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}


	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}


