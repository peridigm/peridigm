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

#include <iostream>
#include <cmath>
#include <map>
#include <set>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Array.h"
#include "PdZoltan.h"
#include "BondFilter.h"
#include "NeighborhoodList.h"
#include "QuickGrid.h"

#include "Field.h"
#include "PdutMpiFixture.h"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

using std::shared_ptr;

using std::cout;
using std::set;
using std::map;

int dimension_answer=3;
size_t globalNumPoints_answer=2;
size_t numPoints_answer=2;
int sizeNeighborhoodList_answer=4;
int gids_answer[] = {0,1};
int neighborhood_answer[] = {1,1,1,0};
int neighborhoodPtr_answer[] = {0,2};

using QUICKGRID::QuickGridMeshGenerationIterator;

QUICKGRID::QuickGridData getGrid(const string json_filename, int numProcs, int myRank) {
	shared_ptr<QuickGridMeshGenerationIterator> g = QUICKGRID::getMeshGenerator(numProcs,json_filename);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, *g);

	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}


TEUCHOS_UNIT_TEST(ReLoadBalance, TwoPointReloadBalanceTest) {

        Teuchos::RCP<Epetra_Comm> comm;

         #ifdef HAVE_MPI
              
              comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
         #else
              comm = rcp(new Epetra_SerialComm);
         #endif

         int numProcs = comm->NumProc();
         int myRank   = comm->MyPID();



        const int* gids;
        const int* neighborhood;
        const int* neighborhoodPtr;
        
	const string json_file="input_files/ut_twoPointReLoadBalance.json";
	QUICKGRID::QuickGridData decomp_1 = getGrid(json_file, numProcs, myRank);

        std::cout << "\ttwoPointReloadBalance: UNBALANCED MESH assert_grid()\n" << "\n";

        TEST_ASSERT(dimension_answer==decomp_1.dimension);
	TEST_ASSERT(globalNumPoints_answer==decomp_1.globalNumPoints);
	TEST_ASSERT(numPoints_answer==decomp_1.numPoints);
	TEST_ASSERT(sizeNeighborhoodList_answer==decomp_1.sizeNeighborhoodList);

        gids=decomp_1.myGlobalIDs.get();
	for(size_t n=0;n<numPoints_answer;n++,gids++)
            TEST_ASSERT(gids_answer[n]==*gids);

        neighborhood = decomp_1.neighborhood.get();
	for(int n=0;n<sizeNeighborhoodList_answer;n++,neighborhood++)
	    TEST_ASSERT(neighborhood_answer[n]==*neighborhood);

        neighborhoodPtr = decomp_1.neighborhoodPtr.get();
	for(size_t n=0;n<numPoints_answer;n++,neighborhoodPtr++)
	    TEST_ASSERT(neighborhoodPtr_answer[n]==*neighborhoodPtr);
        
        
	

	/*
	 * Reload balance
	 */
	QUICKGRID::QuickGridData decomp_2 = PDNEIGH::getLoadBalancedDiscretization(decomp_1);

        std::cout << "\ttwoPointReloadBalance: LOADBALANCED MESH assert_grid()\n" << "\n";

        TEST_ASSERT(dimension_answer==decomp_1.dimension);
	TEST_ASSERT(globalNumPoints_answer==decomp_1.globalNumPoints);
	TEST_ASSERT(numPoints_answer==decomp_1.numPoints);
	TEST_ASSERT(sizeNeighborhoodList_answer==decomp_1.sizeNeighborhoodList);


        gids=decomp_1.myGlobalIDs.get();
	for(size_t n=0;n<numPoints_answer;n++,gids++)
            TEST_ASSERT(gids_answer[n]==*gids);

        neighborhood = decomp_1.neighborhood.get();
	for(int n=0;n<sizeNeighborhoodList_answer;n++,neighborhood++)
	    TEST_ASSERT(neighborhood_answer[n]==*neighborhood);

        neighborhoodPtr = decomp_1.neighborhoodPtr.get();
	for(size_t n=0;n<numPoints_answer;n++,neighborhoodPtr++)
	    TEST_ASSERT(neighborhoodPtr_answer[n]==*neighborhoodPtr);

}


int main
(
		int argc,
		char* argv[]
)
{

	Teuchos::GlobalMPISession mpiSession(&argc, &argv);
	// Initialize UTF

	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}

