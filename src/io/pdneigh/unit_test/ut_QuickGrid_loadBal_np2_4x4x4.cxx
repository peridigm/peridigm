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


#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "PdutMpiFixture.h"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include <iostream>

using namespace Pdut;
using std::tr1::shared_ptr;
using std::cout;

static size_t myRank;
static size_t numProcs;
const size_t nx = 4;
const size_t ny = 4;
const size_t nz = 4;
const double xStart = 0.0;
const double xLength = 100.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = 0.0;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;

QUICKGRID::QuickGridData getGrid() {
	// This scale factor pushes the horizon just over the line so that 2 cells are included
	//  in the half neighborhood
	double SCALE=1.4;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	// This rload-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}

TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np2_4x4x4, p0) {


        if(myRank != 0){
           return;
        }
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
	TEST_ASSERT(32 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=0 end should be on this processor
	 */
	int ids[] = {0,1,4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33,36,37,40,41,44,45,48,49,52,53,56,57,60,61};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P0 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
//	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}


}

TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np2_4x4x4, p1) {

	

        if(myRank != 1){
           return;
        }


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
	TEST_ASSERT(32 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {34,35,38,39,42,43,46,47,50,51,54,55,58,59,62,63,2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P1 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
//	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
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
		std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np2_4x4x4 only makes sense on 2 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 2 ./ut_QuickGrid_loadBal_np2_4x4x4" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF

        return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}


