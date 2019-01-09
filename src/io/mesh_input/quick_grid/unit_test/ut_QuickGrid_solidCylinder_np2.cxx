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
#include "../QuickGridData.h"
#include "Field.h"
#include "PdZoltan.h"
#include "PdutMpiFixture.h"
#include <valarray>
#include <iostream>
#include <cmath>


using std::shared_ptr;
using std::shared_ptr;

using namespace Pdut;

using std::cout;


static size_t myRank;
static size_t numProcs;
const double cylinderRadius = 1.0;
const size_t numRings = 2;
const double zStart = 0.0;
const double cylinderLength = 1.0;
const double horizon = .6;


QUICKGRID::TensorProductSolidCylinder getMeshGenerator(){

	QUICKGRID::TensorProductSolidCylinder meshGen(numProcs,cylinderRadius,numRings,zStart,cylinderLength);
	return meshGen;
}


TEUCHOS_UNIT_TEST( QuickGrid_solidCylinder_np2, RunTest) {

	QUICKGRID::TensorProductSolidCylinder meshGen = getMeshGenerator();
	TEST_ASSERT(57==meshGen.getNumGlobalCells());
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, meshGen);
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	TEST_ASSERT(57==decomp.globalNumPoints);
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
		std::cerr << "Unit test runtime ERROR: ut_QuickGrid_solidCylinder_np2 only makes sense on 2 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 2 ./ut_QuickGrid_solidCylinder_np2" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
