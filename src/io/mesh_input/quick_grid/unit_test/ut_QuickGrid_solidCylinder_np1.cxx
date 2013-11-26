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
#include "Field.h"
#include "../QuickGrid.h"
#include <valarray>
#include <iostream>
#include <cmath>



using std::tr1::shared_ptr;
using std::cout;

static size_t myRank;
static size_t numProcs;
const double cylinderRadius = 1.0;
const size_t numRings = 2;
const double zStart = 0.0;
const double cylinderLength = 1.0;


QUICKGRID::TensorProductSolidCylinder getMeshGenerator(){

	QUICKGRID::TensorProductSolidCylinder meshGen(numProcs,cylinderRadius,numRings,zStart,cylinderLength);
	return meshGen;
}


TEUCHOS_UNIT_TEST( QuickGrid_solidCylinder_np1, RunTest) {

	QUICKGRID::TensorProductSolidCylinder meshGen = getMeshGenerator();
	TEST_ASSERT(57==meshGen.getNumGlobalCells());
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, meshGen);
	int numPoints = decomp.numPoints;
	TEST_ASSERT(57==numPoints);

	/*
	 * Assert that there are 3 core points
	 */
	double *x = decomp.myX.get();
	double dz = cylinderLength/3;
	double z0 = zStart + dz/2;
	const double tolerance = 1.0e-15;
	TEST_FLOATING_EQUALITY(*(x+0),0.0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+1),0.0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+2),z0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+19*3+0),0.0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+19*3+1),0.0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+19*3+2),z0+dz,tolerance);
	TEST_FLOATING_EQUALITY(*(x+2*19*3+0),0.0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+2*19*3+1),0.0,tolerance);
	TEST_FLOATING_EQUALITY(*(x+2*19*3+2),z0+2.0*dz,tolerance);

	/*
	 * Check volume of core of 3 core points
	 */
	double rI = cylinderRadius/numRings/sqrt(M_PI);
	double vol = M_PI*rI*rI*dz;
	double *v = decomp.cellVolume.get();
	TEST_FLOATING_EQUALITY(*(v+0),vol,tolerance);
	TEST_FLOATING_EQUALITY(*(v+19),vol,tolerance);
	TEST_FLOATING_EQUALITY(*(v+2*19),vol,tolerance);

	/*
	 * Check global assigned "global ids"
	 */
	int *ids = decomp.myGlobalIDs.get();
	for(int i=0;ids!=decomp.myGlobalIDs.get()+numPoints;ids++,i++){
		TEST_ASSERT(*ids==i);
	}
}



int main
(
		int argc,
		char* argv[]
)
{

	// These are static (file scope) variables
	myRank = 0;
	numProcs = 1;

	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
