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
#include "quick_grid/QuickGrid.h"
#include "Sortable.h"
#include "../NeighborhoodList.h"
#include <set>
#include <utility>
#include <map>
#include <vector>
#include <cstdlib>
#include <set>
#include <time.h>


using UTILITIES::CartesianComponent;
using namespace PDNEIGH;
using std::tr1::shared_ptr;


const size_t numProcs=1;
const size_t myRank=0;

const int nx = 4;
const int ny = nx;
const double lX = 1.0;
const double lY = lX;
const double lZ = 10.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const int nz = (int)(lZ * nx / lX);
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;
const double horizon=1.01*sqrt(pow(lX/nx,2)+pow(lY/ny,2)+pow(lZ/nz,2));

TEUCHOS_UNIT_TEST(AxisAlignedBoundaryPoints, AxisAlignedMinimumTest) {

	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,QUICKGRID::SphericalNorm);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	CartesianComponent axis = UTILITIES::Z;
	std::tr1::shared_ptr<double> xPtr = decomp.myX;
	size_t numPoints = decomp.numPoints;
	TEST_ASSERT(numCells==numPoints);

	/*
	 * points at z minimum end
	 */
	std::set<int> answerIds;
	for(int i=0;i<32;i++){
		answerIds.insert(i);
	}
	std::set<int>::iterator setEnd = answerIds.end();

	/*
	 * This finds 2 planes of points (x-y plane) at the minimum value of z-end of the bar
	 */
	Array<int> bcIds = UTILITIES::getPointsAxisAlignedMinimum(axis,xPtr,numPoints,horizon);
	TEST_ASSERT(32==bcIds.get_size());
	for(int *ids = bcIds.get();ids!=bcIds.end();ids++)
		TEST_ASSERT(setEnd != answerIds.find(*ids));

}

TEUCHOS_UNIT_TEST(AxisAlignedBoundaryPoints, AxisAlignedMaximumTest) {

	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,QUICKGRID::SphericalNorm);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	CartesianComponent axis = UTILITIES::Z;
	std::tr1::shared_ptr<double> xPtr = decomp.myX;
	size_t numPoints = decomp.numPoints;
	TEST_ASSERT(numCells==numPoints);

	/*
	 * points at z-maximum end
	 */
	std::set<int> answerIds;
	for(size_t i=608;i<numPoints;i++){
		answerIds.insert(i);
	}
	std::set<int>::iterator setEnd = answerIds.end();
	/*
	 * This finds 2 planes of points (x-y plane) at the maximum value of z-end of the bar
	 */
	Array<int> bcIds = UTILITIES::getPointsAxisAlignedMaximum(axis,xPtr,numPoints,horizon);
	TEST_ASSERT(32==bcIds.get_size());
	for(int *ids = bcIds.get();ids!=bcIds.end();ids++){
		TEST_ASSERT(setEnd != answerIds.find(*ids));
	}

}

int main( int argc, char* argv[] ) {
  
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
