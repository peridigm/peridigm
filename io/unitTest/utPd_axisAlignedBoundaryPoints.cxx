/*! \file utPd_axisAlignedBoundaryPoints.cxx */

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
#include "PdNeighborhood.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdNeighborhood.h"
#include "Field.h"
#include "PdVTK.h"
#include <set>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <utility>
#include <map>
#include <vector>
#include <cstdlib>
#include <set>
#include <time.h>
#include <tr1/memory>


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using namespace Field_NS;
using std::tr1::shared_ptr;
using namespace boost::unit_test;


const int numProcs=1;
const int myRank=0;

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
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=1.01*sqrt(pow(lX/nx,2)+pow(lY/ny,2)+pow(lZ/nz,2));


void axisAlignedMinimum() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::Z;
	std::tr1::shared_ptr<double> xPtr = decomp.myX;
	int numPoints = decomp.numPoints;
	BOOST_CHECK(numCells==numPoints);

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
	Pd_shared_ptr_Array<int> bcIds = PdNeighborhood::getPointsAxisAlignedMinimum(axis,xPtr,numPoints,horizon);
	BOOST_CHECK(32==bcIds.getSize());
	for(int *ids = bcIds.get();ids!=bcIds.end();ids++)
		BOOST_CHECK(setEnd != answerIds.find(*ids));

}

void axisAlignedMaximum() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::Z;
	std::tr1::shared_ptr<double> xPtr = decomp.myX;
	int numPoints = decomp.numPoints;
	BOOST_CHECK(numCells==numPoints);

	/*
	 * points at z-maximum end
	 */
	std::set<int> answerIds;
	for(int i=608;i<numPoints;i++){
		answerIds.insert(i);
	}
	std::set<int>::iterator setEnd = answerIds.end();
	/*
	 * This finds 2 planes of points (x-y plane) at the maximum value of z-end of the bar
	 */
	Pd_shared_ptr_Array<int> bcIds = PdNeighborhood::getPointsAxisAlignedMaximum(axis,xPtr,numPoints,horizon);
	BOOST_CHECK(32==bcIds.getSize());
	for(int *ids = bcIds.get();ids!=bcIds.end();ids++){
		BOOST_CHECK(setEnd != answerIds.find(*ids));
	}

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPd_axisAlignedBoundaryPoints" );
	proc->add(BOOST_TEST_CASE( &axisAlignedMinimum ));
	proc->add(BOOST_TEST_CASE( &axisAlignedMaximum ));
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
