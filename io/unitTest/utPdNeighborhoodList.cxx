/*! \file utPdNeighborhoodList.cxx */

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
#include "PdBondFilter.h"
#include <Teuchos_RCP.hpp>
#include <set>

using namespace PdQuickGrid;
using namespace PdNeighborhood;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using Teuchos::RCP;
using PdBondFilter::BondFilter;

const int nx = 1;
const int ny = 1;
const int nz = 11;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = 0.0;
const double zLength = 1.0;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;



void axialBarLinearSpacing() {

	int numProcs=1;
	int myRank=0;
	double SCALE=3.1;
	double horizon = SCALE*zSpec.getCellSize();
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	int numPoints = decomp.numPoints;
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);

	int n0[]  = {0,1,2,3};
	int n1[]  = {0,1,2,3,4};
	int n2[]  = {0,1,2,3,4,5};
	int n3[]  = {0,1,2,3,4,5,6};
	int n4[]  = {1,2,3,4,5,6,7};
	int n5[]  = {2,3,4,5,6,7,8};
	int n6[]  = {3,4,5,6,7,8,9};
	int n7[]  = {4,5,6,7,8,9,10};
	int n8[]  = {5,6,7,8,9,10};
	int n9[]  = {6,7,8,9,10};
	int n10[] = {7,8,9,10};
	int* nPtr[]  = {n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10};
	int s[] = {4,5,6,7,7,7,7,7,6,5,4};

	/*
	 * Construct neighborhood
	 * NOTE THAT: Since this is a serial test, numPoints = numOverlapPoints; and x = xOverlap
	 */
	RCP<BondFilter> bondFilterPtr(new PdBondFilter::BondFilterDefault(true));
	NeighborhoodList list = PdNeighborhood::getNeighborhoodList(horizon,decomp.numPoints,decomp.numPoints,decomp.myX,decomp.myX,*bondFilterPtr);
	BOOST_CHECK(list.getNumPoints() == numPoints);
	int size = 0;
	for(int n=0;n<numPoints;n++)
		size+=list.getNumNeigh(n);
	BOOST_CHECK(65 == size);

	for(int n=0;n<numPoints;n++){
		size=list.getNumNeigh(n);
		BOOST_CHECK(s[n] == size);
		// skip first entry -- thats the size (numNeigh)
		const int *neigh = list.getNeighborhood(n)+1;
		std::set<int> found(neigh,neigh+size);
		int* N = nPtr[n];
		for(int j=0;j<size;j++){
			BOOST_CHECK(found.end()!=found.find(N[j]));
		}
	}

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPdNeighborhoodList" );
	proc->add(BOOST_TEST_CASE( &axialBarLinearSpacing ));
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
