/*! \file utPdVtkKdTree.cxx */

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
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdutMpiFixture.h"
#include "PdVTK.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include <iostream>
#include <functional>
#include <valarray>
#include <set>

using namespace PdQuickGrid;
using namespace Pdut;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;

using std::valarray;
using std::binary_function;


struct Distance : public binary_function< valarray<double>, valarray<double>, double > {
	double operator()(const valarray<double>& u, const valarray<double>& v){
		double dx = v[0]-u[0];
		double dy = v[1]-u[1];
		double dz = v[2]-u[2];
		return sqrt(dx*dx+dy*dy+dz*dz);
	}
};

struct Dot : public binary_function< valarray<double>, valarray<double>, double > {
	double operator()(const valarray<double>& u, const valarray<double>& v){
		return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	}
};

struct Cross : public binary_function< valarray<double>, valarray<double>, valarray<double> > {
	valarray<double> operator()(const valarray<double>& u, const valarray<double>& v){
		double r0=-u[2]* v[1] + u[1] * v[2];
		double r1= u[2]* v[0] - u[0] * v[2];
		double r2=-u[1]* v[0] + u[0] * v[1];
		valarray<double> r(3); r[0]=r0;r[1]=r1;r[2]=r2;
		return r;
	}
};

const double xStart = 0.0;
const double xLength = 3.0;
const double yStart = 0.0;
const double yLength = 3.0;
const double zStart = 0.0;
const double zLength = 3.0;


void knownNeighborhood()
{
	const int nx = 3;
	const int ny = 3;
	const int nz = 1;
	const PdQPointSet1d xSpec(nx,xStart,xLength);
	const PdQPointSet1d ySpec(ny,yStart,yLength);
	const PdQPointSet1d zSpec(nz,zStart,zLength);
	const int numCells = nx*ny*nz;
	const double SCALE=1.05*sqrt(2);
	const double horizon = SCALE*xSpec.getCellSize();


	int myRank = 0;
	int numProcs = 1;

	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	/*
	 * problem dimension is 3
	 */
	BOOST_CHECK(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	BOOST_CHECK(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	BOOST_CHECK(numCells == myNumPoints);

	std::tr1::shared_ptr<double> yPtr = decomp.myX;
	vtkSmartPointer<vtkUnstructuredGrid> gBalanced = PdVTK::getGrid(yPtr,myNumPoints);

	vtkKdTreePointLocator* kdTree = vtkKdTreePointLocator::New();
	kdTree->SetDataSet(gBalanced);


	/*
	 * Neighborhood answers
	 */
	int numNeigh[] = {3,5,3,5,8,5,3,5,3};
	int n0[] = {1,4,3};
	int n1[] = {0,2,5,4,3};
	int n2[] = {1,4,5};
	int n3[] = {0,1,4,7,6};
	int n4[] = {0,1,2,3,5,6,7,8};
	int n5[] = {2,1,4,7,8};
	int n6[] = {3,4,7};
	int n7[] = {6,3,4,5,8};
	int n8[] = {7,4,5};
	int* list[] = {n0,n1,n2,n3,n4,n5,n6,n7,n8};

	double *x = decomp.myX.get();
	for(int c=0;c<myNumPoints;c++,x+=3){

		vtkIdList* kdTreeList = vtkIdList::New();
		kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
		BOOST_CHECK(numNeigh[c] == kdTreeList->GetNumberOfIds()-1);
		std::set<int> answers(list[c],list[c]+numNeigh[c]);
		std::set<int>::iterator end = answers.end();
		for(int i=0;i<kdTreeList->GetNumberOfIds();i++){
			int uid = kdTreeList->GetId(i);
			if(uid == c) continue;
			BOOST_CHECK(end != answers.find(uid));
		}
		kdTreeList->Delete();
	}

	kdTree->Delete();
}

void sphericalNeighborhood()
{
	const int nx = 3;
	const int ny = 3;
	const int nz = 1;
	const PdQPointSet1d xSpec(nx,xStart,xLength);
	const PdQPointSet1d ySpec(ny,yStart,yLength);
	const PdQPointSet1d zSpec(nz,zStart,zLength);
	const int numCells = nx*ny*nz;

	const double SCALE=1.05*sqrt(2);
	const double horizon = SCALE*xSpec.getCellSize();


	int myRank = 0;
	int numProcs = 1;

	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec, PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	/*
	 * problem dimension is 3
	 */
	BOOST_CHECK(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	BOOST_CHECK(numCells == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	BOOST_CHECK(numCells == myNumPoints);


	/*
	 * Neighborhood answers
	 */
	int numNeigh[] = {3,5,3,5,8,5,3,5,3};
	int n0[] = {1,4,3};
	int n1[] = {0,2,5,4,3};
	int n2[] = {1,4,5};
	int n3[] = {0,1,4,7,6};
	int n4[] = {0,1,2,3,5,6,7,8};
	int n5[] = {2,1,4,7,8};
	int n6[] = {3,4,7};
	int n7[] = {6,3,4,5,8};
	int n8[] = {7,4,5};
	int* list[] = {n0,n1,n2,n3,n4,n5,n6,n7,n8};

	int *neighborhoodList = decomp.neighborhood.get();
	for(int c=0;c<myNumPoints;c++){

		BOOST_CHECK(numNeigh[c] == (*neighborhoodList));
		neighborhoodList++;
		std::set<int> answers(list[c],list[c]+numNeigh[c]);
		std::set<int>::iterator end = answers.end();
		for(int i=0;i<numNeigh[c];i++){
			int uid = *neighborhoodList; neighborhoodList++;
			BOOST_CHECK(end != answers.find(uid));
		}
	}

}

void compareKdTreeWithSphericalNorm()
{
	const int nx = 10;
	const int ny = 10;
	const int nz = 10;
	const PdQPointSet1d xSpec(nx,xStart,xLength);
	const PdQPointSet1d ySpec(ny,yStart,yLength);
	const PdQPointSet1d zSpec(nz,zStart,zLength);
	const int numCells = nx*ny*nz;
	const double SCALE=1.05*sqrt(3);
	const double horizon = SCALE*xSpec.getCellSize();


	int myRank = 0;
	int numProcs = 1;

	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec, PdQuickGrid::SphericalNormFunction);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	/*
	 * problem dimension is 3
	 */
	BOOST_CHECK(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	BOOST_CHECK(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	BOOST_CHECK(numCells == myNumPoints);

	std::tr1::shared_ptr<double> yPtr = decomp.myX;
	vtkSmartPointer<vtkUnstructuredGrid> gBalanced = PdVTK::getGrid(yPtr,myNumPoints);

	vtkKdTreePointLocator* kdTree = vtkKdTreePointLocator::New();
	kdTree->SetDataSet(gBalanced);

	double *x = decomp.myX.get();
	int *neighborhoodList = decomp.neighborhood.get();
	for(int c=0;c<myNumPoints;c++,x+=3){
		int numNeigh = *neighborhoodList; neighborhoodList++;
		vtkIdList* kdTreeList = vtkIdList::New();
		kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
		BOOST_CHECK(numNeigh == kdTreeList->GetNumberOfIds()-1);
		std::set<int> answers(neighborhoodList,neighborhoodList+numNeigh);
		std::set<int>::iterator end = answers.end();
		for(int i=0;i<numNeigh;i++){
			int uid = kdTreeList->GetId(i);
			if(uid == c) continue;
			BOOST_CHECK(end != answers.find(uid));
		}
		neighborhoodList+=numNeigh;
		kdTreeList->Delete();
	}

	kdTree->Delete();
}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "utPdVtkTree" );
	proc->add(BOOST_TEST_CASE( &knownNeighborhood ));
	proc->add(BOOST_TEST_CASE( &sphericalNeighborhood ));
	proc->add(BOOST_TEST_CASE( &compareKdTreeWithSphericalNorm ));
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
