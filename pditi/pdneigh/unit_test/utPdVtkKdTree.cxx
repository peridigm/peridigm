/*
 * utPdVtkKdTree.cxx
 *
 *  Created on: Feb 17, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include <iostream>
#include <set>

#include "vtk/PdVTK.h"
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "../BondFilter.h"

using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;


const double xStart = 0.0;
const double xLength = 3.0;
const double yStart = 0.0;
const double yLength = 3.0;
const double zStart = 0.0;
const double zLength = 3.0;


void knownNeighborhood()
{
	const size_t nx = 3;
	const size_t ny = 3;
	const size_t nz = 1;
	const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
	const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
	const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
	const size_t numCells = nx*ny*nz;
	const double SCALE=1.05*sqrt(2);
	const double horizon = SCALE*xSpec.getCellSize();


	size_t myRank = 0;
	size_t numProcs = 1;

	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

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
	const size_t nx = 3;
	const size_t ny = 3;
	const size_t nz = 1;
	const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
	const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
	const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
	const size_t numCells = nx*ny*nz;

	const double SCALE=1.05*sqrt(2);
	const double horizon = SCALE*xSpec.getCellSize();


	size_t myRank = 0;
	size_t numProcs = 1;

	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec, QUICKGRID::SphericalNorm);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

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
	const size_t nx = 10;
	const size_t ny = 10;
	const size_t nz = 10;
	const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
	const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
	const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
	const size_t numCells = nx*ny*nz;
	const double SCALE=1.05*sqrt(3);
	const double horizon = SCALE*xSpec.getCellSize();


	size_t myRank = 0;
	size_t numProcs = 1;

	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec, QUICKGRID::SphericalNormFunction);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

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
	test_suite* proc = BOOST_TEST_SUITE( "utPdVtkKdTree" );
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
