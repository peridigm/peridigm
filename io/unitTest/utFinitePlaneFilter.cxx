/*
 * utFinitePlaneFilter.cxx
 *
 *  Created on: Dec 23, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdNeighborhood.h"
#include "PdBondFilter.h"
#include "PdVTK.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include <Teuchos_RCP.hpp>
#include <iostream>

using Teuchos::RCP;
using namespace PdQuickGrid;
using namespace PdBondFilter;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;

static int myRank = 0;
static int numProcs = 1;
const int nx = 2;
const int ny = 2;
const int nz = 2;
/*
 * NOTE THAT THIS MAKES edges of cube of length 1.0;
 */
const double cube_edge_length=2.0;
const double xStart = -cube_edge_length/nx/2.0;
const double xLength = cube_edge_length;
const double yStart = -cube_edge_length/ny/2.0;
const double yLength = cube_edge_length;
const double zStart = -cube_edge_length/nz/2.0;
const double zLength = cube_edge_length;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double SCALE=1.1*sqrt(3);
const double horizon = SCALE*xSpec.getCellSize();

PdGridData getGrid() {

	/*
	 * 2x2x2 Grid of points/cells
	 */
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	// This load-balances
	decomp = getLoadBalancedDiscretization(decomp);
	return decomp;
}

FinitePlane getCase_1a(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=-1.0/sqrt2;n[1]=1.0/sqrt2;n[2]=0.0;
	double r0[3]; r0[0]=0.0; r0[1]=0.0; r0[2]=0.0;
	double ua[3]; ua[0]=1.0/sqrt2; ua[1]=1.0/sqrt2;ua[2]=0.0;
	double a=1.0, b=1.0;
	return FinitePlane(n,r0,ua,a,b);
}


void case_1a() {

	PdGridData decomp = getGrid();
	FinitePlane plane = getCase_1a();
	RCP<BondFilter> filterPtr=RCP<BondFilter>(new FinitePlaneFilter(plane));

	/*
	 * Create KdTree; Since this is serial xOwned = xOverlap and numOwned = numOverlap
	 */
	std::tr1::shared_ptr<double> xOwnedPtr = decomp.myX;
	std::tr1::shared_ptr<double> xOverlapPtr = decomp.myX;
	size_t numOverlapPoints = decomp.numPoints;
	vtkSmartPointer<vtkUnstructuredGrid> overlapGrid = PdVTK::getGrid(xOverlapPtr,numOverlapPoints);
	vtkKdTreePointLocator* kdTree = vtkKdTreePointLocator::New();
	kdTree->SetDataSet(overlapGrid);

	{
		/*
		 * SANITY CHECK on Expected coordinates and IDs
		 */
		int ids[] = {0,1,2,3,4,5,6,7};
		double x[]={
				0.0,0.0,0.0,
				1.0,0.0,0.0,
				0.0,1.0,0.0,
				1.0,1.0,0.0,
				0.0,0.0,1.0,
				1.0,0.0,1.0,
				0.0,1.0,1.0,
				1.0,1.0,1.0
		};
		BOOST_CHECK(8==numOverlapPoints);
		for(int j=0;j<numOverlapPoints;j++){
			BOOST_CHECK(decomp.myGlobalIDs.get()[j]==ids[j]);
//			cout << "decomp.myX.get()[j*3+0], decomp.myX.get()[j*3+1], decomp.myX.get()[j*3+2] = "
//					<< decomp.myX.get()[j*3+0] << ", "
//					<< decomp.myX.get()[j*3+1] << ", "
//					<< decomp.myX.get()[j*3+2] << std::endl;
			BOOST_CHECK(decomp.myX.get()[j*3+0]==x[j*3+0]);
			BOOST_CHECK(decomp.myX.get()[j*3+1]==x[j*3+1]);
			BOOST_CHECK(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}

	}

	{
		/*
		 * look at neighborhood of id = 0
		 */
		int id=0;
		vtkIdList* kdTreeList = vtkIdList::New();
		/*
		 * Note that list returned includes this point *
		 */
		double *x = decomp.myX.get()+3*id;
		kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
		size_t listSize = filterPtr->filterListSize(kdTreeList,x,id,decomp.myX.get());
//		cout << "id, listSize = " << id << ", " << listSize << std::endl;
		kdTreeList->Delete();
	}

	{
		/*
		 * look at neighborhood of id = 0
		 */
		int id=1;
		vtkIdList* kdTreeList = vtkIdList::New();
		/*
		 * Note that list returned includes this point *
		 */
		double *x = decomp.myX.get()+3*id;
		kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
		size_t listSize = filterPtr->filterListSize(kdTreeList,x,id,decomp.myX.get());
		cout << "id, listSize = " << id << ", " << listSize << std::endl;
		kdTreeList->Delete();
	}

//	decomp = createAndAddNeighborhood(decomp,horizon,filterPtr);


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "utFinitePlaneFilter" );
	proc->add(BOOST_TEST_CASE( &case_1a ));
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
