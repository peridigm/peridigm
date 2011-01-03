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

	{
		/*
		 * SANITY CHECK on Expected coordinates and IDs
		 */
		size_t numOverlapPoints = decomp.numPoints;
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
			BOOST_CHECK(decomp.myX.get()[j*3+0]==x[j*3+0]);
			BOOST_CHECK(decomp.myX.get()[j*3+1]==x[j*3+1]);
			BOOST_CHECK(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}

	}

	return decomp;
}

FinitePlane getCase_1a(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=-1.0/sqrt2;n[1]=1.0/sqrt2;n[2]=0.0;
	double r0[3]; r0[0]=0.0; r0[1]=0.0; r0[2]=0.0;
	double ua[3]; ua[0]=1.0/sqrt2; ua[1]=1.0/sqrt2;ua[2]=0.0;
	double a=sqrt(2), b=1.0;
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

	/*
	 * ANSWERS for each ID
	 * list size for each point
	 */
	// known local ids
	size_t ids[] = {0,1,2,3,4,5,6,7};
	bool markForExclusion[8];
	// Expected: filter should evaluate this list size for each id
	size_t size[] = {4,2,2,4,4,2,2,4};
	// Expected: filter should return these flags for each local id
	bool n0[]={1,1,1,0,0,1,1,0};
	bool n1[]={1,1,1,1,1,0,1,1};
	bool n2[]={1,1,1,1,1,1,0,1};
	bool n3[]={0,1,1,1,0,1,1,0};
	bool n4[]={0,1,1,0,1,1,1,0};
	bool n5[]={1,0,1,1,1,1,1,1};
	bool n6[]={1,1,0,1,1,1,1,1};
	bool n7[]={0,1,1,0,0,1,1,1};
	bool * expectedFlags[] = {n0,n1,n2,n3,n4,n5,n6,n7};
	{
		for(size_t i=0;i<8;i++){
			/*
			 * look at neighborhood of id = 0
			 */
			size_t id=ids[i];
			vtkIdList* kdTreeList = vtkIdList::New();
			/*
			 * Note that list returned includes this point *
			 */
			double *x = decomp.myX.get()+3*id;
			kdTree->FindPointsWithinRadius(horizon, x, kdTreeList);
			/*
			 * Use filter to decide length of neighborhood list (will include extra entry for 'number of neighbors: numNeigh')
			 */
			size_t listSize = filterPtr->filterListSize(kdTreeList,x,id,decomp.myX.get());
			BOOST_CHECK(listSize==size[id]);

			/*
			 * Now determine which points are included
			 */
			filterPtr->filterBonds(kdTreeList,x,id,decomp.myX.get(),markForExclusion);
			bool *flags = expectedFlags[i];
			/*
			 * Assert flags
			 */
			for(int j=0;j<8;j++){
//				cout << "filter flag, expected flag = " << *(markForExclusion+j) << ", " << *(flags+j) << endl;
				BOOST_CHECK(*(flags+j)==*(markForExclusion+j));
			}

			// delete tree list for this point
			kdTreeList->Delete();
		}
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
