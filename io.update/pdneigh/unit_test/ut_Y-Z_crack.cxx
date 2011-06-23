/*
 * ut_Y-Z_crack.cxx
 *
 *  Created on: Mar 1, 2011
 *      Author: jamitch
 *
 * REPLACES: utPd_Y-Z_crack.cxx
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "../BondFilter.h"
#include <iostream>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace PdBondFilter;
using namespace PDNEIGH;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

const size_t numProcs = 1;
const size_t myRank = 0;

const size_t nx = 4;
const size_t ny = 4;
const size_t nz = 3;
const size_t numCells = nx*ny*nz;
const double xStart = -2.5;
const double xLength = 5.0;
const double yStart = -2.5;
const double yLength = 5.0;
const double zStart = -0.5;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const double dx = xSpec.getCellSize();
const double dy = ySpec.getCellSize();
const double dz = zSpec.getCellSize();
const double _cellVolume = dx*dy*dz;
/*
 * function prototype
 */
FinitePlane getYZ_CrackPlane();

/*
 * Horizon
 */
const double horizon=1.1*sqrt( (2.0*dx)*(2.0*dx) +(2.0*dy)*(2.0*dy) +(2.0*dz)*(2.0*dz) );
/*
 * This demonstrates how the first and second coordinate
 * along an axis are computed
 */
//const double x0 = xStart+xSpec.getCellSize()/2.0;
//const double x1 = x0 + xSpec.getCellSize();
//const double y0 = yStart+ySpec.getCellSize()/2.0;
//const double y1 = y0 + ySpec.getCellSize();
//const double z0 = zStart+zSpec.getCellSize()/2.0;
//const double z1 = z0 + zSpec.getCellSize();


QUICKGRID::QuickGridData getGrid() {

	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData gridData =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	gridData=getLoadBalancedDiscretization(gridData);
	return gridData;
}


FinitePlane getYZ_CrackPlane() {

	/*
	 * Crack runs along y-axis
	 * Crack length along y-axis is 1.5 * dy
	 * Crack runs from bottom to top of plate in z-dir (length = 2 * dz)
	 */

	/*
	 * Lower left hand corner of crack plane when viewing down
	 * normal in the +dir
	 */
	const double x0 = xStart+xSpec.getCellSize()/2.0+1.5*dx;
	const double y0 = yStart+ySpec.getCellSize()/2.0;
	const double z0 = zStart+zSpec.getCellSize()/2.0;


	double n[3]; n[0]=-1.0;n[1]=0.0;n[2]=0.0;
	double r0[3]; r0[0]=x0; r0[1]=y0; r0[2]=z0;
	double ub[3]; ub[0]=0; ub[1]=1.0;ub[2]=0.0;
	double b=1.5*dy, a=2.0*dz;
	return FinitePlane(n,r0,ub,b,a);
}

void printNeighborhood(int numNeigh, int* neigh){
	for(int i=0;i<numNeigh;i++,neigh++){
		cout << ", " << *neigh;
	}
	cout << endl;
}

void assertNeighborhood(){
	QUICKGRID::QuickGridData gridData = getGrid();
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	int *neigh = list.get_neighborhood().get();

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 */
	BOOST_CHECK(numCells == gridData.numPoints);
	BOOST_CHECK(numCells == list.get_num_owned_points());

	/*
	 * Expected neighborhood data
	 *
	 * TO DO
	 * FINISH remaining points -- TEDIOUS
	 * The first 7 points are done below.  Need to finish all 48 points.
	 */
	int N[] = {21,21,21,21,29,26,26};
	int n0[] = { 1, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21, 24, 25, 28, 32, 33, 36, 37, 40, 41, 44 };
	int n1[] = { 0, 4, 5, 8, 9, 12, 13, 16, 17, 20, 21, 24, 25, 29, 32, 33, 36, 37, 40, 41, 45 };
	int n2[] = { 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 30, 34, 35, 38, 39, 42, 43, 46 };
	int n3[] = { 2, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 31, 34, 35, 38, 39, 42, 43, 47 };
	int n4[] = { 0, 1, 5, 8, 9, 10, 12, 13, 14, 16, 17, 20, 21, 24, 25, 26, 28, 29, 30, 32, 33, 36, 37, 40, 41, 42, 44, 45, 46 };
	int n5[] = { 0, 1, 4, 8, 9, 12, 13, 14, 16, 17, 20, 21, 24, 25, 28, 29, 30, 32, 33, 36, 37, 40, 41, 44, 45, 46 };
	int n6[] = { 2, 3, 7, 10, 11, 13, 14, 15, 18, 19, 22, 23, 26, 27, 29, 30, 31, 34, 35, 38, 39, 42, 43, 45, 46, 47 };
	int* NN[] = {n0,n1,n2,n3,n4,n5,n6};
	int NUMPOINTS = 7;


	for(int j=0;j<NUMPOINTS;j++){

		int gid=j;
		int numNeighAnswer = N[gid];
		int numNeigh = *neigh; neigh++;
//		cout << "gid, numNeigh = " << gid << ", " << numNeigh << endl;
//		printNeighborhood(numNeigh,neigh);
		BOOST_CHECK(numNeighAnswer==numNeigh);
		int *neighAnswer = NN[gid];
		for(int n=0;n<numNeigh;n++){
			BOOST_CHECK(*(neighAnswer+n)==*(neigh+n));
		}

		/*
		 * Move to next point
		 */
		neigh = neigh+numNeigh;
	}

//{
//		int gid=7;
//		int numNeigh = *neigh; neigh++;
//		cout << "gid, numNeigh = " << gid << ", " << numNeigh << endl;
//		printNeighborhood(numNeigh,neigh);
//
//	}


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "ut_Y-Z_crack" );
	proc->add(BOOST_TEST_CASE( &assertNeighborhood ));
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
#ifdef HAVE_MPI
		MPI_Init(&argc, &argv);
#endif
	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
#ifdef HAVE_MPI
	MPI_Finalize();
#endif

}
