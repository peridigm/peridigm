/*
 * ut_Y-Z_crack_jacobian_np4.cxx
 *
 *  Created on: Mar 1, 2011
 *      Author: jamitch
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

#include "PdutMpiFixture.h"
#include "vtk/PdVTK.h"
#include <iostream>
#include <set>
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

using namespace Pdut;
using std::cout;
using std::set;

static int myRank;
static int numProcs;
const int nx = 4;
const int ny = 4;
const int nz = 3;
const int numCells = nx*ny*nz;
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
	/*
	 * Write file for debugging
	 */
	/*
	const FieldSpec myRankSpec(FieldSpec::DEFAULT_FIELDTYPE,FieldSpec::SCALAR,"MyRank");
	Field<double> X(COORD3D,gridData.myX,gridData.numPoints);
	Field<int> rankField(myRankSpec,gridData.numPoints);
	rankField.setValue(myRank);

	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(gridData.myX.get(), gridData.numPoints);
	PdVTK::writeField(grid,X);
	PdVTK::writeField(grid,rankField);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer= PdVTK::getWriter("ut_Y-Z_crack_jacobian_np4.pvtu", numProcs, myRank, PdVTK::vtkASCII);
	PdVTK::write(writer,grid);
	*/

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

void assertNeighborhood_p0(){
	QUICKGRID::QuickGridData gridData = getGrid();
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane,true));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	int *neigh = list.get_neighborhood().get();
	int *gids = gridData.myGlobalIDs.get();

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	BOOST_CHECK(12 == gridData.numPoints);
	BOOST_CHECK(12 == list.get_num_owned_points());
	/*
	 * GIDS ON THIS PROCESSOR
	 *
	 */
	int GIDS[] = {2,3,6,7,18,19,22,23,34,35,38,39};
	/*
	 * Expected neighborhood data
	 *
	 * TO DO
	 * FINISH remaining points -- TEDIOUS
	 * The first 3 points are done below.  Need to finish all 12 points on this processor
	 */
	int N[] = {22,22,27,30,22};
	int n2[]  = { 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 30, 34, 35, 38, 39, 42, 43, 46 };
	int n3[]  = { 2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 31, 34, 35, 38, 39, 42, 43, 47 };
	int n6[]  = { 2, 3, 6, 7, 10, 11, 13, 14, 15, 18, 19, 22, 23, 26, 27, 29, 30, 31, 34, 35, 38, 39, 42, 43, 45, 46, 47 };
	int n7[]  = { 2, 3, 6, 7, 18, 19, 22, 23, 34, 35, 38, 39, 25, 29,  9, 13, 41, 45, 42, 43, 46, 47, 10, 11, 14, 15, 26, 27, 30, 31 };
	int n18[] = { 2, 3, 6,  7, 18, 19, 22, 23, 34, 35, 38, 39, 42, 43, 46, 10, 11, 14, 26, 27, 30, 31 };
	int* NN[] = {n2,n3,n6,n7,n18};
	int NUMPOINTS = 5;


	for(int j=0;j<NUMPOINTS;j++,gids++){

		BOOST_CHECK(GIDS[j]==*gids);
		int numNeighAnswer = N[j];
		int numNeigh = *neigh; neigh++;
//		cout << "rank, gid, numNeigh = " << myRank << ", " << *gids << ", " << numNeigh << endl;
//		printNeighborhood(numNeigh,neigh);
		BOOST_CHECK(numNeighAnswer==numNeigh);
//		int *neighAnswer = NN[j];
		set<int> neighAnswer(NN[j],NN[j]+numNeigh);
		set<int>::iterator end = neighAnswer.end();
		for(int n=0;n<numNeigh;n++){
			BOOST_CHECK(end != neighAnswer.find(*(neigh+n)));
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

void assertNeighborhood_p1(){
	QUICKGRID::QuickGridData gridData = getGrid();
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	BOOST_CHECK(12 == gridData.numPoints);
	BOOST_CHECK(12 == list.get_num_owned_points());

}

void assertNeighborhood_p2(){
	QUICKGRID::QuickGridData gridData = getGrid();
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	BOOST_CHECK(12 == gridData.numPoints);
	BOOST_CHECK(12 == list.get_num_owned_points());

}

void assertNeighborhood_p3(){
	QUICKGRID::QuickGridData gridData = getGrid();
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	BOOST_CHECK(12 == gridData.numPoints);
	BOOST_CHECK(12 == list.get_num_owned_points());
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	if(0 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "ut_Y-Z_crack_jacobian_np4_p0" );
		proc->add(BOOST_TEST_CASE( &assertNeighborhood_p0 ));
		framework::master_test_suite().add( proc );
		return success;
	}
	if(1 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "ut_Y-Z_crack_jacobian_np4_p1" );
		proc->add(BOOST_TEST_CASE( &assertNeighborhood_p1 ));
		framework::master_test_suite().add( proc );
		return success;
	}
	if(2 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "ut_Y-Z_crack_jacobian_np4_p2" );
		proc->add(BOOST_TEST_CASE( &assertNeighborhood_p2 ));
		framework::master_test_suite().add( proc );
		return success;
	}
	if(3 == myRank){
		test_suite* proc = BOOST_TEST_SUITE( "ut_Y-Z_crack_jacobian_np4_p3" );
		proc->add(BOOST_TEST_CASE( &assertNeighborhood_p3 ));
		framework::master_test_suite().add( proc );
		return success;
	}
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
	// Initialize MPI and timer
	PdutMpiFixture myMpi = PdutMpiFixture(argc,argv);

	// These are static (file scope) variables
	myRank = myMpi.rank;
	numProcs = myMpi.numProcs;
	/**
	 * This test only make sense for numProcs == 4
	 */
	if(4 != numProcs){
		std::cerr << "Unit test runtime ERROR: ut_Y-Z_crack_jacobian_np4 only makes sense on 4 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 4 ./ut_Y-Z_crack_jacobian_np4" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}



	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
