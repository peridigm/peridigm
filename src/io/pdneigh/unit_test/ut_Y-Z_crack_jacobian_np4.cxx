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
#include "Teuchos_GlobalMPISession.hpp"
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "../BondFilter.h"

#include "PdutMpiFixture.h"
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

using std::cout;
using std::endl;
using std::set;


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

bool init = false;
Teuchos::RCP<Epetra_Comm> comm;   


void initialize(){
         
       #ifdef HAVE_MPI
              
              comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
       #else
              comm = rcp(new Epetra_SerialComm);
       #endif
     
       
       init = true;
}

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




QUICKGRID::QuickGridData getGrid(int numProcs, int myRank) {

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

TEUCHOS_UNIT_TEST(Y_Z_crack_jacobian_np4, AssertNeighborhood_p0Test) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 4);

       if(numProcs != 4){
          std::cerr << "Unit test runtime ERROR:  ut_Y-Z_crack_jacobian_np4 only makes sense on 4 processors." << std::endl;
          return;
      }

      

       if(myRank == 0){

	QUICKGRID::QuickGridData gridData = getGrid( numProcs, myRank);
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
	TEST_ASSERT(12 == gridData.numPoints);
	TEST_ASSERT(12 == list.get_num_owned_points());
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

		TEST_ASSERT(GIDS[j]==*gids);
		int numNeighAnswer = N[j];
		int numNeigh = *neigh; neigh++;
//		cout << "rank, gid, numNeigh = " << myRank << ", " << *gids << ", " << numNeigh << endl;
//		printNeighborhood(numNeigh,neigh);
		TEST_ASSERT(numNeighAnswer==numNeigh);
//		int *neighAnswer = NN[j];
		set<int> neighAnswer(NN[j],NN[j]+numNeigh);
		set<int>::iterator end = neighAnswer.end();
		for(int n=0;n<numNeigh;n++){
			TEST_ASSERT(end != neighAnswer.find(*(neigh+n)));
		}

		/*
		 * Move to next point
		 */
		neigh = neigh+numNeigh;
	}

    }

//{
//		int gid=7;
//		int numNeigh = *neigh; neigh++;
//		cout << "gid, numNeigh = " << gid << ", " << numNeigh << endl;
//		printNeighborhood(numNeigh,neigh);
//
//	}


}


TEUCHOS_UNIT_TEST(Y_Z_crack_jacobian_np4, AssertNeighborhood_p1Test) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 4);

       if(numProcs != 4){
          std::cerr << "Unit test runtime ERROR:  ut_Y-Z_crack_jacobian_np4 only makes sense on 4 processors." << std::endl;
          return;
      }

      

       if(myRank == 1){

	QUICKGRID::QuickGridData gridData = getGrid(numProcs, myRank);
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	TEST_ASSERT(12 == gridData.numPoints);
	TEST_ASSERT(12 == list.get_num_owned_points());

      }

}

TEUCHOS_UNIT_TEST(Y_Z_crack_jacobian_np4, AssertNeighborhood_p2Test) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 4);

       if(numProcs != 4){
          std::cerr << "Unit test runtime ERROR:  ut_Y-Z_crack_jacobian_np4 only makes sense on 4 processors." << std::endl;
          return;
      }

      

       if(myRank == 2){


	QUICKGRID::QuickGridData gridData = getGrid( numProcs, myRank);
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	TEST_ASSERT(12 == gridData.numPoints);
	TEST_ASSERT(12 == list.get_num_owned_points());

     }

}


TEUCHOS_UNIT_TEST(Y_Z_crack_jacobian_np4, AssertNeighborhood_p3Test) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 4);

       if(numProcs != 4){
          std::cerr << "Unit test runtime ERROR:  ut_Y-Z_crack_jacobian_np4 only makes sense on 4 processors." << std::endl;
          return;
      }

      

       if(myRank == 3){

	QUICKGRID::QuickGridData gridData = getGrid( numProcs, myRank);
	FinitePlane crackPlane=getYZ_CrackPlane();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * There are a total of 48 points = nx * ny * nz = 4 * 4 * 3
	 * Because of this mesh, each processor should get 12 points
	 */
	TEST_ASSERT(12 == gridData.numPoints);
	TEST_ASSERT(12 == list.get_num_owned_points());
      }
}



int main
(
		int argc,
		char* argv[]
)
{
	Teuchos::GlobalMPISession mpiSession(&argc, &argv);


	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
