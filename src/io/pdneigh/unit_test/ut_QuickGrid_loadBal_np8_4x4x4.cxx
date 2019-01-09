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

#include <iostream>
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif




using std::shared_ptr;
using std::cout;


const int nx = 4;
const int ny = 4;
const int nz = 4;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = 0.0;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;

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



QUICKGRID::QuickGridData getGrid( int numProcs, int myRank ) {
	// This scale factor pushes the horizon just over the line so that 1 cells are included
	//  in the half neighborhood
	double SCALE=1.4;
	double horizon = SCALE*xSpec.getCellSize();
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	// This load-balances above mesh
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}

TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p0) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 0){
          

        QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(0 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=0 end should be on this processor
	 */
	int ids[] = {0,1,4,5,16,17,20,21};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P0 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
//	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

    }
}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p1) {

        
       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 1){



	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(1 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {8,9,12,13,24,25,28,29};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	//std::cout << "P1 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
	//	std::cout << gIdsPtr[p] << ", ";
	}
	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}
    
    }

}

TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p2) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 2){  




	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(2 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {34,35,38,39,50,51,54,55};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P2 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
	//	std::cout << gIdsPtr[p] << ", ";
	}
	//std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

     }

}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p3) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 3){

       

	QUICKGRID::QuickGridData decomp = getGrid(numProcs, myRank);
	TEST_ASSERT(3 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = { 26, 27, 30, 31, 10, 11, 14, 15};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P3 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
	//std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

     }

}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p4) {


       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
       }

      

        if(myRank == 4){
        

	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(4 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {2, 3, 6, 7, 18, 19, 22, 23};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P4 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
	//std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

    }

}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p5) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 5){

       

	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(5 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {40, 41, 44, 45, 56, 57, 60, 61};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P5 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
//	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

     }

}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p6) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 6){

	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(6 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {48, 49, 52, 53, 32, 33, 36, 37};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P6 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
//	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}

     }

}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np8_4x4x4, p7) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 8);

       if(numProcs != 8){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np8_4x4x4 only makes sense on 8 processors." << std::endl;
          return;
      }

      

        if(myRank == 7){

	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
	TEST_ASSERT(7 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == decomp.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == decomp.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = decomp.numPoints;
	TEST_ASSERT(8 == myNumPoints);

	/*
	 * Assert global ids on this processor; All points at x=L end should be on this processor
	 */
	int ids[] = {58, 59, 62, 63, 42, 43, 46, 47};
	shared_ptr<int> gIds = decomp.myGlobalIDs;
	int *gIdsPtr = gIds.get();
//	std::cout << "P6 gids = ";
	for(int p=0;p<myNumPoints;p++){
		TEST_ASSERT(gIdsPtr[p]==ids[p]);
//		std::cout << gIdsPtr[p] << ", ";
	}
//	std::cout << std::endl;

	// assert cell volumes
	double cellVolume = xSpec.getCellSize() * ySpec.getCellSize() * zSpec.getCellSize();
	const double tolerance = 1.0e-15;
	double *v = decomp.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,cellVolume,tolerance);
	}
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


