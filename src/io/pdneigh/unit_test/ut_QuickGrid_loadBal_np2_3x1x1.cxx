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


#include "PdZoltan.h"
#include "QuickGrid.h"
#include "NeighborhoodList.h"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
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

const size_t nx = 3;
const size_t ny = 1;
const size_t nz = 1;
const double xStart = -1.5;
const double xLength = 3.0;
const double yStart = -0.5;
const double yLength = 1.0;
const double zStart = -0.5;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const double _cellVolume = xSpec.getCellSize()*ySpec.getCellSize()*zSpec.getCellSize();
const double x1 = xStart+xSpec.getCellSize()/2.0;
const double x2 = x1 + xSpec.getCellSize();
const double x3 = x2 + xSpec.getCellSize();
const double y = yStart + ySpec.getCellSize()/2.0;
const double z = zStart + zSpec.getCellSize()/2.0;
const size_t numCells = nx*ny*nz;

static int _neighborList[] = {
		1,1,                     /* numNeigh, neighbors */
		2,0,2,                   /* numNeigh, neighbors */
		1,1                      /* numNeigh, neighbors */
};
static int _neighborListSizeP0 = 5;
static int _neighborListSizeP1 = 2;

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


QUICKGRID::QuickGridData getGrid( int numProcs, int myRank) {
	double horizon = 1.1;

	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData gridData =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	gridData=PDNEIGH::getLoadBalancedDiscretization(gridData);
	return gridData;
}

TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np2_3x1x1, p0) {

       if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 2);

       if(numProcs != 2){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np2_3x1x1 only makes sense on 2 processors." << std::endl;
          return;
      }

      

    if(myRank == 0){
           
	QUICKGRID::QuickGridData gridData = getGrid( numProcs, myRank);

	TEST_ASSERT(0 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == gridData.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == gridData.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = gridData.numPoints;

	/*
	 * Zoltan load balances this such that 2 points end up on P0
	 */
	TEST_ASSERT(2 == myNumPoints);

	/*
	 * assert length of neighborhood list on this processor
	 */
	TEST_ASSERT( _neighborListSizeP0 == gridData.sizeNeighborhoodList );

	/*
	 * Assert global ids on this processor
	 * Assert neighborhood list
	 */
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int *neighborhoodList = gridData.neighborhood.get();
	int *_neighAns = _neighborList;
	int start = 0;
	for(size_t id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
		TEST_ASSERT( *gIdsPtr == (int)id );
		int numNeigh = *_neighAns;

		TEST_ASSERT( numNeigh == *neighborhoodList ); _neighAns++; neighborhoodList++;
		for(int i=0;i<numNeigh;i++){
			TEST_ASSERT( *_neighAns == *neighborhoodList ); _neighAns++; neighborhoodList++;
		}
		/*
		 * coordinates
		 */
	}

	/*
	 * Assert coordinates
	 */
	int id=0;
	const double tolerance = 1.0e-15;
	double *r = gridData.myX.get();
	TEST_FLOATING_EQUALITY(r[3*id+0],x1,tolerance);
	TEST_FLOATING_EQUALITY(r[3*id+1],y,tolerance);
	TEST_FLOATING_EQUALITY(r[3*id+2],z,tolerance);
	id=1;
	TEST_FLOATING_EQUALITY(r[3*id+0],x2,tolerance);
	TEST_FLOATING_EQUALITY(r[3*id+1],y,tolerance);
	TEST_FLOATING_EQUALITY(r[3*id+2],z,tolerance);

	// assert cell volumes
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,_cellVolume,tolerance);
	}

     
      }

}


TEUCHOS_UNIT_TEST(QuickGrid_loadBal_np2_3x1x1, p1) {


      if (!init) initialize();

       int numProcs = comm->NumProc();
       int myRank   = comm->MyPID();


       TEST_COMPARE(numProcs, ==, 2);

       if(numProcs != 2){
          std::cerr << "Unit test runtime ERROR: ut_QuickGrid_loadBal_np2_3x1x1 only makes sense on 2 processors." << std::endl;
         return;
      }



        // TEST_COMPARE(myRank, ==, 1);

     if(myRank == 1){
           
   
        QUICKGRID::QuickGridData gridData = getGrid( numProcs, myRank);

	TEST_ASSERT(1 == myRank);
	/*
	 * problem dimension is 3
	 */
	TEST_ASSERT(3 == gridData.dimension);

	/*
	 * Total number of cells in test
	 */
	TEST_ASSERT(nx*ny*nz == gridData.globalNumPoints);

	/*
	 * Number of cells on this processor
	 */
	int myNumPoints = gridData.numPoints;

	/*
	 * Zoltan load balances this such that 1 points end up on P1
	 */
	TEST_ASSERT(1 == myNumPoints);

	/*
	 * assert length of neighborhood list on this processor
	 */
	TEST_ASSERT( _neighborListSizeP1 == gridData.sizeNeighborhoodList );

	/*
	 * Assert global ids on this processor
	 * Assert neighborhood list
	 */
	shared_ptr<int> gIds = gridData.myGlobalIDs;
	int *gIdsPtr = gIds.get();
	int *neighborhoodList = gridData.neighborhood.get();
	int *_neighAns = &_neighborList[5];
	int start = 2;
	for(size_t id=start;id<gridData.numPoints+start;id++,gIdsPtr++){
		TEST_ASSERT( *gIdsPtr == (int)id );
		int numNeigh = *_neighAns;
		TEST_ASSERT( numNeigh == *neighborhoodList ); _neighAns++; neighborhoodList++;
		for(int i=0;i<numNeigh;i++){
			TEST_ASSERT( *_neighAns == *neighborhoodList ); _neighAns++; neighborhoodList++;
		}
	}

	/*
	 * Assert coordinates
	 */
	int id=0; // this is a local id
	const double tolerance = 1.0e-15;
	double *r = gridData.myX.get();
        TEST_FLOATING_EQUALITY(r[3*id+0],x3,tolerance);
	TEST_FLOATING_EQUALITY(r[3*id+1],y,tolerance);
	TEST_FLOATING_EQUALITY(r[3*id+2],z,tolerance);

	// assert cell volumes
	double *v = gridData.cellVolume.get();
	double *end = v+myNumPoints;
	for(; v != end ; v++){
		TEST_FLOATING_EQUALITY(*v,_cellVolume,tolerance);
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
