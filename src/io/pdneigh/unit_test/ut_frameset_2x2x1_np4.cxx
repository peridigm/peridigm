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
<<<<<<< HEAD
=======
#include "Teuchos_GlobalMPISession.hpp"
>>>>>>> FETCH_HEAD
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "../BondFilter.h"

#include "zoltan.h"
#include "Sortable.h"
#include "mpi.h"
#include <valarray>
#include <iostream>
#include <cmath>
#include <map>
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


using std::tr1::shared_ptr;
using std::cout;
using std::set;
using std::map;


using namespace QUICKGRID;
using namespace PDNEIGH;

const int nx = 2;
const int ny = nx;
const double lX = 2.0;
const double lY = lX;
const double lZ = 1.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const int nz = 1;
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;
const double horizon=1.1;
using std::cout;
using std::endl;



QUICKGRID::QuickGridData getGrid( int numProcs, int myRank) {
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	// This load-balances
	decomp = getLoadBalancedDiscretization(decomp);
	return decomp;
}


shared_ptr< std::set<int> > constructFrame(PDNEIGH::NeighborhoodList& list) {
	shared_ptr< std::set<int> > frameSet = UTILITIES::constructFrameSet(list.get_num_owned_points(),list.get_owned_x(),list.get_frameset_buffer_size());
	/*
	 * There is only 1 point on each processor and by default it must show up in the frameset
	 */
	
	return frameSet;
}

TEUCHOS_UNIT_TEST(Frameset_2x2x1_np4, CreateNeighborhoodTest) {

<<<<<<< HEAD
	QUICKGRID::QuickGridData decomp = getGrid();
=======
        

       
        shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));

        int numProcs = comm->NumProc();
        int myRank = comm->MyPID();

        TEST_COMPARE(numProcs, ==, 4);

        if(numProcs != 4){
           std::cerr << "Unit test runtime ERROR: ut_frameset_2x2x1_np4 only makes sense on 4 processors." << std::endl;
           return;
        }
       



	QUICKGRID::QuickGridData decomp = getGrid( numProcs, myRank);
        
>>>>>>> FETCH_HEAD
	TEST_ASSERT(1==decomp.numPoints);
	TEST_ASSERT(4==decomp.globalNumPoints);
	shared_ptr<BondFilter> bondFilterPtr(new PdBondFilter::BondFilterDefault(true));
	//shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,2.0*horizon,bondFilterPtr);
	shared_ptr< std::set<int> > frameSet = constructFrame(list);
        TEST_ASSERT(1==frameSet->size());
	TEST_ASSERT(5==list.get_size_neighborhood_list());
<<<<<<< HEAD
=======

       
>>>>>>> FETCH_HEAD
	/*
	 * Neighborhood of every point should have every other point
	 */
	int n[] = {0,1,2,3};
	set<int> neighSet(n,n+4);
	int *neighborhood = list.get_neighborhood().get();
	int numNeigh = *neighborhood;
	TEST_ASSERT(4==numNeigh);
	neighborhood++;
	for(int i=0;i<numNeigh;i++,neighborhood++){
		TEST_ASSERT(neighSet.end()!=neighSet.find(*neighborhood));
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
