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
#include <Epetra_SerialComm.h>
#include <Teuchos_RCP.hpp>
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "quick_grid/QuickGridData.h"
#include "../NeighborhoodList.h"
#include "../BondFilter.h"
#include "../../Peridigm_ZoltanSearchTree.hpp"
#include <iostream>

using namespace PdBondFilter;
using std::tr1::shared_ptr;
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
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double SCALE=1.1*sqrt(3);
const double horizon = SCALE*xSpec.getCellSize();

// known local ids
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

template<class T>
struct NonDeleter{
	void operator()(T* d) {}
};

QUICKGRID::QuickGridData getGrid() {

	/*
	 * 2x2x2 Grid of points/cells
	 */
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

	return decomp;
}

FinitePlane getCase_1a(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=-1.0/sqrt2;n[1]=1.0/sqrt2;n[2]=0.0;
	double r0[3]; r0[0]=0.0; r0[1]=0.0; r0[2]=0.0;
	double ub[3]; ub[0]=1.0/sqrt2; ub[1]=1.0/sqrt2;ub[2]=0.0;
	double b=sqrt(2), a=1.0;
	return FinitePlane(n,r0,ub,b,a);
}

FinitePlane getCase_1b(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=1.0/sqrt2;n[1]=1.0/sqrt2;n[2]=0.0;
	double r0[3]; r0[0]=0.0; r0[1]=1.0; r0[2]=0.0;
	double ub[3]; ub[0]=1.0/sqrt2; ub[1]=-1.0/sqrt2;ub[2]=0.0;
	double b=sqrt(2), a=1.0;
	return FinitePlane(n,r0,ub,b,a);
}

FinitePlane getCase_2a(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=-1.0/sqrt2;n[1]=0.0; n[2]=1.0/sqrt2;
	double r0[3]; r0[0]=0.0; r0[1]=0.0; r0[2]=0.0;
	double ub[3]; ub[0]=0.0; ub[1]=1.0;ub[2]=0.0;
	double b=1.0, a=sqrt2;
	return FinitePlane(n,r0,ub,b,a);
}

FinitePlane getCase_2b(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=1.0/sqrt2;n[1]=0.0; n[2]=1.0/sqrt2;
	double r0[3]; r0[0]=1.0; r0[1]=1.0; r0[2]=0.0;
	double ub[3]; ub[0]=0.0; ub[1]=-1.0;ub[2]=0.0;
	double b=1.0, a=sqrt2;
	return FinitePlane(n,r0,ub,b,a);
}

FinitePlane getCase_3a(){
	double sqrt2=sqrt(2.0);
	double n[3]; n[0]=0.0; n[1]=-1.0/sqrt2; n[2]=1.0/sqrt2;
	double r0[3]; r0[0]=0.0; r0[1]=0.0; r0[2]=0.0;
	double ub[3]; ub[0]=0.0; ub[1]=1.0/sqrt2;ub[2]=1.0/sqrt2;
	double a=1.0, b=sqrt2;
	return FinitePlane(n,r0,ub,b,a);
}


TEUCHOS_UNIT_TEST(FinitePlaneFilter, Case_1a) {

       
	QUICKGRID::QuickGridData decomp = getGrid();
        int numOverlapPoints = decomp.numPoints;

        /*
         * SANITY CHECK on Expected coordinates and IDs
         */

        TEST_ASSERT(8==numOverlapPoints);
		for(int j=0;j<numOverlapPoints;j++){
			TEST_ASSERT(decomp.myGlobalIDs.get()[j]==ids[j]);
			TEST_ASSERT(decomp.myX.get()[j*3+0]==x[j*3+0]);
			TEST_ASSERT(decomp.myX.get()[j*3+1]==x[j*3+1]);
			TEST_ASSERT(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}

	FinitePlane plane = getCase_1a();
	shared_ptr<BondFilter> filterPtr=shared_ptr<BondFilter>(new FinitePlaneFilter(plane));

	/*
	 * Create KdTree; Since this is serial xOwned = xOverlap and numOwned = numOverlap
	 */
	double* xOverlapPtr = decomp.myX.get();
	

    PeridigmNS::ZoltanSearchTree searchTree(numOverlapPoints, xOverlapPtr);

	/*
	 * ANSWERS for each ID
	 * list size for each point
	 */
	
	bool markForExclusion[8];
	// Expected: filter should evaluate this list size for each id
	int size[] = {4,2,2,4,4,2,2,4};
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
		for(int i=0;i<8;i++){
			/*
			 * look at neighborhood of id = 0
			 */
			int id=ids[i];
            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point *
			 */
			double *x = decomp.myX.get()+3*id;
			searchTree.FindPointsWithinRadius(x, horizon, treeList);

			/*
			 * Now determine which points are included
			 */
			filterPtr->filterBonds(treeList,x,id,decomp.myX.get(),markForExclusion);
			bool *flags = expectedFlags[i];
			/*
			 * Assert flags
			 */
			for(int j=0;j<8;j++){
//				cout << "filter flag, expected flag = " << *(markForExclusion+j) << ", " << *(flags+j) << endl;
				TEST_ASSERT(*(flags+j)==*(markForExclusion+j));
			}
		}
	}

    Epetra_SerialComm serialComm;
    std::tr1::shared_ptr<const Epetra_Comm> comm(&serialComm,NonDeleter<const Epetra_Comm>());
	PDNEIGH::NeighborhoodList neighList(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);

	/*
	 * Assert neighbors
	 */
	int *neigh = neighList.get_neighborhood().get();
	for(int n=0;n<8;n++){
		/*
		 * Assert number of neighbors
		 */
		int numNeigh = *neigh; neigh++;
		/*
		 * Note that we subtract 1 here because the list
		 * size is always '1' greater than the number of
		 * neighbors in order to store 'num neighbors'
		 * in the list.
		 */
		TEST_ASSERT((static_cast<int>(size[n])-1)==numNeigh);

		/*
		 * Expected
		 */
		bool *flags = expectedFlags[n];
		for(int j=0;j<numNeigh;j++,neigh++){
          int id = *neigh;
          TEST_ASSERT((flags+id));
		}

	}


}


TEUCHOS_UNIT_TEST(FinitePlaneFilter, Case_1b) {

  QUICKGRID::QuickGridData decomp = getGrid();

  int numOverlapPoints = decomp.numPoints;

        /*
         * SANITY CHECK on Expected coordinates and IDs
         */

        TEST_ASSERT(8==numOverlapPoints);
		for(int j=0;j<numOverlapPoints;j++){
			TEST_ASSERT(decomp.myGlobalIDs.get()[j]==ids[j]);
			TEST_ASSERT(decomp.myX.get()[j*3+0]==x[j*3+0]);
			TEST_ASSERT(decomp.myX.get()[j*3+1]==x[j*3+1]);
			TEST_ASSERT(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}
	FinitePlane plane = getCase_1b();
        Teuchos::RCP<BondFilter> filterPtr=Teuchos::rcp<BondFilter>(new FinitePlaneFilter(plane));

	/*
	 * Create KdTree; Since this is serial xOwned = xOverlap and numOwned = numOverlap
	 */
	double* xOverlapPtr = decomp.myX.get();
	

    PeridigmNS::ZoltanSearchTree searchTree(numOverlapPoints, xOverlapPtr);

	/*
	 * ANSWERS for each ID
	 * list size for each point
	 */
	
	bool markForExclusion[8];
	// Expected: filter should evaluate this list size for each id
	int size[] = {2,4,4,2,2,4,4,2};
	// Expected: filter should return these flags for each local id
	bool n0[]={1,1,1,1,0,1,1,1};
	bool n1[]={1,1,0,1,1,0,0,1};
	bool n2[]={1,0,1,1,1,0,0,1};
	bool n3[]={1,1,1,1,1,1,1,0};
	bool n4[]={0,1,1,1,1,1,1,1};
	bool n5[]={1,0,0,1,1,1,0,1};
	bool n6[]={1,0,0,1,1,0,1,1};
	bool n7[]={1,1,1,0,1,1,1,1};
	bool * expectedFlags[] = {n0,n1,n2,n3,n4,n5,n6,n7};
	int numCheck=8;
	{
		for(int i=0;i<numCheck;i++){
			/*
			 * look at neighborhood of id = 0
			 */
			int id=ids[i];
            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point *
			 */
			double *x = decomp.myX.get()+3*id;
			searchTree.FindPointsWithinRadius(x, horizon, treeList);

			/*
			 * Now determine which points are included
			 */
			filterPtr->filterBonds(treeList,x,id,decomp.myX.get(),markForExclusion);
			bool *flags = expectedFlags[i];
			/*
			 * Assert flags
			 */
			for(int j=0;j<8;j++){
//				cout << "filter flag, expected flag = " << *(markForExclusion+j) << ", " << *(flags+j) << endl;
				TEST_ASSERT(*(flags+j)==*(markForExclusion+j));
			}
		}
	}

    Epetra_SerialComm serialComm;
    std::tr1::shared_ptr<const Epetra_Comm> comm(&serialComm,NonDeleter<const Epetra_Comm>());
	PDNEIGH::NeighborhoodList neighList(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);

	/*
	 * Assert neighbors
	 */
	int *neigh = neighList.get_neighborhood().get();
	for(int n=0;n<numCheck;n++){
		/*
		 * Assert number of neighbors
		 */
		int numNeigh = *neigh; neigh++;
		/*
		 * Note that we subtract 1 here because the list
		 * size is always '1' greater than the number of
		 * neighbors in order to store 'num neighbors'
		 * in the list.
		 */
		TEST_ASSERT((size[n]-1)==numNeigh);

		/*
		 * Expected
		 */
		bool *flags = expectedFlags[n];
		for(int j=0;j<numNeigh;j++,neigh++){
			int id = *neigh;
			TEST_ASSERT((flags+id));
		}

	}


}

TEUCHOS_UNIT_TEST(FinitePlaneFilter, Case_2a) {

  QUICKGRID::QuickGridData decomp = getGrid();
  int numOverlapPoints = decomp.numPoints;

        /*
         * SANITY CHECK on Expected coordinates and IDs
         */

        TEST_ASSERT(8==numOverlapPoints);
		for(int j=0;j<numOverlapPoints;j++){
			TEST_ASSERT(decomp.myGlobalIDs.get()[j]==ids[j]);
			TEST_ASSERT(decomp.myX.get()[j*3+0]==x[j*3+0]);
			TEST_ASSERT(decomp.myX.get()[j*3+1]==x[j*3+1]);
			TEST_ASSERT(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}
	FinitePlane plane = getCase_2a();
    Teuchos::RCP<BondFilter> filterPtr=Teuchos::rcp<BondFilter>(new FinitePlaneFilter(plane));

	/*
	 * Create KdTree; Since this is serial xOwned = xOverlap and numOwned = numOverlap
	 */
	double* xOverlapPtr = decomp.myX.get();
	

    PeridigmNS::ZoltanSearchTree searchTree(numOverlapPoints, xOverlapPtr);

	/*
	 * ANSWERS for each ID
	 * list size for each point
	 */
	// known local ids
	
	bool markForExclusion[8];
	// Expected: filter should evaluate this list size for each id
	int size[] = {4,2,4,2,2,4,2,4};
	// Expected: filter should return these flags for each local id
	bool n0[]={1,1,0,1,1,0,1,0};
	bool n1[]={1,1,1,0,1,1,1,1};
	bool n2[]={0,1,1,1,1,0,1,0};
	bool n3[]={1,0,1,1,1,1,1,1};
	bool n4[]={1,1,1,1,1,1,0,1};
	bool n5[]={0,1,0,1,1,1,1,0};
	bool n6[]={1,1,1,1,0,1,1,1};
	bool n7[]={0,1,0,1,1,0,1,1};
	bool * expectedFlags[] = {n0,n1,n2,n3,n4,n5,n6,n7};
	{
		for(int i=0;i<8;i++){
			/*
			 * look at neighborhood of id = 0
			 */
			int id=ids[i];
            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point *
			 */
			double *x = decomp.myX.get()+3*id;
			searchTree.FindPointsWithinRadius(x, horizon, treeList);

			/*
			 * Now determine which points are included
			 */
			filterPtr->filterBonds(treeList,x,id,decomp.myX.get(),markForExclusion);
			bool *flags = expectedFlags[i];
			/*
			 * Assert flags
			 */
			for(int j=0;j<8;j++){
//				cout << "filter flag, expected flag = " << *(markForExclusion+j) << ", " << *(flags+j) << endl;
				TEST_ASSERT(*(flags+j)==*(markForExclusion+j));
			}
		}
	}

    Epetra_SerialComm serialComm;
    std::tr1::shared_ptr<const Epetra_Comm> comm(&serialComm,NonDeleter<const Epetra_Comm>());
	PDNEIGH::NeighborhoodList neighList(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);

	/*
	 * Assert neighbors
	 */
	int *neigh = neighList.get_neighborhood().get();
	for(int n=0;n<8;n++){
		/*
		 * Assert number of neighbors
		 */
		int numNeigh = *neigh; neigh++;
		/*
		 * Note that we subtract 1 here because the list
		 * size is always '1' greater than the number of
		 * neighbors in order to store 'num neighbors'
		 * in the list.
		 */
		TEST_ASSERT((size[n]-1)==numNeigh);

		/*
		 * Expected
		 */
		bool *flags = expectedFlags[n];
		for(int j=0;j<numNeigh;j++,neigh++){
			int id = *neigh;
			TEST_ASSERT((flags+id));
		}

	}

}

TEUCHOS_UNIT_TEST(FinitePlaneFilter, Case_2b) {


	QUICKGRID::QuickGridData decomp = getGrid();

         int numOverlapPoints = decomp.numPoints;

        /*
         * SANITY CHECK on Expected coordinates and IDs
         */

        TEST_ASSERT(8==numOverlapPoints);
		for(int j=0;j<numOverlapPoints;j++){
		        TEST_ASSERT(decomp.myGlobalIDs.get()[j]==ids[j]);
			TEST_ASSERT(decomp.myX.get()[j*3+0]==x[j*3+0]);
			TEST_ASSERT(decomp.myX.get()[j*3+1]==x[j*3+1]);
			TEST_ASSERT(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}
	FinitePlane plane = getCase_2b();
    Teuchos::RCP<BondFilter> filterPtr=Teuchos::rcp<BondFilter>(new FinitePlaneFilter(plane));

	/*
	 * Create KdTree; Since this is serial xOwned = xOverlap and numOwned = numOverlap
	 */
	double* xOverlapPtr = decomp.myX.get();

    PeridigmNS::ZoltanSearchTree searchTree(numOverlapPoints, xOverlapPtr);

	/*
	 * ANSWERS for each ID
	 * list size for each point
	 */
	
	bool markForExclusion[8];
	// Expected: filter should evaluate this list size for each id
	int size[] = {2,4,2,4,4,2,4,2};
	// Expected: filter should return these flags for each local id
	bool n0[]={1,1,0,1,1,1,1,1};
	bool n1[]={1,1,1,0,0,1,0,1};
	bool n2[]={0,1,1,1,1,1,1,1};
	bool n3[]={1,0,1,1,0,1,0,1};
	bool n4[]={1,0,1,0,1,1,0,1};
	bool n5[]={1,1,1,1,1,1,1,0};
	bool n6[]={1,0,1,0,0,1,1,1};
	bool n7[]={1,1,1,1,1,0,1,1};
	bool * expectedFlags[] = {n0,n1,n2,n3,n4,n5,n6,n7};
	int numCheck=8;
	{
		for(int i=0;i<numCheck;i++){
			/*
			 * look at neighborhood of id = 0
			 */
			int id=ids[i];
            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point *
			 */
			double *x = decomp.myX.get()+3*id;
			searchTree.FindPointsWithinRadius(x, horizon, treeList);

			/*
			 * Now determine which points are included
			 */
			filterPtr->filterBonds(treeList,x,id,decomp.myX.get(),markForExclusion);
			bool *flags = expectedFlags[i];
			/*
			 * Assert flags
			 */
			for(int j=0;j<8;j++){
//				cout << "filter flag, expected flag = " << *(markForExclusion+j) << ", " << *(flags+j) << endl;
				TEST_ASSERT(*(flags+j)==*(markForExclusion+j));
			}
		}
	}

    Epetra_SerialComm serialComm;
    std::tr1::shared_ptr<const Epetra_Comm> comm(&serialComm,NonDeleter<const Epetra_Comm>());
	PDNEIGH::NeighborhoodList neighList(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);

	/*
	 * Assert neighbors
	 */
	int *neigh = neighList.get_neighborhood().get();
	for(int n=0;n<numCheck;n++){
		/*
		 * Assert number of neighbors
		 */
		int numNeigh = *neigh; neigh++;
		/*
		 * Note that we subtract 1 here because the list
		 * size is always '1' greater than the number of
		 * neighbors in order to store 'num neighbors'
		 * in the list.
		 */
		TEST_ASSERT((size[n]-1)==numNeigh);

		/*
		 * Expected
		 */
		bool *flags = expectedFlags[n];
		for(int j=0;j<numNeigh;j++,neigh++){
			int id = *neigh;
			TEST_ASSERT((flags+id));
		}

	}


}

TEUCHOS_UNIT_TEST(FinitePlaneFilter, Case_3a) {

	QUICKGRID::QuickGridData decomp = getGrid();

         int numOverlapPoints = decomp.numPoints;

        /*
         * SANITY CHECK on Expected coordinates and IDs
         */

        TEST_ASSERT(8==numOverlapPoints);
		for(int j=0;j<numOverlapPoints;j++){
			TEST_ASSERT(decomp.myGlobalIDs.get()[j]==ids[j]);
			TEST_ASSERT(decomp.myX.get()[j*3+0]==x[j*3+0]);
			TEST_ASSERT(decomp.myX.get()[j*3+1]==x[j*3+1]);
			TEST_ASSERT(decomp.myX.get()[j*3+2]==x[j*3+2]);
		}
	FinitePlane plane = getCase_3a();
    Teuchos::RCP<BondFilter> filterPtr=Teuchos::rcp<BondFilter>(new FinitePlaneFilter(plane));

	/*
	 * Create KdTree; Since this is serial xOwned = xOverlap and numOwned = numOverlap
	 */
	double* xOverlapPtr = decomp.myX.get();
	

    PeridigmNS::ZoltanSearchTree searchTree(numOverlapPoints, xOverlapPtr);

	/*
	 * ANSWERS for each ID
	 * list size for each point
	 */
	
	bool markForExclusion[8];
	// Expected: filter should evaluate this list size for each id
	int size[] = {4,4,2,2,2,2,4,4};
	// Expected: filter should return these flags for each local id
	bool n0[]={1,0,1,1,1,1,0,0};
	bool n1[]={0,1,1,1,1,1,0,0};
	bool n2[]={1,1,1,0,1,1,1,1};
	bool n3[]={1,1,0,1,1,1,1,1};
	bool n4[]={1,1,1,1,1,0,1,1};
	bool n5[]={1,1,1,1,0,1,1,1};
	bool n6[]={0,0,1,1,1,1,1,0};
	bool n7[]={0,0,1,1,1,1,0,1};
	bool * expectedFlags[] = {n0,n1,n2,n3,n4,n5,n6,n7};
	{
		for(int i=0;i<8;i++){
			/*
			 * look at neighborhood of id = 0
			 */
			int id=ids[i];
            std::vector<int> treeList;
			/*
			 * Note that list returned includes this point *
			 */
			double *x = decomp.myX.get()+3*id;
			searchTree.FindPointsWithinRadius(x, horizon, treeList);

			/*
			 * Now determine which points are included
			 */
			filterPtr->filterBonds(treeList,x,id,decomp.myX.get(),markForExclusion);
			bool *flags = expectedFlags[i];
			/*
			 * Assert flags
			 */
			for(int j=0;j<8;j++){
//				cout << "filter flag, expected flag = " << *(markForExclusion+j) << ", " << *(flags+j) << endl;
				TEST_ASSERT(*(flags+j)==*(markForExclusion+j));
			}
		}
	}

    Epetra_SerialComm serialComm;
    std::tr1::shared_ptr<const Epetra_Comm> comm(&serialComm,NonDeleter<const Epetra_Comm>());
	PDNEIGH::NeighborhoodList neighList(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);

	/*
	 * Assert neighbors
	 */
	int *neigh = neighList.get_neighborhood().get();
	for(int n=0;n<8;n++){
		/*
		 * Assert number of neighbors
		 */
		int numNeigh = *neigh; neigh++;
		/*
		 * Note that we subtract 1 here because the list
		 * size is always '1' greater than the number of
		 * neighbors in order to store 'num neighbors'
		 * in the list.
		 */
		TEST_ASSERT((size[n]-1)==numNeigh);

		/*
		 * Expected
		 */
		bool *flags = expectedFlags[n];
		for(int j=0;j<numNeigh;j++,neigh++){
			int id = *neigh;
			TEST_ASSERT((flags+id));
		}

	}


}



int main
(
		int argc,
		char* argv[]
)
{


	// Initialize UTF
	 return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
