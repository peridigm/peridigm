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
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include <iostream>


using std::tr1::shared_ptr;


QUICKGRID::QuickGridData setUp(){
	/*
	 * This setup is identical to the setup in utPdQuickGridHorizon.cxx
	 *     function call "CellsPerProcessor3D_smallNeighborhoodSerialTest_NumProcs_1"
	 * That function tests the setup rigorously
	 *
	 */

	// use this spec along the x and y axes
	size_t numCells = 3;
	double xStart = 1.0;
	double xLength=1.0;
	QUICKGRID::Spec1D spec(numCells,xStart,xLength);

	// this creates a different mesh along the z-axis
	size_t numCellsZ=2;
	QUICKGRID::Spec1D zSpec(numCellsZ,xStart,xLength);

	// This scale factor pushes the horizon just over the line so that 1 cell is included
	//  in the half neighborhood
	double SCALE=1.0;
	double horizon = SCALE*spec.getCellSize();
	int numProc=1;

	QUICKGRID::TensorProduct3DMeshGenerator cellIter(numProc,horizon,spec,spec,zSpec);
	QUICKGRID::QuickGridData pdGridDataProc0 = cellIter.allocatePdGridData();
	std::pair<QUICKGRID::Cell3D,QUICKGRID::QuickGridData> pdGridData = cellIter.beginIterateProcs(pdGridDataProc0);
	QUICKGRID::QuickGridData gridData = pdGridData.second;

	return gridData;


}


TEUCHOS_UNIT_TEST(ZoltanQuery_pointSizeInBytes, SmallNeighborhood) {

	/* This test is about exercising the zoltan call back
	 * function:
	 * void zoltanQuery_pointSizeInBytes
(
		void *pdGridData,
		int numGids,
		int numLids,
		int numPoints,
		ZOLTAN_ID_PTR zoltanGlobalIds,
		ZOLTAN_ID_PTR zoltanLocalIds,
		int *sizes,
		int *ierr
)
	 *
	 */
	QUICKGRID::QuickGridData gridData = setUp();
	int numGids = 1;
	int numLids = 1;
	int numPoints = 9;

	/*
	 * The above setup has 2 slabs of 9 nodes each
	 * Use the call function as if we wanted to ship out 1 of the slabs
	 *
	 */
	ZOLTAN_ID_TYPE exportLocalIds[] = {9, 10, 11, 12, 13, 14, 15, 16, 17};
	// Global ids are not used in the size function
	ZOLTAN_ID_TYPE *exportGlobalIds = NULL;
	int sizes[] = {0,0,0,0,0,0,0,0,0};
	int ierr[] = {0};
	PDNEIGH::zoltanQuery_pointSizeInBytes(&gridData,numGids,numLids,numPoints,exportGlobalIds,exportLocalIds,sizes,ierr);

	// assert sizes
	// these are the known number of neighbors for each point being exported
	int numNeighbors[]={ 7,11,7,11,17,11,7,11,7};
	               // coordinates + volume + numNeigh + neighbors
	int sizeAnswers[] = {0,0,0,0,0,0,0,0,0};
	int ptrAnswers[] = {0,0,0,0,0,0,0,0,0};
	int dimension=3;
	int ptrStart = 4*7 + 4 + 4*11 + 4 + 17*1 + 1; // = 98
	ptrAnswers[0] = ptrStart;
	for(int p=1;p<numPoints;p++){
		ptrAnswers[p] = ptrAnswers[p-1] + numNeighbors[p-1] + 1;
	}

	for(int p=0;p<numPoints;p++){
		sizeAnswers[p]=(dimension+1)*sizeof(double)+(1+numNeighbors[p])*sizeof(int);
	}
	int *neighPtr = gridData.neighborhoodPtr.get();
	for(int p=0;p<9;p++){
		TEST_ASSERT( ptrAnswers[p] == neighPtr[exportLocalIds[p]] );
		TEST_ASSERT( sizeAnswers[p] == sizes[p] );
	}
	std::cout << std::endl;
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

