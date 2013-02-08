/*! \file Peridigm_ProximitySearch.cpp */
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

#include <Epetra_Import.h>
#include "Peridigm_ProximitySearch.hpp"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "pdneigh/PdZoltan.h"
#include "pdneigh/BondFilter.h"
#include "pdneigh/NeighborhoodList.h"

using namespace std;

void PeridigmNS::ProximitySearch::NeighborListToEpetraVector(int neighborListSize, const int* neighborList, Epetra_BlockMap& map, Teuchos::RCP<Epetra_Vector>& epetraVector)
{
  // Copy the neighbor list information into an Epetra_Vector with variable-length elements

  // Epetra_BlockMap and Epetra_Map do not allow elements of size zero
  // So, any point with zero neighbors must be excluded from the Epetra_Vector

  // The following are needed to construct the Epetra_BlockMap
  int numGlobalElements(-1);
  int numMyElements(0);
  vector<int> myGlobalElements;
  vector<int> elementSizeList;
  int indexBase(0);

  // Record the number of elements in the neighborList, which could be larger
  // than numMyElements in the event that there are points with zero neighbors

  int numMyElementsInNeighborList(0);

  // Determine the number of elements and resize the arrays
  int neighborListIndex(0), numNeighbors(0);
  while(neighborListIndex < neighborListSize){
    numNeighbors = neighborList[neighborListIndex++];
    neighborListIndex += numNeighbors;
    numMyElementsInNeighborList += 1;
    if(numNeighbors > 0)
      numMyElements += 1;
  }
  myGlobalElements.resize(numMyElements);
  elementSizeList.resize(numMyElements);

  // Fill the arrays
  neighborListIndex = 0;
  numNeighbors = 0;
  int numElementsWithZeroNeighbors(0);
  for(int i=0 ; i<numMyElementsInNeighborList ; ++i){
    numNeighbors = neighborList[neighborListIndex++];
    neighborListIndex += numNeighbors;
    if(numNeighbors > 0){
      myGlobalElements[i-numElementsWithZeroNeighbors] =  map.GID(i);
      elementSizeList[i-numElementsWithZeroNeighbors] = numNeighbors;
    }
    else{
      numElementsWithZeroNeighbors += 1;
    }
  }

  // Create the map
  Epetra_BlockMap variableElementSizeMap(numGlobalElements, numMyElements, &myGlobalElements[0], &elementSizeList[0], indexBase, map.Comm());

  // Create the vector
  epetraVector = Teuchos::rcp(new Epetra_Vector(variableElementSizeMap));

  // Fill the epetraVector with the neighbor list

}

void PeridigmNS::ProximitySearch::GlobalProximitySearch(Epetra_Vector& x,
                                                        double searchRadius,
                                                        std::vector<int>& neighborGlobalIdList)
{
  // Copy information from the Epetra_Vector into a QUICKGRID::Data object
  const Epetra_BlockMap& originalMap = x.Map();
  int dimension = 3;
  int numMyElements = originalMap.NumMyElements();
  int numGlobalElements = originalMap.NumGlobalElements();
  int* globalIds = originalMap.MyGlobalElements();
  double* coordinates;
  x.ExtractView(&coordinates);
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(numMyElements, dimension);
  decomp.globalNumPoints = numGlobalElements;
  int* decompGlobalIds = decomp.myGlobalIDs.get();
  double* decompVolumes = decomp.cellVolume.get();
  double* decompCoordinates = decomp.myX.get();
  for(int i=0 ; i<numMyElements ; ++i){
    decompGlobalIds[i] = globalIds[i];
    decompVolumes[i] = 1.0; // use dummy value
    decompCoordinates[i*3] = coordinates[i*3];
    decompCoordinates[i*3+1] = coordinates[i*3+1];
    decompCoordinates[i*3+2] = coordinates[i*3+2];
  }

  // Rebalance the decomp (RCB decomposition via Zoltan)
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
 
  // Execute neighbor search
  std::tr1::shared_ptr<PdBondFilter::BondFilter> bondFilterPtr(new PdBondFilter::BondFilterDefault(false));
  std::tr1::shared_ptr<const Epetra_Comm> commSp(&originalMap.Comm(), NonDeleter<const Epetra_Comm>());
  PDNEIGH::NeighborhoodList list(commSp,
                                 decomp.zoltanPtr.get(),
                                 decomp.numPoints,
                                 decomp.myGlobalIDs,
                                 decomp.myX,
                                 searchRadius,
                                 bondFilterPtr);

  // Perform parallel bookkeeping to get the neighbor list that is currently distributed
  // according to the decomposition of decomp into the decomposition of the original
  // Epetra_Vector x.

  int sizeNeighborList = list.get_size_neighborhood_list();
  int* neighborList = list.get_neighborhood().get();

  // // Copy the neighbor list (in the rebalanced decomposition) into an Epetra_Vector
  // Epetra_BlockMap neighborListRebalancedDecomp(decomp.globalNumPoints,
  //                                              decomp.numPoints,
  //                                              decomp.myGlobalIDs,
  //                                              elementSizeList,
  //                                              0,
  //                                              originalMap.Comm());

  // Create an Epetra_Vector in the original configuration that can recieve the neighbor list

  // Communicate the number of neighbors
  Epetra_BlockMap rebalancedOneDimensionalMap(decomp.globalNumPoints, decomp.numPoints, decomp.myGlobalIDs.get(), 1, 0, originalMap.Comm());
  Epetra_Vector numberOfNeighborsRebalancedDecomp(rebalancedOneDimensionalMap);
  double* numberOfNeighborsRebalancedDecompPtr;
  numberOfNeighborsRebalancedDecomp.ExtractView(&numberOfNeighborsRebalancedDecompPtr);
  int neighborListIndex(0), numNeighbors(0);
  for(int i=0 ; i<numberOfNeighborsRebalancedDecomp.MyLength() ; ++i){
    numNeighbors = neighborList[neighborListIndex];
    numberOfNeighborsRebalancedDecompPtr[i] = static_cast<double>(numNeighbors);
    neighborListIndex += numNeighbors;
  }
  int* myGlobalElementsOriginalMap = originalMap.MyGlobalElements();
  Epetra_BlockMap originalOneDimensionalMap(originalMap.NumGlobalElements(), originalMap.NumMyElements(), myGlobalElementsOriginalMap, 1, 0, originalMap.Comm());
  Epetra_Vector numberOfNeighborsOriginalDecomp(originalOneDimensionalMap);
  Epetra_Import oneDimensionalImporter(originalOneDimensionalMap, rebalancedOneDimensionalMap);
  numberOfNeighborsOriginalDecomp.Import(numberOfNeighborsRebalancedDecomp, oneDimensionalImporter, Insert);

  // Communicate the neighbor list
  // The fact that Epetra_BlockMaps do not allow zero-length elements adds some complexity here...
  // int numMyElementsUpperBound = rebalancedOneDimensionalMap->NumMyElements();
  // int numGlobalElements = -1; 
  // int numMyElements = 0;
  // int* rebalancedOneDimensionalMapGlobalElements = rebalancedOneDimensionalMap->MyGlobalElements();
  // int* myGlobalElements = new int[numMyElementsUpperBound];
  // int* elementSizeList = new int[numMyElementsUpperBound];
  // int numPointsWithZeroNeighbors = 0;
  // for(int i=0 ; i<numMyElementsUpperBound ; ++i){
  //   int numBonds = (int)( (*rebalancedNumberOfBonds)[i] );
  //   if(numBonds > 0){
  //     numMyElements++;
  //     myGlobalElements[i-numPointsWithZeroNeighbors] = rebalancedOneDimensionalMapGlobalElements[i];
  //     elementSizeList[i-numPointsWithZeroNeighbors] = numBonds;
  //   }
  //   else{
  //     numPointsWithZeroNeighbors++;
  //   }
  // }
  // int indexBase = 0;
  // Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap = 
  //   Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, *peridigmComm));
  // delete[] myGlobalElements;
  // delete[] elementSizeList;



  //   int globalID = oneDimensionalMap->GID(i);
  //   int bondMapLocalID = bondMap->LID(globalID);
  //   if(bondMapLocalID != -1)
  //     (*numberOfBonds)[i] = (double)( bondMap->ElementSize(i) );
  // }
  // Teuchos::RCP<Epetra_Vector> rebalancedNumberOfBonds = Teuchos::rcp(new Epetra_Vector(*rebalancedOneDimensionalMap));
  // rebalancedNumberOfBonds->Import(*numberOfBonds, *oneDimensionalMapToRebalancedOneDimensionalMapImporter, Insert);


  // Create a list of neighbors for

  // decomp.neighborhood=list.get_neighborhood();
  // decomp.sizeNeighborhoodList=list.get_size_neighborhood_list();
  // decomp.neighborhoodPtr=list.get_neighborhood_ptr();

}

