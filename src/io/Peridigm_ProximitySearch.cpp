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
#include "Peridigm_Discretization.hpp"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "pdneigh/PdZoltan.h"
#include "pdneigh/BondFilter.h"
#include "pdneigh/NeighborhoodList.h"

using namespace std;

void PeridigmNS::ProximitySearch::RebalanceNeighborhoodList(Teuchos::RCP<const Epetra_BlockMap> currentOwnedMap,     /* input  */
                                                            Teuchos::RCP<const Epetra_BlockMap> currentOverlapMap,   /* input  */
                                                            int currentNeighborListSize,                             /* input  */
                                                            const int* currentNeighborList,                          /* input  */
                                                            Teuchos::RCP<const Epetra_BlockMap> targetOwnedMap,      /* input  */
                                                            Teuchos::RCP<Epetra_BlockMap> targetOverlapMap,          /* output */
                                                            int& targetNeighborListSize,                             /* output */
                                                            int* targetNeighborList)                                 /* output (allocated within function) */
{
  // Communicate the number of neighbors into the target decomposition

  int currentNumOwnedPoints = currentOwnedMap->NumMyElements();
  int targetNumOwnedPoints = targetOwnedMap->NumMyElements();

  Epetra_Vector currentNumberOfNeighbors(*currentOwnedMap);
  Epetra_Vector targetNumberOfNeighbors(*targetOwnedMap);

  double* currentNumberOfNeighborsPtr;
  currentNumberOfNeighbors.ExtractView(&currentNumberOfNeighborsPtr);

  double* targetNumberOfNeighborsPtr;
  targetNumberOfNeighbors.ExtractView(&targetNumberOfNeighborsPtr);

  int neighborListIndex(0), numNeighbors(0), index(0);
  while(neighborListIndex < currentNeighborListSize){
    numNeighbors = currentNeighborList[neighborListIndex++];
    currentNumberOfNeighborsPtr[index++] = numNeighbors;
    neighborListIndex += numNeighbors;
  }

  Epetra_Import oneDimensionalImporter(*currentOwnedMap, *targetOwnedMap);
  targetNumberOfNeighbors.Import(currentNumberOfNeighbors, oneDimensionalImporter, Insert);

  // Create Epetra_Vectors with variable-length elements to store the neighborhood data in
  // both the current and target decompositions

  // Epetra_BlockMap and Epetra_Map do not allow elements of size zero
  // So, within this function, if a point has no neighbors we'll create
  // a neighbor list of length one and fill it with the value -1

  int numGlobalElements(-1);
  int numMyElements(0);
  int* myGlobalElements(0);
  vector<int> elementSizeList;
  int indexBase(0);

  // Create a vector with variable-length elements and fill it with the global ids of each element's neighbors

  numMyElements = currentNumOwnedPoints;
  myGlobalElements = currentOwnedMap->MyGlobalElements();
  elementSizeList.resize(currentNumOwnedPoints);
  for(int i=0 ; i<currentNumOwnedPoints ; ++i){
    elementSizeList[i] = currentNumberOfNeighbors[i];
    if(elementSizeList[i] == 0)
      elementSizeList[i] = 1;
  }

  Epetra_BlockMap currentNeighborMap(numGlobalElements, numMyElements, myGlobalElements, &elementSizeList[0], indexBase, currentOwnedMap->Comm());
  Epetra_Vector currentNeighbors(currentNeighborMap);
  double* currentNeighborsPtr;
  currentNeighbors.ExtractView(&currentNeighborsPtr);

  neighborListIndex = 0;
  int epetraVectorIndex(0);
  for(int i=0 ; i<currentNumOwnedPoints ; ++i){
    numNeighbors = currentNeighborList[neighborListIndex++];
    if(numNeighbors == 0)
      currentNeighborsPtr[epetraVectorIndex++] = -1;
    for(int j=0 ; j<numNeighbors ; ++j)
      currentNeighborsPtr[epetraVectorIndex++] = currentOverlapMap->GID( currentNeighborList[neighborListIndex++] );
  }

  // Create a vector with variable-length elements in the target configuration to recieve the neighborhood information

  numMyElements = targetNumOwnedPoints;
  myGlobalElements = targetOwnedMap->MyGlobalElements();
  elementSizeList.resize(targetNumOwnedPoints);
  for(int i=0 ; i<targetNumOwnedPoints ; ++i){
    elementSizeList[i] = targetNumberOfNeighborsPtr[i];
    if(elementSizeList[i] == 0)
      elementSizeList[i] = 1;
  }

  Epetra_BlockMap targetNeighborMap(numGlobalElements, numMyElements, myGlobalElements, &elementSizeList[0], indexBase, targetOwnedMap->Comm());
  Epetra_Vector targetNeighbors(targetNeighborMap);
  double* targetNeighborsPtr;
  targetNeighbors.ExtractView(&targetNeighborsPtr);

  // Import the neighborhood data

  Epetra_Import neighborhoodImporter(currentNeighborMap, targetNeighborMap);
  targetNeighbors.Import(currentNeighbors, neighborhoodImporter, Insert);

  // Create a taret overlap map
  int localId, globalId;
  set<int> offProcessorElements;
  for(int i=0 ; i<targetNeighbors.MyLength() ; ++i){
    globalId = targetNeighborsPtr[i];
    localId = targetOwnedMap->LID(globalId);
    if(localId == -1)
      offProcessorElements.insert(globalId);
  }
  int numOffProcessorElements = static_cast<int>( offProcessorElements.size() );

  vector<int> targetOverlapGlobalIds(targetNumOwnedPoints + numOffProcessorElements);
  int* targetOwnedGlobalIds = targetOwnedMap->MyGlobalElements();
  for(int i=0 ; i<targetNumOwnedPoints ; ++i)
    targetOverlapGlobalIds[i] = targetOwnedGlobalIds[i];
  index = 0;
  for(set<int>::iterator it=offProcessorElements.begin() ; it!=offProcessorElements.end() ; it++){
    targetOverlapGlobalIds[targetNumOwnedPoints+index] = *it;
    index += 1;
  }

  // Create the target overlap map

  targetOverlapMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements,
                                                      static_cast<int>( targetOverlapGlobalIds.size() ),
                                                      &targetOverlapGlobalIds[0],
                                                      1,
                                                      0,
                                                      targetOwnedMap->Comm()));

  // Allocate the target neighbor list

  targetNeighborListSize = 0;
  for(int i=0 ; i<targetNumOwnedPoints ; ++i)
    targetNeighborListSize += 1 + targetNumberOfNeighborsPtr[i];
  targetNeighborList = new int[targetNeighborListSize];

  // Fill the target neighbor list

  neighborListIndex = 0;
  for(int i=0 ; i<targetNumOwnedPoints ; ++i){
    int numNeighbors = targetNumberOfNeighborsPtr[i];
    targetNeighborList[neighborListIndex++] = numNeighbors;
    int firstPointInElement = targetNeighborMap.FirstPointInElement(i);
    for(int j=0 ; j<numNeighbors ; ++j){
      globalId = targetNeighborsPtr[firstPointInElement + j];
      localId = targetOverlapMap->LID(globalId);
      targetNeighborList[neighborListIndex++] = localId;
    }
  }
}

void PeridigmNS::ProximitySearch::GlobalProximitySearch(Epetra_Vector& x,                         /* input  */
                                                        double searchRadius,                      /* input  */
                                                        Teuchos::RCP<Epetra_BlockMap> overlapMap, /* output */
                                                        int& neighborListSize,                    /* output */
                                                        int* neighborList)                        /* output (allocated within function) */

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

  dimension = 1;
  Teuchos::RCP<Epetra_BlockMap> currentOwnedMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOwnedMap(originalMap.Comm(), decomp, dimension)));
  Teuchos::RCP<Epetra_BlockMap> currentOverlapMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOverlapMap(originalMap.Comm(), decomp, dimension)));
  Teuchos::RCP<const Epetra_BlockMap> targetOwnedMap = Teuchos::rcpFromRef(originalMap);

  RebalanceNeighborhoodList(currentOwnedMap,                   /* input  */
                            currentOverlapMap,                 /* input  */
                            list.get_size_neighborhood_list(), /* input  */
                            list.get_neighborhood().get(),     /* input  */
                            targetOwnedMap,                    /* input  */
                            overlapMap,                        /* output */
                            neighborListSize,                  /* output */
                            neighborList);                     /* output (allocated within function) */

  // CHECK TO MAKE SURE DECOMP CLEANS UP MEMORY WHEN DELETED
}

