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

#include <set>
#include <Epetra_Import.h>
#include "Peridigm_ProximitySearch.hpp"
#include "QuickGrid.h"
#include "PdZoltan.h"
#include "NeighborhoodList.h"

using std::vector;
using std::map;
using std::pair;
using std::cout;
using std::set;

void PeridigmNS::ProximitySearch::RebalanceNeighborhoodList(Teuchos::RCP<const Epetra_BlockMap> currentOwnedMap,     /* input  */
                                                            Teuchos::RCP<const Epetra_BlockMap> currentOverlapMap,   /* input  */
                                                            int currentNeighborListSize,                             /* input  */
                                                            const int* currentNeighborList,                          /* input  */
                                                            Teuchos::RCP<const Epetra_BlockMap> targetOwnedMap,      /* input  */
                                                            Teuchos::RCP<Epetra_BlockMap>& targetOverlapMap,         /* output */
                                                            int& targetNeighborListSize,                             /* output */
                                                            int*& targetNeighborList)                                /* output (allocated within function) */
{
  // Communicate the number of neighbors into the target decomposition

  int currentNumOwnedPoints = currentOwnedMap->NumMyElements();
  int targetNumOwnedPoints = targetOwnedMap->NumMyElements();

  Epetra_Vector currentNumberOfNeighbors(*currentOwnedMap);
  Epetra_Vector targetNumberOfNeighbors(*targetOwnedMap);

  double* currentNumberOfNeighborsPtr;
  currentNumberOfNeighbors.ExtractView(&currentNumberOfNeighborsPtr);

  int neighborListIndex(0), numNeighbors(0), index(0);
  while(neighborListIndex < currentNeighborListSize){
    numNeighbors = currentNeighborList[neighborListIndex++];
    currentNumberOfNeighborsPtr[index++] = numNeighbors;
    neighborListIndex += numNeighbors;
  }

  Epetra_Import oneDimensionalImporter(*targetOwnedMap, *currentOwnedMap);
  targetNumberOfNeighbors.Import(currentNumberOfNeighbors, oneDimensionalImporter, Insert);

  double* targetNumberOfNeighborsPtr;
  targetNumberOfNeighbors.ExtractView(&targetNumberOfNeighborsPtr);

  // Create Epetra_Vectors with variable-length elements to store the neighborhood data in
  // both the current and target decompositions

  // Epetra_BlockMap and Epetra_Map do not allow elements of size zero.
  // So, within this function, if a point has no neighbors we'll create
  // a neighbor list of length one and fill it with the value -1.

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
    numNeighbors = static_cast<int>( currentNumberOfNeighbors[i] );
    if(numNeighbors == 0)
      elementSizeList[i] = 1;
    else
      elementSizeList[i] = numNeighbors;
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
    numNeighbors = static_cast<int>( targetNumberOfNeighborsPtr[i] );
    if(numNeighbors == 0)
      elementSizeList[i] = 1;
    else
      elementSizeList[i] = numNeighbors;
  }

  Epetra_BlockMap targetNeighborMap(numGlobalElements, numMyElements, myGlobalElements, &elementSizeList[0], indexBase, targetOwnedMap->Comm());
  Epetra_Vector targetNeighbors(targetNeighborMap);
  double* targetNeighborsPtr;
  targetNeighbors.ExtractView(&targetNeighborsPtr);

  // Import the neighborhood data

  Epetra_Import neighborhoodImporter(targetNeighborMap, currentNeighborMap);
  targetNeighbors.Import(currentNeighbors, neighborhoodImporter, Insert);

  // Create a taret overlap map
  int localId, globalId;
  set<int> offProcessorElements;
  for(int i=0 ; i<targetNeighbors.MyLength() ; ++i){
    globalId = static_cast<int>( targetNeighborsPtr[i] );
    // Note that in this list a globalId of -1 denotes zero neighbors
    if(globalId != -1){
      localId = targetOwnedMap->LID(globalId);
      if(localId == -1)
        offProcessorElements.insert(globalId);
    }
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
    targetNeighborListSize += 1 + static_cast<int>( targetNumberOfNeighborsPtr[i] );
  targetNeighborList = new int[targetNeighborListSize];

  // Fill the target neighbor list
  neighborListIndex = 0;
  int firstPointInElement;
  for(int i=0 ; i<targetNumOwnedPoints ; ++i){
    numNeighbors = static_cast<int>( targetNumberOfNeighborsPtr[i] );
    targetNeighborList[neighborListIndex++] = numNeighbors;
    firstPointInElement = targetNeighborMap.FirstPointInElement(i);
    for(int j=0 ; j<numNeighbors ; ++j){
      globalId = static_cast<int>( targetNeighborsPtr[firstPointInElement + j] );
      localId = targetOverlapMap->LID(globalId);
      targetNeighborList[neighborListIndex++] = localId;
    }
  }
}

void PeridigmNS::ProximitySearch::GlobalProximitySearch(Teuchos::RCP<Epetra_Vector> x,                                              /* input  */
                                                        Teuchos::RCP<Epetra_Vector> searchRadii,                                    /* input  */
                                                        Teuchos::RCP<Epetra_BlockMap>& overlapMap,                                  /* output */
                                                        int& neighborListSize,                                                      /* output */
                                                        int*& neighborList,                                                         /* output (allocated within function) */
                                                        std::vector< std::shared_ptr<PdBondFilter::BondFilter> > bondFilters,       /* optional input */
                                                        double radiusAddition)                                                      /* optional input */

{
  // The proximity search does not appear to function properly if any of the search radii are set to zero
  for(int i=0 ; i<searchRadii->MyLength() ; ++i){
    TEUCHOS_TEST_FOR_TERMINATION((*searchRadii)[i] <= 0.0, "\n****Error:  PeridigmNS::ProximitySearch::GlobalProximitySearch(), search radii must be greater than or equal to zero.\n");
  }

  // Copy information from the Epetra_Vector into a QUICKGRID::Data object
  const Epetra_BlockMap& originalMap = x->Map();
  int dimension = 3;
  int numMyElements = originalMap.NumMyElements();
  int numGlobalElements = originalMap.NumGlobalElements();
  int* globalIds = originalMap.MyGlobalElements();
  double* coordinates;
  x->ExtractView(&coordinates);
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(numMyElements, dimension);
  decomp.globalNumPoints = numGlobalElements;
  int* decompGlobalIds = decomp.myGlobalIDs.get();
  double* decompVolumes = decomp.cellVolume.get();
  double* decompCoordinates = decomp.myX.get();
  for(int i=0 ; i<numMyElements ; ++i){
    decompGlobalIds[i] = globalIds[i];
    decompVolumes[i] = 1.0;  // use dummy value
    decompCoordinates[i*3]   = coordinates[i*3];
    decompCoordinates[i*3+1] = coordinates[i*3+1];
    decompCoordinates[i*3+2] = coordinates[i*3+2];
  }

  // Rebalance the decomp (RCB decomposition via Zoltan)
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

  // Obtain a horizon list in the rebalanced decomposition
  Epetra_BlockMap rebalancedMap(decomp.globalNumPoints, decomp.numPoints, decomp.myGlobalIDs.get(), 1, 0, originalMap.Comm());
  Epetra_Import searchRadiiImporter(rebalancedMap, searchRadii->Map());
  Teuchos::RCP<Epetra_Vector> rebalancedSearchRadii = Teuchos::rcp(new Epetra_Vector(rebalancedMap));
  rebalancedSearchRadii->Import(*searchRadii, searchRadiiImporter, Insert);

  // Add in radiusAddition, if any.
  // Note that only rebalancedSearchRadii is modified.
  // The input vector, searchRadii, is unchanged
  if(radiusAddition != 0.0){
    for(int i=0 ; i<rebalancedSearchRadii->MyLength() ; ++i)
      (*rebalancedSearchRadii)[i] += radiusAddition;
  }

  // The initial implementation supports only a single bond filter
  if(bondFilters.size() == 0){
    std::shared_ptr<PdBondFilter::BondFilter> bondFilterPtr(new PdBondFilter::BondFilterDefault(false));
    bondFilters.push_back(bondFilterPtr);
  }
 
  // Execute neighbor search
  std::shared_ptr<const Epetra_Comm> commSp(&originalMap.Comm(), NonDeleter<const Epetra_Comm>());

  PDNEIGH::NeighborhoodList list(commSp,
                                 decomp.zoltanPtr.get(),
                                 decomp.numPoints,
                                 decomp.myGlobalIDs,
                                 decomp.myX,
                                 rebalancedSearchRadii,
                                 bondFilters);

  // The neighbor search is complete, but needs to be brought back into the initial decomposition

  int listNumOwned = list.get_num_owned_points();
  int listNumShared = list.get_num_shared_points();
  int listNumTotal = listNumOwned + listNumShared;
  vector<int> overlapIds(listNumTotal);
  int* listOwnedIds = list.get_owned_gids().get();
  int* listSharedIds = list.get_shared_gids().get();
  for(int i=0 ; i<listNumOwned ; ++i)
    overlapIds[i] = listOwnedIds[i];
  for(int i=0 ; i<listNumShared ; ++i)
    overlapIds[i+listNumOwned] = listSharedIds[i];

  Teuchos::RCP<Epetra_BlockMap> currentOwnedMap = Teuchos::rcp(new Epetra_BlockMap(-1, listNumOwned, &overlapIds[0], 1, 0, originalMap.Comm()));
  Teuchos::RCP<Epetra_BlockMap> currentOverlapMap = Teuchos::rcp(new Epetra_BlockMap(-1, listNumTotal, &overlapIds[0], 1, 0, originalMap.Comm()));
  Teuchos::RCP<const Epetra_BlockMap> targetOwnedMap = Teuchos::rcp(new Epetra_BlockMap(-1, originalMap.NumMyElements(), originalMap.MyGlobalElements(), 1, 0, originalMap.Comm()));

  RebalanceNeighborhoodList(currentOwnedMap,                     /* input  */
                            currentOverlapMap,                   /* input  */
                            list.get_size_neighborhood_list(),   /* input  */
                            list.get_local_neighborhood().get(), /* input  */
                            targetOwnedMap,                      /* input  */
                            overlapMap,                          /* output */
                            neighborListSize,                    /* output */
                            neighborList);                       /* output (allocated within function) */
}

