/*! \file Peridigm_BlockBase.cpp */

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

#include "Peridigm_BlockBase.hpp"
#include "Peridigm_Field.hpp"
#include <vector>
#include <set>

using namespace std;

PeridigmNS::BlockBase::BlockBase(std::string blockName_, int blockID_, Teuchos::ParameterList& blockParams_)
  : blockName(blockName_), blockID(blockID_), blockParams(blockParams_)
{}

void PeridigmNS::BlockBase::initialize(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOwnedVectorPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOverlapVectorPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarBondMap,
                                   Teuchos::RCP<const Epetra_Vector>   globalBlockIds,
                                   Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData)
{
  createMapsFromGlobalMaps(globalOwnedScalarPointMap,
                           globalOverlapScalarPointMap,
                           globalOwnedVectorPointMap,
                           globalOverlapVectorPointMap,
                           globalOwnedScalarBondMap,
                           globalBlockIds,
                           globalNeighborhoodData);

  neighborhoodData = createNeighborhoodDataFromGlobalNeighborhoodData(globalOverlapScalarPointMap,
                                                                      globalNeighborhoodData);
}

void PeridigmNS::BlockBase::importData(const Epetra_Vector& source, int fieldId, PeridigmField::Step step, Epetra_CombineMode combineMode)
{
  if(dataManager->hasData(fieldId, step)){

    // scalar data
    if(source.Map().ElementSize() == 1){
      if(oneDimensionalImporter.is_null())
        oneDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapScalarPointMap(), source.Map()));
      dataManager->getData(fieldId, step)->Import(source, *oneDimensionalImporter, combineMode);
    }

    // vector data
    else if(source.Map().ElementSize() == 3){
      if(threeDimensionalImporter.is_null())
        threeDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapVectorPointMap(), source.Map()));
      dataManager->getData(fieldId, step)->Import(source, *threeDimensionalImporter, combineMode);
    }
  }
}

void PeridigmNS::BlockBase::exportData(Epetra_Vector& target, int fieldId, PeridigmField::Step step, Epetra_CombineMode combineMode)
{
  if(dataManager->hasData(fieldId, step)){

    // scalar data
    if(target.Map().ElementSize() == 1){
      if(oneDimensionalImporter.is_null())
        oneDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapScalarPointMap(), target.Map()));
      target.Export(*(dataManager->getData(fieldId, step)), *oneDimensionalImporter, combineMode);  
    }

    // vector data
    else if(target.Map().ElementSize() == 3){
      if(threeDimensionalImporter.is_null())
        threeDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapVectorPointMap(), target.Map()));
      target.Export(*(dataManager->getData(fieldId, step)), *threeDimensionalImporter, combineMode);  
    }
  }
}

void PeridigmNS::BlockBase::createMapsFromGlobalMaps(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
                                                     Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                                     Teuchos::RCP<const Epetra_BlockMap> globalOwnedVectorPointMap,
                                                     Teuchos::RCP<const Epetra_BlockMap> globalOverlapVectorPointMap,
                                                     Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarBondMap,
                                                     Teuchos::RCP<const Epetra_Vector>   globalBlockIds,
                                                     Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData)
{
  double* globalBlockIdsPtr;
  globalBlockIds->ExtractView(&globalBlockIdsPtr);

  // Create a list of all the on-processor elements that are part of this block

  vector<int> IDs;
  IDs.reserve(globalOverlapScalarPointMap->NumMyElements()); // upper bound
  vector<int> bondIDs;
  bondIDs.reserve(globalOverlapScalarPointMap->NumMyElements());
  vector<int> bondElementSize;
  bondElementSize.reserve(globalOwnedScalarPointMap->NumMyElements());

  for(int iLID=0 ; iLID<globalOwnedScalarPointMap->NumMyElements() ; ++iLID){
    if(globalBlockIdsPtr[iLID] == blockID) {
      int globalID = globalOwnedScalarPointMap->GID(iLID);
      IDs.push_back(globalID);
    }
  }

  // Record the size of these elements in the bond map
  // Note that if an element has no bonds, it has no entry in the bondMap
  // So, the bond map and the scalar map can have a different number of entries (different local IDs)

  for(int iLID=0 ; iLID<globalOwnedScalarBondMap->NumMyElements() ; ++iLID){
    int globalID = globalOwnedScalarBondMap->GID(iLID);
    int localID = globalOwnedScalarPointMap->LID(globalID);
    if(globalBlockIdsPtr[localID] == blockID){
      bondIDs.push_back(globalID);
      bondElementSize.push_back(globalOwnedScalarBondMap->ElementSize(iLID));
    }
  }

  // Create the owned scalar point map, the owned vector point map, and the owned scalar bond map

  int numGlobalElements = -1;
  int numMyElements = IDs.size();
  int* myGlobalElements = 0;
  if(numMyElements > 0)
    myGlobalElements = &IDs.at(0);
  int elementSize = 1;
  int indexBase = 0;
  ownedScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  elementSize = 3;
  ownedVectorPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  numMyElements = bondElementSize.size();
  myGlobalElements = 0;
  int* elementSizeList = 0;
  if(numMyElements > 0){
    myGlobalElements = &bondIDs.at(0);
    elementSizeList = &bondElementSize.at(0);
  }
  ownedScalarBondMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, globalOwnedScalarPointMap->Comm()));

  // Create a list of nodes that need to be ghosted (both across material boundaries and across processor boundaries)
  set<int> ghosts;

  // Check the neighborhood list for things that need to be ghosted
  int* const globalNeighborhoodList = globalNeighborhoodData->NeighborhoodList();
  int globalNeighborhoodListIndex = 0;
  for(int iLID=0 ; iLID<globalNeighborhoodData->NumOwnedPoints() ; ++iLID){
    int numNeighbors = globalNeighborhoodList[globalNeighborhoodListIndex++];
    if(globalBlockIdsPtr[iLID] == blockID) {
      for(int i=0 ; i<numNeighbors ; ++i){
        int neighborGlobalID = globalOverlapScalarPointMap->GID( globalNeighborhoodList[globalNeighborhoodListIndex + i] );
        ghosts.insert(neighborGlobalID);
      }
    }
    globalNeighborhoodListIndex += numNeighbors;
  }

  // Remove entries from ghosts that are already in IDs
  for(unsigned int i=0 ; i<IDs.size() ; ++i)
    ghosts.erase(IDs[i]);

  // Copy IDs, this is the owned global ID list
  vector<int> ownedIDs(IDs.begin(), IDs.end());

  // Append ghosts to IDs
  // This creates the overlap global ID list
  for(set<int>::iterator it=ghosts.begin() ; it!=ghosts.end() ; ++it)
    IDs.push_back(*it);

  // Create the overlap scalar point map and the overlap vector point map

  numMyElements = IDs.size();
  myGlobalElements = 0;
  if(numMyElements > 0)
    myGlobalElements = &IDs.at(0);
  elementSize = 1;
  overlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  elementSize = 3;
  overlapVectorPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  // Invalidate the importers
  oneDimensionalImporter = Teuchos::RCP<Epetra_Import>();
  threeDimensionalImporter = Teuchos::RCP<Epetra_Import>();
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> PeridigmNS::BlockBase::createNeighborhoodDataFromGlobalNeighborhoodData(Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                                                                                               Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData)
{
  int numOwnedPoints = ownedScalarPointMap->NumMyElements();
  int* ownedPointGlobalIDs = ownedScalarPointMap->MyGlobalElements();

  vector<int> ownedIDs(numOwnedPoints);
  vector<int> neighborhoodList;
  vector<int> neighborhoodPtr(numOwnedPoints);

  int* const globalNeighborhoodList = globalNeighborhoodData->NeighborhoodList();
  int* const globalNeighborhoodPtr = globalNeighborhoodData->NeighborhoodPtr();

  // Create the neighborhoodList and neighborhoodPtr for this block.
  // All the IDs in the neighborhoodList and neighborhoodPtr are local IDs into 
  // the block-specific overlap map.

  for(int i=0 ; i<numOwnedPoints ; ++i){
    neighborhoodPtr[i] = (int)(neighborhoodList.size());
    int globalID = ownedPointGlobalIDs[i];
    ownedIDs[i] = overlapScalarPointMap->LID(globalID);
    int globalNeighborhoodListIndex = globalNeighborhoodPtr[globalOverlapScalarPointMap->LID(globalID)];
    int numNeighbors = globalNeighborhoodList[globalNeighborhoodListIndex++];
    neighborhoodList.push_back(numNeighbors);
    for(int j=0 ; j<numNeighbors ; ++j){
      int globalNeighborID = globalOverlapScalarPointMap->GID(globalNeighborhoodList[globalNeighborhoodListIndex++]);
      neighborhoodList.push_back( overlapScalarPointMap->LID(globalNeighborID) );
    }
  }

  // create the NeighborhoodData for this block

  Teuchos::RCP<PeridigmNS::NeighborhoodData> blockNeighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
  blockNeighborhoodData->SetNumOwned(ownedIDs.size());
  if(ownedIDs.size() > 0){
    memcpy(blockNeighborhoodData->OwnedIDs(), 
           &ownedIDs.at(0),
           ownedIDs.size()*sizeof(int));
  }
  if(neighborhoodPtr.size() > 0){
    memcpy(blockNeighborhoodData->NeighborhoodPtr(), 
           &neighborhoodPtr.at(0),
           neighborhoodPtr.size()*sizeof(int));
  }
  blockNeighborhoodData->SetNeighborhoodListSize(neighborhoodList.size());
  if(neighborhoodList.size() > 0){
    memcpy(blockNeighborhoodData->NeighborhoodList(),
           &neighborhoodList.at(0),
           neighborhoodList.size()*sizeof(int));
  }

  return blockNeighborhoodData;
}

void PeridigmNS::BlockBase::initializeDataManager(vector<int> fieldIds)
{
  // The material model must be set prior to initializing the data manager.
  // Note that not all the maps are strictly required, so these conditions could be relaxed somewhat.
  TEUCHOS_TEST_FOR_EXCEPT_MSG(ownedScalarPointMap.is_null() ||
                              ownedVectorPointMap.is_null() ||
                              overlapScalarPointMap.is_null() ||
                              overlapVectorPointMap.is_null() ||
                              ownedScalarBondMap.is_null(),
                              "\n**** Maps must be set prior to calling BlockBase::initializeDataManager()\n");
  
  dataManager = Teuchos::rcp(new PeridigmNS::DataManager);

  dataManager->setMaps(ownedScalarPointMap,
                       overlapScalarPointMap,
                       ownedVectorPointMap,
                       overlapVectorPointMap,
                       ownedScalarBondMap);
  
  // remove duplicates
  sort(fieldIds.begin(), fieldIds.end());
  vector<int>::iterator newEnd = unique(fieldIds.begin(), fieldIds.end());
  fieldIds.erase(newEnd, fieldIds.end());

  // Allocate data in the data manager
  dataManager->allocateData(fieldIds);
}
