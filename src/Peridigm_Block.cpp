/*! \file Peridigm_Block.cpp */

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

#include "Peridigm_Block.hpp"
#include <vector>
#include <set>

void PeridigmNS::Block::createMapsFromGlobalMaps(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
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

  std::vector<int> IDs;
  IDs.reserve(globalOverlapScalarPointMap->NumMyElements()); // upper bound
  std::vector<int> bondIDs;
  bondIDs.reserve(globalOverlapScalarPointMap->NumMyElements());
  std::vector<int> bondElementSize;
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
  int elementSize = 1;
  int indexBase = 0;
  ownedScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &IDs[0], elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  elementSize = 3;
  ownedVectorPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &IDs[0], elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  numMyElements = bondElementSize.size();
  ownedScalarBondMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &bondIDs[0], &bondElementSize[0], indexBase, globalOwnedScalarPointMap->Comm()));

  // Create a list of nodes that need to be ghosted (both across material boundaries and across processor boundaries)
  std::set<int> ghosts;
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
  std::vector<int> ownedIDs(IDs.begin(), IDs.end());

  // Append ghosts to IDs
  // This creates the overlap global ID list
  for(std::set<int>::iterator it=ghosts.begin() ; it!=ghosts.end() ; ++it)
    IDs.push_back(*it);

  // Create the overlap scalar point map and the overlap vector point map

  numMyElements = IDs.size();
  elementSize = 1;
  overlapScalarPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &IDs[0], elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  elementSize = 3;
  overlapVectorPointMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, &IDs[0], elementSize, indexBase, globalOwnedScalarPointMap->Comm()));

  // Create the NeighborhoodData for this block
  // All the IDs in the neighborhood list are local IDs into the new block-specific overlap map

  std::vector<int> neighborhoodList;
  std::vector<int> neighborhoodPtr(ownedIDs.size());
  int* const globalNeighborhoodPtr = globalNeighborhoodData->NeighborhoodPtr();
  for(unsigned int i=0 ; i<ownedIDs.size() ; ++i){
    neighborhoodPtr[i] = (int)(neighborhoodList.size());
    int globalNeighborhoodListIndex = globalNeighborhoodPtr[globalOverlapScalarPointMap->LID(ownedIDs[i])];
    int numNeighbors = globalNeighborhoodList[globalNeighborhoodListIndex++];
    neighborhoodList.push_back(numNeighbors);
    for(int j=0 ; j<numNeighbors ; ++j){
      int globalNeighborID = globalOverlapScalarPointMap->GID(globalNeighborhoodList[globalNeighborhoodListIndex++]);
      neighborhoodList.push_back( overlapScalarPointMap->LID(globalNeighborID) );
    }
  }

  // convert ownedIDs to local IDs into the new overlap map
  for(unsigned int i=0 ; i<ownedIDs.size() ; ++i)
    ownedIDs[i] = overlapScalarPointMap->LID(ownedIDs[i]);

  // create the NeighborhoodData for this block

  neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
  neighborhoodData->SetNumOwned(ownedIDs.size());
  memcpy(neighborhoodData->OwnedIDs(), 
 		 &ownedIDs[0],
 		 ownedIDs.size()*sizeof(int));
  memcpy(neighborhoodData->NeighborhoodPtr(), 
 		 &neighborhoodPtr[0],
 		 neighborhoodPtr.size()*sizeof(int));
  neighborhoodData->SetNeighborhoodListSize(neighborhoodList.size());
  memcpy(neighborhoodData->NeighborhoodList(),
 		 &neighborhoodList[0],
 		 neighborhoodList.size()*sizeof(int));
}

void PeridigmNS::Block::initializeDataManager(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs)
{
  // The material model and maps must all be set prior to initializing the data manager (and the contact model as well, if there is one).
  // Note that not all the maps are strictly required, so these conditions could be relaxed somewhat.
  TEST_FOR_EXCEPT_MSG(materialModel.is_null(), "\n**** Material model must be set via Block::setMaterialModel() prior to calling Block::initializeDataManager()\n");
  TEST_FOR_EXCEPT_MSG(ownedScalarPointMap.is_null() ||
                      ownedVectorPointMap.is_null() ||
                      overlapScalarPointMap.is_null() ||
                      overlapVectorPointMap.is_null() ||
                      ownedScalarBondMap.is_null(),
                      "\n**** Maps must be set prior to calling Block::initializeDataManager()\n");
  
  // Collect all the required field specs
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  // Specs pass in via fieldSpecs argument (if any)
  if(!fieldSpecs.is_null()){
    specs->insert(specs->end(), fieldSpecs->begin(), fieldSpecs->end());
  }
  // Material model specs
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > materialModelSpecs = materialModel->VariableSpecs();
  specs->insert(specs->end(), materialModelSpecs->begin(), materialModelSpecs->end());
  // Contact model specs (if any)
  if(!contactModel.is_null()){
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > contactModelSpecs = contactModel->VariableSpecs();
    specs->insert(specs->end(), contactModelSpecs->begin(), contactModelSpecs->end());
  }
  
  dataManager = Teuchos::rcp(new PeridigmNS::DataManager);

  dataManager->setMaps(ownedScalarPointMap,
                       overlapScalarPointMap,
                       ownedVectorPointMap,
                       overlapVectorPointMap,
                       ownedScalarBondMap);
  
  // Allocate data in the data manager
  dataManager->allocateData(specs);
}

void PeridigmNS::Block::initializeMaterialModel()
{
  TEST_FOR_EXCEPT_MSG(materialModel.is_null(),
                      "\n**** Material model must be set via Block::setMaterialModel() prior to calling Block::initializeMaterialModel()\n");
  TEST_FOR_EXCEPT_MSG(neighborhoodData.is_null(),
                      "\n**** Neighborhood data must be set via Block::setNeighborhoodData() prior to calling Block::initializeMaterialModel()\n");
  TEST_FOR_EXCEPT_MSG(dataManager.is_null(),
                      "\n**** DataManager must be initialized via Block::initializeDataManager() prior to calling Block::initializeMaterialModel()\n");

  double dt = 0.0;
  materialModel->initialize(dt,
                            neighborhoodData->NumOwnedPoints(),
                            neighborhoodData->OwnedIDs(),
                            neighborhoodData->NeighborhoodList(),
                            *dataManager);
}

void PeridigmNS::Block::importData(const Epetra_Vector& source, Field_NS::FieldSpec spec, Field_ENUM::Step step, Epetra_CombineMode combineMode)
{
  if(dataManager->hasData(spec, step)){

    // scalar data
    if(source.Map().ElementSize() == 1){
      if( oneDimensionalImporter.is_null() || !oneDimensionalImporter->SourceMap().SameAs(source.Map()) )
        oneDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapScalarPointMap(), source.Map()));
      dataManager->getData(spec, step)->Import(source, *oneDimensionalImporter, combineMode);
    }

    // vector data
    else if(source.Map().ElementSize() == 3){
      if( threeDimensionalImporter.is_null() || !threeDimensionalImporter->SourceMap().SameAs(source.Map()) )
        threeDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapVectorPointMap(), source.Map()));
      dataManager->getData(spec, step)->Import(source, *threeDimensionalImporter, combineMode);
    }
  }
}

void PeridigmNS::Block::exportData(Epetra_Vector& target, Field_NS::FieldSpec spec, Field_ENUM::Step step, Epetra_CombineMode combineMode)
{
  if(dataManager->hasData(spec, step)){

    // scalar data
    if(target.Map().ElementSize() == 1){
      if( oneDimensionalImporter.is_null() || !oneDimensionalImporter->SourceMap().SameAs(target.Map()) )
        oneDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapScalarPointMap(), target.Map()));
      target.Export(*(dataManager->getData(spec, step)), *oneDimensionalImporter, combineMode);  
    }

    // vector data
    else if(target.Map().ElementSize() == 3){
      if( threeDimensionalImporter.is_null() || !threeDimensionalImporter->SourceMap().SameAs(target.Map()) )
        threeDimensionalImporter = Teuchos::rcp(new Epetra_Import(*dataManager->getOverlapVectorPointMap(), target.Map()));
      target.Export(*(dataManager->getData(spec, step)), *threeDimensionalImporter, combineMode);  
    }
  }
}
