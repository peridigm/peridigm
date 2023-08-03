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
#include "Peridigm_Field.hpp"
#include <vector>
#include <set>
#include <sstream>

using std::vector;

void PeridigmNS::Block::initialize(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOwnedVectorPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOverlapVectorPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarBondMap,
                                   Teuchos::RCP<const Epetra_Vector>   globalBlockIds,
                                   Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData)
{
  // Initialize the base class to create maps and neighborhood data
  BlockBase::initialize(globalOwnedScalarPointMap,
                        globalOverlapScalarPointMap,
                        globalOwnedVectorPointMap,
                        globalOverlapVectorPointMap,
                        globalOwnedScalarBondMap,
                        globalBlockIds,
                        globalNeighborhoodData);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(materialModel.is_null(),
                              "\n**** Material model must be set via Block::setMaterialModel() prior to calling Block::initialize()\n");

  // Collect all the required field Ids
  vector<int> fieldIds;
  // Ids passed in via setAuxiliaryFieldIds(), if any
  fieldIds.insert(fieldIds.end(), auxiliaryFieldIds.begin(), auxiliaryFieldIds.end());
  // Material model field Ids
  vector<int> materialModelFieldIds = materialModel->FieldIds();
  fieldIds.insert(fieldIds.end(), materialModelFieldIds.begin(), materialModelFieldIds.end());
  // Damage model field Ids (if any)
  if(!damageModel.is_null()){
    vector<int> damageModelFieldIds = damageModel->FieldIds();
    fieldIds.insert(fieldIds.end(), damageModelFieldIds.begin(), damageModelFieldIds.end());
  }

  BlockBase::initializeDataManager(fieldIds);
}

void PeridigmNS::Block::initializeMaterialModel(double timeStep)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(materialModel.is_null(),
                      "\n**** Material model must be set via Block::setMaterialModel() prior to calling Block::initializeMaterialModel()\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(neighborhoodData.is_null(),
                      "\n**** Neighborhood data must be set via Block::setNeighborhoodData() prior to calling Block::initializeMaterialModel()\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(dataManager.is_null(),
                      "\n**** DataManager must be initialized via Block::initializeDataManager() prior to calling Block::initializeMaterialModel()\n");

  materialModel->initialize(timeStep,
                            neighborhoodData->NumOwnedPoints(),
                            neighborhoodData->OwnedIDs(),
                            neighborhoodData->NeighborhoodList(),
                            *dataManager);
}

void PeridigmNS::Block::initializeDamageModel(double timeStep)
{
  if(damageModel.is_null())
    return;

  TEUCHOS_TEST_FOR_EXCEPT_MSG(neighborhoodData.is_null(),
                      "\n**** Neighborhood data must be set via Block::setNeighborhoodData() prior to calling Block::initializeDamageModel()\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(dataManager.is_null(),
                      "\n**** DataManager must be initialized via Block::initializeDataManager() prior to calling Block::initializeDamageModel()\n");

  damageModel->initialize(timeStep,
                          neighborhoodData->NumOwnedPoints(),
                          neighborhoodData->OwnedIDs(),
                          neighborhoodData->NeighborhoodList(),
                          *dataManager);
}

PeridigmNS::DataManagerSynchronizer& PeridigmNS::DataManagerSynchronizer::self() {
  static DataManagerSynchronizer dataManagerSynchronizer;
  return dataManagerSynchronizer;
}

void PeridigmNS::DataManagerSynchronizer::initialize(Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap,
                                                     Teuchos::RCP<const Epetra_BlockMap> threeDimensionalMap) {
  scalarScratch = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  scalarSum     = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  vectorScratch = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  vectorSum     = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
}

void PeridigmNS::DataManagerSynchronizer::setFieldIdsToSynchronizeAfterInitialize(std::vector<int>& fieldIds)
{
  fieldIdsToSychAfterInitialize.insert(fieldIdsToSychAfterInitialize.end(), fieldIds.begin(), fieldIds.end());
  // Remove duplicates
  std::sort( fieldIdsToSychAfterInitialize.begin(), fieldIdsToSychAfterInitialize.end() );
  fieldIdsToSychAfterInitialize.erase( unique( fieldIdsToSychAfterInitialize.begin(),
                                               fieldIdsToSychAfterInitialize.end() ),
                                       fieldIdsToSychAfterInitialize.end() );
}

void PeridigmNS::DataManagerSynchronizer::setFieldIdsToSynchronizeAfterPrecompute(std::vector<int>& fieldIds)
{
  fieldIdsToSychAfterPrecompute.insert(fieldIdsToSychAfterPrecompute.end(), fieldIds.begin(), fieldIds.end());
  // Remove duplicates
  std::sort( fieldIdsToSychAfterPrecompute.begin(), fieldIdsToSychAfterPrecompute.end() );
  fieldIdsToSychAfterPrecompute.erase( unique( fieldIdsToSychAfterPrecompute.begin(),
                                               fieldIdsToSychAfterPrecompute.end() ),
                                       fieldIdsToSychAfterPrecompute.end() );
}

void PeridigmNS::DataManagerSynchronizer::checkFieldValidity(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks)
{
  std::vector<int> fIds = fieldIdsToSychAfterInitialize;
  fIds.insert(fIds.end(), fieldIdsToSychAfterPrecompute.begin(), fieldIdsToSychAfterPrecompute.end());

  for (auto fieldId : fIds) {

    FieldSpec fieldSpec = PeridigmNS::FieldManager::self().getFieldSpec(fieldId);

    PeridigmField::Relation relation = fieldSpec.getRelation();
    TEUCHOS_TEST_FOR_TERMINATION((relation != PeridigmField::NODE && relation != PeridigmField::ELEMENT),
                                "**** Error in DataManagerSynchronizer::checkFieldValidity():  Parallel synchronization available only for node and elemet variables.\n");

    PeridigmField::Length length = fieldSpec.getLength();
    TEUCHOS_TEST_FOR_TERMINATION((length != PeridigmField::SCALAR && length != PeridigmField::VECTOR),
                                "**** Error in DataManagerSynchronizer::checkFieldValidity():  Parallel synchronization available only for scalar and vector variables.\n");

    PeridigmField::Temporal temporal = fieldSpec.getTemporal();
    PeridigmField::Step step = PeridigmField::STEP_NONE;
    if (temporal == PeridigmField::TWO_STEP) {
      step = PeridigmField::STEP_NP1;
    }

    for(std::vector<PeridigmNS::Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      bool hasField = blockIt->getDataManager()->hasData(fieldId, step);
      if (!hasField) {
        std::stringstream ss;
        ss << "**** Error in DataManagerSynchronizer::checkFieldValidity():  Field ";
        ss << fieldSpec.getLabel() << " not present in block " << blockIt->getName() << ".\n";
        TEUCHOS_TEST_FOR_TERMINATION(!hasField, ss.str());
      }
    }
  }
}

void PeridigmNS::DataManagerSynchronizer::synchronizeDataAfterInitialize(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks)
{
  synchronize(blocks, fieldIdsToSychAfterInitialize);
}

void PeridigmNS::DataManagerSynchronizer::synchronizeDataAfterPrecompute(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks)
{
  synchronize(blocks, fieldIdsToSychAfterPrecompute);
}

void PeridigmNS::DataManagerSynchronizer::synchronize(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks,
                                                      std::vector<int> fieldIds)
{
  if (fieldIds.size() == 0)
    return;

  std::vector<PeridigmNS::Block>::iterator blockIt;

  for (auto fieldId : fieldIds) {

    FieldSpec fieldSpec = PeridigmNS::FieldManager::self().getFieldSpec(fieldId);

    PeridigmField::Temporal temporal = fieldSpec.getTemporal();
    PeridigmField::Step step = PeridigmField::STEP_NONE;
    if (temporal == PeridigmField::TWO_STEP) {
      step = PeridigmField::STEP_NP1;
    }

    PeridigmField::Length length = fieldSpec.getLength();
    if (length == PeridigmField::SCALAR) {
      scalarSum->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scalarScratch->PutScalar(0.0);
        blockIt->exportData(scalarScratch, fieldId, step, Add);
        scalarSum->Update(1.0, *(scalarScratch), 1.0);
      }
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(scalarSum, fieldId, step, Insert);
      }
    }
    else if(length == PeridigmField::VECTOR) {
      vectorSum->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        vectorScratch->PutScalar(0.0);
        blockIt->exportData(vectorScratch, fieldId, step, Add);
        vectorSum->Update(1.0, *(vectorScratch), 1.0);
      }
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(vectorSum, fieldId, step, Insert);
      }
    }
    else {
      TEUCHOS_TEST_FOR_TERMINATION((length != PeridigmField::SCALAR && length != PeridigmField::VECTOR),
                                  "**** Error in DataManagerSynchronizer::synchronize():  Parallel synchronization available only for scalar and vector variables.\n");
    }
  }
}
