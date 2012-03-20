/*! \file Peridigm_DataManager.cpp */

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

#include <Teuchos_Exceptions.hpp>
#include <Epetra_Import.h>
#include "Peridigm_DataManager.hpp"

void PeridigmNS::DataManager::allocateData(Teuchos::RCP< std::vector<Field_NS::FieldSpec> > specs)
{
  fieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>(*specs));
 
  // remove duplicates
  sort(fieldSpecs->begin(), fieldSpecs->end());
  std::vector<Field_NS::FieldSpec>::iterator newEnd = unique(fieldSpecs->begin(), fieldSpecs->end());
  fieldSpecs->erase(newEnd, fieldSpecs->end());

  statelessScalarPointFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statelessVectorPointFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statelessScalarBondFieldSpecs  = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulScalarPointFieldSpecs  = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulVectorPointFieldSpecs  = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulScalarBondFieldSpecs   = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  for(unsigned int i=0; i<fieldSpecs->size() ; ++i){
    Field_NS::FieldSpec& spec = (*fieldSpecs)[i];
    // scalar point data
    if(spec.getLength() == Field_ENUM::SCALAR && spec.getRelation() == Field_ENUM::ELEMENT){
      if(spec.get_temporal() == Field_ENUM::CONSTANT)
        statelessScalarPointFieldSpecs->push_back(spec);
      else
        statefulScalarPointFieldSpecs->push_back(spec);
    }
    // vector point data
    else if(spec.getLength() == Field_ENUM::VECTOR3D && spec.getRelation() == Field_ENUM::NODE){
      if(spec.get_temporal() == Field_ENUM::CONSTANT)
        statelessVectorPointFieldSpecs->push_back(spec);
      else
        statefulVectorPointFieldSpecs->push_back(spec);
    }
    // scalar bond data
    else if(spec.getLength() == Field_ENUM::SCALAR && spec.getRelation() == Field_ENUM::BOND){
      if(spec.get_temporal() == Field_ENUM::CONSTANT)
        statelessScalarBondFieldSpecs->push_back(spec);
      else
        statefulScalarBondFieldSpecs->push_back(spec);
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::RangeError, 
                         "PeridigmNS::DataManager::allocateData, invalid FieldSpec!");
    }
  }

  // make sure maps exist before trying to create states
  if(statelessScalarPointFieldSpecs->size() + statefulScalarPointFieldSpecs->size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(overlapScalarPointMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate scalar data with no map (forget setMaps()?).");
  if(statelessVectorPointFieldSpecs->size() + statefulVectorPointFieldSpecs->size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(overlapVectorPointMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate vector data with no map (forget setMaps()?).");
  if(statelessScalarBondFieldSpecs->size() + statefulScalarBondFieldSpecs->size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(ownedScalarBondMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate bond data with no map (forget setMaps()?).");

  // create the states
  if(statelessScalarPointFieldSpecs->size() + statelessVectorPointFieldSpecs->size() + statelessScalarBondFieldSpecs->size() > 0){
    stateNONE = Teuchos::rcp(new State);
    if(statelessScalarPointFieldSpecs->size() > 0)
      stateNONE->allocateScalarPointData(statelessScalarPointFieldSpecs, overlapScalarPointMap);
    if(statelessVectorPointFieldSpecs->size() > 0)
      stateNONE->allocateVectorPointData(statelessVectorPointFieldSpecs, overlapVectorPointMap);
    if(statelessScalarBondFieldSpecs->size() > 0)
      stateNONE->allocateScalarBondData(statelessScalarBondFieldSpecs, ownedScalarBondMap);
  }
  if(statefulScalarPointFieldSpecs->size() + statefulVectorPointFieldSpecs->size() + statefulScalarBondFieldSpecs->size() > 0){
    stateN = Teuchos::rcp(new State);
    stateNP1 = Teuchos::rcp(new State);
    if(statefulScalarPointFieldSpecs->size() > 0){
      stateN->allocateScalarPointData(statefulScalarPointFieldSpecs, overlapScalarPointMap);
      stateNP1->allocateScalarPointData(statefulScalarPointFieldSpecs, overlapScalarPointMap);
    }
    if(statefulVectorPointFieldSpecs->size() > 0){
      stateN->allocateVectorPointData(statefulVectorPointFieldSpecs, overlapVectorPointMap);
      stateNP1->allocateVectorPointData(statefulVectorPointFieldSpecs, overlapVectorPointMap);
    }   
    if(statefulScalarBondFieldSpecs->size() > 0){
      stateN->allocateScalarBondData(statefulScalarBondFieldSpecs, ownedScalarBondMap);
      stateNP1->allocateScalarBondData(statefulScalarBondFieldSpecs, ownedScalarBondMap);
    }
  }
}

void PeridigmNS::DataManager::scatterToGhosts()
{
  // goal:
  //   for each global ID, copy values from the processor that owns it to the
  //   processors that have a ghosted copy.
  // approach:
  //   1) create non-overlap multivector
  //   2) copy locally-owned data to it
  //   3) scatter back to the overlap multivector

  for(int iState=0 ; iState<3 ; ++iState){

    Teuchos::RCP<State> state = stateNONE;
    if(iState == 1)
      state = stateN;
    else if(iState == 2)
      state = stateNP1;

    // process scalar data
    Teuchos::RCP<Epetra_MultiVector> overlapScalarPointMultiVector = state->getScalarPointMultiVector();
    if(!overlapScalarPointMultiVector.is_null()){

      TEUCHOS_TEST_FOR_EXCEPTION(overlapScalarPointMap.is_null() || ownedScalarPointMap.is_null(), Teuchos::NullReferenceError,
                         "Error in PeridigmNS::DataManager::scatterToGhosts(), inconsistent scalar maps.");

      // create an owned (non-overlap) multivector
      int numVectors = overlapScalarPointMultiVector->NumVectors();
      Teuchos::RCP<Epetra_MultiVector> ownedScalarPointMultiVector = Teuchos::rcp(new Epetra_MultiVector(*ownedScalarPointMap, numVectors));

      // copy data from the overlap vector into the owned (non-overlap) vector
      for(int iVec=0 ; iVec<numVectors ; ++iVec){
        double* overlapScalarPointData = (*overlapScalarPointMultiVector)[iVec];
        double* ownedScalarPointData = (*ownedScalarPointMultiVector)[iVec];
        for(int iLID=0 ; iLID<ownedScalarPointMap->NumMyElements() ; ++iLID){
          int globalID = ownedScalarPointMap->GID(iLID);
          int overlapMapLocalID = overlapScalarPointMap->LID(globalID);
          ownedScalarPointData[iLID] = overlapScalarPointData[overlapMapLocalID];
        }
      }

      // scatter the data back from the owned (non-overlap) multivector into the overlap multivector
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*overlapScalarPointMap, *ownedScalarPointMap));
      overlapScalarPointMultiVector->Import(*ownedScalarPointMultiVector, *importer, Insert);
    }

    // process vector data
    Teuchos::RCP<Epetra_MultiVector>  overlapVectorPointMultiVector = state->getVectorPointMultiVector();
    if(!overlapVectorPointMultiVector.is_null()){

      TEUCHOS_TEST_FOR_EXCEPTION(overlapVectorPointMap.is_null() || ownedVectorPointMap.is_null(), Teuchos::NullReferenceError,
                         "Error in PeridigmNS::DataManager::scatterToGhosts(), inconsistent vector maps.");

      // create an owned (non-overlap) multivector
      int numVectors = overlapVectorPointMultiVector->NumVectors();
      Teuchos::RCP<Epetra_MultiVector> ownedVectorPointMultiVector = Teuchos::rcp(new Epetra_MultiVector(*ownedVectorPointMap, numVectors));

      // copy data from the overlap vector into the owned (non-overlap) vector
      for(int iVec=0 ; iVec<numVectors ; ++iVec){
        double* overlapVectorPointData = (*overlapVectorPointMultiVector)[iVec];
        double* ownedVectorPointData = (*ownedVectorPointMultiVector)[iVec];
        for(int iLID=0 ; iLID<ownedVectorPointMap->NumMyElements() ; ++iLID){
          int globalID = ownedVectorPointMap->GID(iLID);
          int overlapMapLocalID = overlapVectorPointMap->LID(globalID);
          ownedVectorPointData[iLID*3]   = overlapVectorPointData[overlapMapLocalID*3];
          ownedVectorPointData[iLID*3+1] = overlapVectorPointData[overlapMapLocalID*3+1];
          ownedVectorPointData[iLID*3+2] = overlapVectorPointData[overlapMapLocalID*3+2];
        }
      }

      // scatter the data back from the owned (non-overlap) multivector into the overlap multivector
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*overlapVectorPointMap, *ownedVectorPointMap));
      overlapVectorPointMultiVector->Import(*ownedVectorPointMultiVector, *importer, Insert);
    }

    // note: bond data is not ghosted, so there's no need to scatter to ghosts.
  }
}

void PeridigmNS::DataManager::rebalance(Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedScalarPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOverlapScalarPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedVectorPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOverlapVectorPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedScalarBondMap)
{
  rebalanceCount++;

  // Rebalance involves importing from the original overlap multivectors to the new overlap multivectors.
  // Elements in the overlap (ghosted) portions of the original multivectors exist on multiple processors,
  // and there is no guarantee that the different processors hold the same values (they won't in general).
  // Therefore, prior to rebalancing, call scatterToGhosts() to get the same values on all processors.
  scatterToGhosts();

  // importers
  Teuchos::RCP<const Epetra_Import> overlapScalarPointImporter;
  if(!overlapScalarPointMap.is_null() || !rebalancedOverlapScalarPointMap.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(overlapScalarPointMap.is_null() || rebalancedOverlapScalarPointMap.is_null(), Teuchos::NullReferenceError,
                       "Error in PeridigmNS::DataManager::rebalance(), inconsistent scalar maps.");
    overlapScalarPointImporter = Teuchos::rcp(new Epetra_Import(*rebalancedOverlapScalarPointMap, *overlapScalarPointMap));
  }
  Teuchos::RCP<const Epetra_Import> overlapVectorPointImporter;
  if(!overlapVectorPointMap.is_null() || !rebalancedOverlapVectorPointMap.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(overlapVectorPointMap.is_null() || rebalancedOverlapVectorPointMap.is_null(), Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::rebalance(), inconsistent vector maps.");
    overlapVectorPointImporter = Teuchos::rcp(new Epetra_Import(*rebalancedOverlapVectorPointMap, *overlapVectorPointMap));
  }
  Teuchos::RCP<const Epetra_Import> ownedScalarBondImporter;
  if(!ownedScalarBondMap.is_null() || !rebalancedOwnedScalarBondMap.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(ownedScalarBondMap.is_null() || rebalancedOwnedScalarBondMap.is_null(), Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::rebalance(), inconsistent bond maps.");
    ownedScalarBondImporter = Teuchos::rcp(new Epetra_Import(*rebalancedOwnedScalarBondMap, *ownedScalarBondMap));
  }
    
  // state NONE
  if(statelessScalarPointFieldSpecs->size() + statelessVectorPointFieldSpecs->size() + statelessScalarBondFieldSpecs->size() > 0){
    Teuchos::RCP<State> rebalancedStateNONE = Teuchos::rcp(new State);
    if(statelessScalarPointFieldSpecs->size() > 0){
      rebalancedStateNONE->allocateScalarPointData(statelessScalarPointFieldSpecs, rebalancedOverlapScalarPointMap);
      rebalancedStateNONE->getScalarPointMultiVector()->Import(*stateNONE->getScalarPointMultiVector(), *overlapScalarPointImporter, Insert);
    }
    if(statelessVectorPointFieldSpecs->size() > 0){
      rebalancedStateNONE->allocateVectorPointData(statelessVectorPointFieldSpecs, rebalancedOverlapVectorPointMap);
      rebalancedStateNONE->getVectorPointMultiVector()->Import(*stateNONE->getVectorPointMultiVector(), *overlapVectorPointImporter, Insert);
    }
    if(statelessScalarBondFieldSpecs->size() > 0){
      rebalancedStateNONE->allocateScalarBondData(statelessScalarBondFieldSpecs, rebalancedOwnedScalarBondMap);
      rebalancedStateNONE->getScalarBondMultiVector()->Import(*stateNONE->getScalarBondMultiVector(), *ownedScalarBondImporter, Insert);
    }
    stateNONE = rebalancedStateNONE;
  }

  // states N and NP1
  if(statefulScalarPointFieldSpecs->size() + statefulVectorPointFieldSpecs->size() + statefulScalarBondFieldSpecs->size() > 0){
    Teuchos::RCP<State> rebalancedStateN = Teuchos::rcp(new State);
    if(statefulScalarPointFieldSpecs->size() > 0){
      rebalancedStateN->allocateScalarPointData(statefulScalarPointFieldSpecs, rebalancedOverlapScalarPointMap);
      rebalancedStateN->getScalarPointMultiVector()->Import(*stateN->getScalarPointMultiVector(), *overlapScalarPointImporter, Insert);
    }
    if(statefulVectorPointFieldSpecs->size() > 0){
      rebalancedStateN->allocateVectorPointData(statefulVectorPointFieldSpecs, rebalancedOverlapVectorPointMap);
      rebalancedStateN->getVectorPointMultiVector()->Import(*stateN->getVectorPointMultiVector(), *overlapVectorPointImporter, Insert);
    }
    if(statefulScalarBondFieldSpecs->size() > 0){
      rebalancedStateN->allocateScalarBondData(statefulScalarBondFieldSpecs, rebalancedOwnedScalarBondMap);
      rebalancedStateN->getScalarBondMultiVector()->Import(*stateN->getScalarBondMultiVector(), *ownedScalarBondImporter, Insert);
    }
    stateN = rebalancedStateN;
    Teuchos::RCP<State> rebalancedStateNP1 = Teuchos::rcp(new State);
    if(statefulScalarPointFieldSpecs->size() > 0){
      rebalancedStateNP1->allocateScalarPointData(statefulScalarPointFieldSpecs, rebalancedOverlapScalarPointMap);
      rebalancedStateNP1->getScalarPointMultiVector()->Import(*stateNP1->getScalarPointMultiVector(), *overlapScalarPointImporter, Insert);
    }
    if(statefulVectorPointFieldSpecs->size() > 0){
      rebalancedStateNP1->allocateVectorPointData(statefulVectorPointFieldSpecs, rebalancedOverlapVectorPointMap);
      rebalancedStateNP1->getVectorPointMultiVector()->Import(*stateNP1->getVectorPointMultiVector(), *overlapVectorPointImporter, Insert);
    }
    if(statefulScalarBondFieldSpecs->size() > 0){
      rebalancedStateNP1->allocateScalarBondData(statefulScalarBondFieldSpecs, rebalancedOwnedScalarBondMap);
      rebalancedStateNP1->getScalarBondMultiVector()->Import(*stateNP1->getScalarBondMultiVector(), *ownedScalarBondImporter, Insert);
    }
    stateNP1 = rebalancedStateNP1;
  }

  // Maps
  ownedScalarPointMap = rebalancedOwnedScalarPointMap;
  overlapScalarPointMap = rebalancedOverlapScalarPointMap;
  ownedVectorPointMap = rebalancedOwnedVectorPointMap;
  overlapVectorPointMap = rebalancedOverlapVectorPointMap;
  ownedScalarBondMap = rebalancedOwnedScalarBondMap;
}

void PeridigmNS::DataManager::copyLocallyOwnedDataFromDataManager(PeridigmNS::DataManager& source)
{
  if(!stateNONE.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(source.getStateNONE().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
    stateNONE->copyLocallyOwnedDataFromState(source.getStateNONE());
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(!source.getStateNONE().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
  }

  if(!stateN.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(source.getStateN().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
    stateN->copyLocallyOwnedDataFromState(source.getStateN());
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(!source.getStateN().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
  }

  if(!stateNP1.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(source.getStateNP1().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
    stateNP1->copyLocallyOwnedDataFromState(source.getStateNP1());
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(!source.getStateNP1().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
  }
}

bool PeridigmNS::DataManager::hasData(Field_NS::FieldSpec fieldSpec, Field_ENUM::Step step)
{
  bool hasData = false;
  if(step == Field_ENUM::STEP_NONE){
    hasData = stateNONE->hasData(fieldSpec);
  }
  else if(step == Field_ENUM::STEP_N){
    hasData = stateN->hasData(fieldSpec);
  }
  else if(step == Field_ENUM::STEP_NP1){
    hasData = stateNP1->hasData(fieldSpec);
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                       "PeridigmNS::DataManager::getData, invalid FieldStep!");
  }
  return hasData;
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::DataManager::getData(Field_NS::FieldSpec fieldSpec, Field_ENUM::Step step)
{
  Teuchos::RCP<Epetra_Vector> data;
  if(step == Field_ENUM::STEP_NONE){
    data = stateNONE->getData(fieldSpec);
  }
  else if(step == Field_ENUM::STEP_N){
    data = stateN->getData(fieldSpec);
  }
  else if(step == Field_ENUM::STEP_NP1){
    data = stateNP1->getData(fieldSpec);
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                       "PeridigmNS::DataManager::getData, invalid FieldStep!");
  }
  return data;
}
