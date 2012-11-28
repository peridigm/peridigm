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
#include "Peridigm_Field.hpp"

using namespace std;

void PeridigmNS::DataManager::allocateData(vector<int> fieldIds)
{
  // remove duplicates
  sort(fieldIds.begin(), fieldIds.end());
  vector<int>::iterator newEnd = unique(fieldIds.begin(), fieldIds.end());
  fieldIds.erase(newEnd, fieldIds.end());

  FieldManager& fieldManager = FieldManager::self();

  for(unsigned int i=0; i<fieldIds.size() ; ++i){

    int fieldId = fieldIds[i];
    PeridigmNS::FieldSpec spec = fieldManager.getFieldSpec(fieldId);

    // scalar point data
    if(spec.getLength() == PeridigmField::SCALAR && spec.getRelation() == PeridigmField::ELEMENT){
      if(spec.getTemporal() == PeridigmField::CONSTANT)
        statelessScalarPointFieldIds.push_back(fieldId);
      else
        statefulScalarPointFieldIds.push_back(fieldId);
    }
    // vector point data
    else if(spec.getLength() == PeridigmField::VECTOR && spec.getRelation() == PeridigmField::NODE){
      if(spec.getTemporal() == PeridigmField::CONSTANT)
        statelessVectorPointFieldIds.push_back(fieldId);
      else
        statefulVectorPointFieldIds.push_back(fieldId);
    }
    // scalar bond data
    else if(spec.getLength() == PeridigmField::SCALAR && spec.getRelation() == PeridigmField::BOND){
      if(spec.getTemporal() == PeridigmField::CONSTANT)
        statelessScalarBondFieldIds.push_back(fieldId);
      else
        statefulScalarBondFieldIds.push_back(fieldId);
    }
    // ignore global data
    else if (spec.getRelation() == PeridigmField::GLOBAL)
      continue;
    else{
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::RangeError, 
                                 "PeridigmNS::DataManager::allocateData, invalid FieldSpec!");
    }

    allFieldIds.push_back(fieldId);
  }

  // make sure maps exist before trying to create states
  if(statelessScalarPointFieldIds.size() + statefulScalarPointFieldIds.size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(overlapScalarPointMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate scalar data with no map (forget setMaps()?).");
  if(statelessVectorPointFieldIds.size() + statefulVectorPointFieldIds.size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(overlapVectorPointMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate vector data with no map (forget setMaps()?).");
  if(statelessScalarBondFieldIds.size() + statefulScalarBondFieldIds.size() > 0)
    TEUCHOS_TEST_FOR_EXCEPTION(ownedScalarBondMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate bond data with no map (forget setMaps()?).");

  // create the states
  if(statelessScalarPointFieldIds.size() + statelessVectorPointFieldIds.size() + statelessScalarBondFieldIds.size() > 0){
    stateNONE = Teuchos::rcp(new State);
    if(statelessScalarPointFieldIds.size() > 0)
      stateNONE->allocateScalarPointData(statelessScalarPointFieldIds, overlapScalarPointMap);
    if(statelessVectorPointFieldIds.size() > 0)
      stateNONE->allocateVectorPointData(statelessVectorPointFieldIds, overlapVectorPointMap);
    if(statelessScalarBondFieldIds.size() > 0)
      stateNONE->allocateScalarBondData(statelessScalarBondFieldIds, ownedScalarBondMap);
  }
  if(statefulScalarPointFieldIds.size() + statefulVectorPointFieldIds.size() + statefulScalarBondFieldIds.size() > 0){
    stateN = Teuchos::rcp(new State);
    stateNP1 = Teuchos::rcp(new State);
    if(statefulScalarPointFieldIds.size() > 0){
      stateN->allocateScalarPointData(statefulScalarPointFieldIds, overlapScalarPointMap);
      stateNP1->allocateScalarPointData(statefulScalarPointFieldIds, overlapScalarPointMap);
    }
    if(statefulVectorPointFieldIds.size() > 0){
      stateN->allocateVectorPointData(statefulVectorPointFieldIds, overlapVectorPointMap);
      stateNP1->allocateVectorPointData(statefulVectorPointFieldIds, overlapVectorPointMap);
    }   
    if(statefulScalarBondFieldIds.size() > 0){
      stateN->allocateScalarBondData(statefulScalarBondFieldIds, ownedScalarBondMap);
      stateNP1->allocateScalarBondData(statefulScalarBondFieldIds, ownedScalarBondMap);
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
  if(statelessScalarPointFieldIds.size() + statelessVectorPointFieldIds.size() + statelessScalarBondFieldIds.size() > 0){
    Teuchos::RCP<State> rebalancedStateNONE = Teuchos::rcp(new State);
    if(statelessScalarPointFieldIds.size() > 0){
      rebalancedStateNONE->allocateScalarPointData(statelessScalarPointFieldIds, rebalancedOverlapScalarPointMap);
      rebalancedStateNONE->getScalarPointMultiVector()->Import(*stateNONE->getScalarPointMultiVector(), *overlapScalarPointImporter, Insert);
    }
    if(statelessVectorPointFieldIds.size() > 0){
      rebalancedStateNONE->allocateVectorPointData(statelessVectorPointFieldIds, rebalancedOverlapVectorPointMap);
      rebalancedStateNONE->getVectorPointMultiVector()->Import(*stateNONE->getVectorPointMultiVector(), *overlapVectorPointImporter, Insert);
    }
    if(statelessScalarBondFieldIds.size() > 0){
      rebalancedStateNONE->allocateScalarBondData(statelessScalarBondFieldIds, rebalancedOwnedScalarBondMap);
      rebalancedStateNONE->getScalarBondMultiVector()->Import(*stateNONE->getScalarBondMultiVector(), *ownedScalarBondImporter, Insert);
    }
    stateNONE = rebalancedStateNONE;
  }

  // states N and NP1
  if(statefulScalarPointFieldIds.size() + statefulVectorPointFieldIds.size() + statefulScalarBondFieldIds.size() > 0){
    Teuchos::RCP<State> rebalancedStateN = Teuchos::rcp(new State);
    if(statefulScalarPointFieldIds.size() > 0){
      rebalancedStateN->allocateScalarPointData(statefulScalarPointFieldIds, rebalancedOverlapScalarPointMap);
      rebalancedStateN->getScalarPointMultiVector()->Import(*stateN->getScalarPointMultiVector(), *overlapScalarPointImporter, Insert);
    }
    if(statefulVectorPointFieldIds.size() > 0){
      rebalancedStateN->allocateVectorPointData(statefulVectorPointFieldIds, rebalancedOverlapVectorPointMap);
      rebalancedStateN->getVectorPointMultiVector()->Import(*stateN->getVectorPointMultiVector(), *overlapVectorPointImporter, Insert);
    }
    if(statefulScalarBondFieldIds.size() > 0){
      rebalancedStateN->allocateScalarBondData(statefulScalarBondFieldIds, rebalancedOwnedScalarBondMap);
      rebalancedStateN->getScalarBondMultiVector()->Import(*stateN->getScalarBondMultiVector(), *ownedScalarBondImporter, Insert);
    }
    stateN = rebalancedStateN;
    Teuchos::RCP<State> rebalancedStateNP1 = Teuchos::rcp(new State);
    if(statefulScalarPointFieldIds.size() > 0){
      rebalancedStateNP1->allocateScalarPointData(statefulScalarPointFieldIds, rebalancedOverlapScalarPointMap);
      rebalancedStateNP1->getScalarPointMultiVector()->Import(*stateNP1->getScalarPointMultiVector(), *overlapScalarPointImporter, Insert);
    }
    if(statefulVectorPointFieldIds.size() > 0){
      rebalancedStateNP1->allocateVectorPointData(statefulVectorPointFieldIds, rebalancedOverlapVectorPointMap);
      rebalancedStateNP1->getVectorPointMultiVector()->Import(*stateNP1->getVectorPointMultiVector(), *overlapVectorPointImporter, Insert);
    }
    if(statefulScalarBondFieldIds.size() > 0){
      rebalancedStateNP1->allocateScalarBondData(statefulScalarBondFieldIds, rebalancedOwnedScalarBondMap);
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

bool PeridigmNS::DataManager::hasData(int fieldId, PeridigmField::Step step)
{
  bool hasData = false;
  if(step == PeridigmField::STEP_NONE){
    hasData = stateNONE->hasData(fieldId);
  }
  else if(step == PeridigmField::STEP_N){
    hasData = stateN->hasData(fieldId);
  }
  else if(step == PeridigmField::STEP_NP1){
    hasData = stateNP1->hasData(fieldId);
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                       "PeridigmNS::DataManager::getData, invalid fieldId and step!");
  }
  return hasData;
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::DataManager::getData(int fieldId, PeridigmField::Step step)
{
  Teuchos::RCP<Epetra_Vector> data;
  if(step == PeridigmField::STEP_NONE){
    data = stateNONE->getData(fieldId);
  }
  else if(step == PeridigmField::STEP_N){
    data = stateN->getData(fieldId);
  }
  else if(step == PeridigmField::STEP_NP1){
    data = stateNP1->getData(fieldId);
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                       "PeridigmNS::DataManager::getData, invalid fieldId and step!");
  }
  return data;
}
