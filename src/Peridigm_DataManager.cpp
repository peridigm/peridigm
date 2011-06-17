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

  // loop over the specs and determine:
  // 1) the number of scalar, vector2d, and vector3d fields
  // 2) the FieldType for each of the data
  // 3) whether the data has one or two states
  statelessScalarFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statelessVectorFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statelessBondFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulScalarFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulVectorFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  statefulBondFieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>());
  for(unsigned int i=0; i<fieldSpecs->size() ; ++i){
    Field_NS::FieldSpec& spec = (*fieldSpecs)[i];
    if(spec.getLength() == Field_NS::FieldSpec::SCALAR){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        statelessScalarFieldSpecs->push_back(spec);
      else
        statefulScalarFieldSpecs->push_back(spec);
    }
    else if(spec.getLength() == Field_NS::FieldSpec::VECTOR3D){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        statelessVectorFieldSpecs->push_back(spec);
      else
        statefulVectorFieldSpecs->push_back(spec);
    }
    else if(spec.getLength() == Field_NS::FieldSpec::BOND){
      if(spec.getStateArchitecture() == Field_NS::FieldSpec::STATELESS)
        statelessBondFieldSpecs->push_back(spec);
      else
        statefulBondFieldSpecs->push_back(spec);
    }
    else{
      TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                         "PeridigmNS::DataManager::allocateData, invalid FieldSpec!");
    }
  }

  // make sure maps exist before trying to create states
  if(statelessScalarFieldSpecs->size() + statefulScalarFieldSpecs->size() > 0)
    TEST_FOR_EXCEPTION(scalarMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate scalar data with no map (forget setMaps()?).");
  if(statelessVectorFieldSpecs->size() + statefulVectorFieldSpecs->size() > 0)
    TEST_FOR_EXCEPTION(vectorMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate vector data with no map (forget setMaps()?).");
  if(statelessBondFieldSpecs->size() + statefulBondFieldSpecs->size() > 0)
    TEST_FOR_EXCEPTION(bondMap == Teuchos::null, Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate bond data with no map (forget setMaps()?).");

  // create the states
  if(statelessScalarFieldSpecs->size() + statelessVectorFieldSpecs->size() + statelessBondFieldSpecs->size() > 0){
    stateNONE = Teuchos::rcp(new State);
    if(statelessScalarFieldSpecs->size() > 0)
      stateNONE->allocateScalarData(statelessScalarFieldSpecs, scalarMap);
    if(statelessVectorFieldSpecs->size() > 0)
      stateNONE->allocateVectorData(statelessVectorFieldSpecs, vectorMap);
    if(statelessBondFieldSpecs->size() > 0)
      stateNONE->allocateBondData(statelessBondFieldSpecs, bondMap);
  }
  if(statefulScalarFieldSpecs->size() + statefulVectorFieldSpecs->size() + statefulBondFieldSpecs->size() > 0){
    stateN = Teuchos::rcp(new State);
    stateNP1 = Teuchos::rcp(new State);
    if(statefulScalarFieldSpecs->size() > 0){
      stateN->allocateScalarData(statefulScalarFieldSpecs, scalarMap);
      stateNP1->allocateScalarData(statefulScalarFieldSpecs, scalarMap);
    }
    if(statefulVectorFieldSpecs->size() > 0){
      stateN->allocateVectorData(statefulVectorFieldSpecs, vectorMap);
      stateNP1->allocateVectorData(statefulVectorFieldSpecs, vectorMap);
    }   
    if(statefulBondFieldSpecs->size() > 0){
      stateN->allocateBondData(statefulBondFieldSpecs, bondMap);
      stateNP1->allocateBondData(statefulBondFieldSpecs, bondMap);
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
    Teuchos::RCP<Epetra_MultiVector> overlapMultiVector = state->getScalarMultiVector();
    if(!overlapMultiVector.is_null()){

      TEST_FOR_EXCEPTION(scalarMap.is_null() || ownedIDScalarMap.is_null(), Teuchos::NullReferenceError,
                         "Error in PeridigmNS::DataManager::scatterToGhosts(), inconsistent scalar maps.");

      // create a non-overlap multivector
      int numVectors = overlapMultiVector->NumVectors();
      Teuchos::RCP<const Epetra_BlockMap> overlapMap = scalarMap;
      Teuchos::RCP<const Epetra_BlockMap> nonOverlapMap = ownedIDScalarMap;
      Teuchos::RCP<Epetra_MultiVector> nonOverlapMultiVector = Teuchos::rcp(new Epetra_MultiVector(*nonOverlapMap, numVectors));

      // copy data from the overlap vector into the non-overlap vector
      for(int iVec=0 ; iVec<numVectors ; ++iVec){
        double* overlapData = (*overlapMultiVector)[iVec];
        double* nonOverlapData = (*nonOverlapMultiVector)[iVec];
        for(int iLID=0 ; iLID<nonOverlapMap->NumMyElements() ; ++iLID){
          int globalID = nonOverlapMap->GID(iLID);
          int overlapMapLocalID = overlapMap->LID(globalID);
          nonOverlapData[iLID] = overlapData[overlapMapLocalID];
        }
      }

      // scatter the data back from the non-overlap multivector into the overlap multivector
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*overlapMap, *nonOverlapMap));
      overlapMultiVector->Import(*nonOverlapMultiVector, *importer, Insert);
    }

    // process vector data
    overlapMultiVector = state->getVectorMultiVector();
    if(!overlapMultiVector.is_null()){

      TEST_FOR_EXCEPTION(vectorMap.is_null() || ownedIDVectorMap.is_null(), Teuchos::NullReferenceError,
                         "Error in PeridigmNS::DataManager::scatterToGhosts(), inconsistent vector maps.");

      // create a non-overlap multivector
      int numVectors = overlapMultiVector->NumVectors();
      Teuchos::RCP<const Epetra_BlockMap> overlapMap = vectorMap;
      Teuchos::RCP<const Epetra_BlockMap> nonOverlapMap = ownedIDVectorMap;
      Teuchos::RCP<Epetra_MultiVector> nonOverlapMultiVector = Teuchos::rcp(new Epetra_MultiVector(*nonOverlapMap, numVectors));

      // copy data from the overlap vector into the non-overlap vector
      for(int iVec=0 ; iVec<numVectors ; ++iVec){
        double* overlapData = (*overlapMultiVector)[iVec];
        double* nonOverlapData = (*nonOverlapMultiVector)[iVec];
        for(int iLID=0 ; iLID<nonOverlapMap->NumMyElements() ; ++iLID){
          int globalID = nonOverlapMap->GID(iLID);
          int overlapMapLocalID = overlapMap->LID(globalID);
          nonOverlapData[iLID*3] = overlapData[overlapMapLocalID*3];
          nonOverlapData[iLID*3+1] = overlapData[overlapMapLocalID*3+1];
          nonOverlapData[iLID*3+2] = overlapData[overlapMapLocalID*3+2];
        }
      }

      // scatter the data back from the non-overlap multivector into the overlap multivector
      Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(*overlapMap, *nonOverlapMap));
      overlapMultiVector->Import(*nonOverlapMultiVector, *importer, Insert);
    }

    // note: bond data is not ghosted, so there's no need to scatter to ghosts.
  }
}

void PeridigmNS::DataManager::rebalance(Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedIDScalarMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedIDVectorMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedScalarMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedVectorMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap)
{
  rebalanceCount++;

  // Rebalance involves importing from the original overlap multivectors to the new overlap multivectors.
  // Elements in the overlap (ghosted) portions of the original multivectors exist on multiple processors,
  // and there is no guarantee that the different processors hold the same values (they won't in general).
  // Therefore, prior to rebalancing, call scatterToGhosts() to get the same values on all processors.
  scatterToGhosts();

  // importers
  Teuchos::RCP<const Epetra_Import> scalarImporter;
  if(!scalarMap.is_null() || !rebalancedScalarMap.is_null()){
    TEST_FOR_EXCEPTION(scalarMap.is_null() || rebalancedScalarMap.is_null(), Teuchos::NullReferenceError,
                       "Error in PeridigmNS::DataManager::rebalance(), inconsistent scalar maps.");
    scalarImporter = Teuchos::rcp(new Epetra_Import(*rebalancedScalarMap, *scalarMap));
  }
  Teuchos::RCP<const Epetra_Import> vectorImporter;
  if(!vectorMap.is_null() || !rebalancedVectorMap.is_null()){
    TEST_FOR_EXCEPTION(vectorMap.is_null() || rebalancedVectorMap.is_null(), Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::rebalance(), inconsistent vector maps.");
    vectorImporter = Teuchos::rcp(new Epetra_Import(*rebalancedVectorMap, *vectorMap));
  }
  Teuchos::RCP<const Epetra_Import> bondImporter;
  if(!bondMap.is_null() || !rebalancedBondMap.is_null()){
    TEST_FOR_EXCEPTION(bondMap.is_null() || rebalancedBondMap.is_null(), Teuchos::NullReferenceError, 
                       "Error in PeridigmNS::DataManager::rebalance(), inconsistent bond maps.");
    bondImporter = Teuchos::rcp(new Epetra_Import(*rebalancedBondMap, *bondMap));
  }
    
  // state NONE
  if(statelessScalarFieldSpecs->size() + statelessVectorFieldSpecs->size() + statelessBondFieldSpecs->size() > 0){
    Teuchos::RCP<State> rebalancedStateNONE = Teuchos::rcp(new State);
    if(statelessScalarFieldSpecs->size() > 0){
      rebalancedStateNONE->allocateScalarData(statelessScalarFieldSpecs, rebalancedScalarMap);
      rebalancedStateNONE->getScalarMultiVector()->Import(*stateNONE->getScalarMultiVector(), *scalarImporter, Insert);
    }
    if(statelessVectorFieldSpecs->size() > 0){
      rebalancedStateNONE->allocateVectorData(statelessVectorFieldSpecs, rebalancedVectorMap);
      rebalancedStateNONE->getVectorMultiVector()->Import(*stateNONE->getVectorMultiVector(), *vectorImporter, Insert);
    }
    if(statelessBondFieldSpecs->size() > 0){
      rebalancedStateNONE->allocateBondData(statelessBondFieldSpecs, rebalancedBondMap);
      rebalancedStateNONE->getBondMultiVector()->Import(*stateNONE->getBondMultiVector(), *bondImporter, Insert);
    }
    stateNONE = rebalancedStateNONE;
  }

  // states N and NP1
  if(statefulScalarFieldSpecs->size() + statefulVectorFieldSpecs->size() + statefulBondFieldSpecs->size() > 0){
    Teuchos::RCP<State> rebalancedStateN = Teuchos::rcp(new State);
    if(statefulScalarFieldSpecs->size() > 0){
      rebalancedStateN->allocateScalarData(statefulScalarFieldSpecs, rebalancedScalarMap);
      rebalancedStateN->getScalarMultiVector()->Import(*stateN->getScalarMultiVector(), *scalarImporter, Insert);
    }
    if(statefulVectorFieldSpecs->size() > 0){
      rebalancedStateN->allocateVectorData(statefulVectorFieldSpecs, rebalancedVectorMap);
      rebalancedStateN->getVectorMultiVector()->Import(*stateN->getVectorMultiVector(), *vectorImporter, Insert);
    }
    if(statefulBondFieldSpecs->size() > 0){
      rebalancedStateN->allocateBondData(statefulBondFieldSpecs, rebalancedBondMap);
      rebalancedStateN->getBondMultiVector()->Import(*stateN->getBondMultiVector(), *bondImporter, Insert);
    }
    stateN = rebalancedStateN;
    Teuchos::RCP<State> rebalancedStateNP1 = Teuchos::rcp(new State);
    if(statefulScalarFieldSpecs->size() > 0){
      rebalancedStateNP1->allocateScalarData(statefulScalarFieldSpecs, rebalancedScalarMap);
      rebalancedStateNP1->getScalarMultiVector()->Import(*stateNP1->getScalarMultiVector(), *scalarImporter, Insert);
    }
    if(statefulVectorFieldSpecs->size() > 0){
      rebalancedStateNP1->allocateVectorData(statefulVectorFieldSpecs, rebalancedVectorMap);
      rebalancedStateNP1->getVectorMultiVector()->Import(*stateNP1->getVectorMultiVector(), *vectorImporter, Insert);
    }
    if(statefulBondFieldSpecs->size() > 0){
      rebalancedStateNP1->allocateBondData(statefulBondFieldSpecs, rebalancedBondMap);
      rebalancedStateNP1->getBondMultiVector()->Import(*stateNP1->getBondMultiVector(), *bondImporter, Insert);
    }
    stateNP1 = rebalancedStateNP1;
  }

  // Maps
  ownedIDScalarMap = rebalancedOwnedIDScalarMap;
  ownedIDVectorMap = rebalancedOwnedIDVectorMap;
  scalarMap = rebalancedScalarMap;
  vectorMap = rebalancedVectorMap;
  bondMap = rebalancedBondMap;
}

void PeridigmNS::DataManager::copyLocallyOwnedDataFromDataManager(PeridigmNS::DataManager& source)
{
  if(!stateNONE.is_null()){
    TEST_FOR_EXCEPTION(source.getStateNONE().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
    stateNONE->copyLocallyOwnedDataFromState(source.getStateNONE());
  }
  else{
    TEST_FOR_EXCEPTION(!source.getStateNONE().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
  }

  if(!stateN.is_null()){
    TEST_FOR_EXCEPTION(source.getStateN().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
    stateN->copyLocallyOwnedDataFromState(source.getStateN());
  }
  else{
    TEST_FOR_EXCEPTION(!source.getStateN().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
  }

  if(!stateNP1.is_null()){
    TEST_FOR_EXCEPTION(source.getStateNP1().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
    stateNP1->copyLocallyOwnedDataFromState(source.getStateNP1());
  }
  else{
    TEST_FOR_EXCEPTION(!source.getStateNP1().is_null(), Teuchos::NullReferenceError, "PeridigmNS::State::copyLocallyOwnedDataFromDataManager() called with incompatible source and target.\n");
  }
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::DataManager::getData(Field_NS::FieldSpec fieldSpec, Field_NS::FieldSpec::FieldStep fieldStep)
{
  Teuchos::RCP<Epetra_Vector> data;
  if(fieldStep == Field_NS::FieldSpec::STEP_NONE){
    data = stateNONE->getData(fieldSpec);
  }
  else if(fieldStep == Field_NS::FieldSpec::STEP_N){
    data = stateN->getData(fieldSpec);
  }
  else if(fieldStep == Field_NS::FieldSpec::STEP_NP1){
    data = stateNP1->getData(fieldSpec);
  }
  else{
    TEST_FOR_EXCEPTION(false, Teuchos::RangeError, 
                       "PeridigmNS::DataManager::getData, invalid FieldStep!");
  }
  return data;
}
