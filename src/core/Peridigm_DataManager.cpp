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
#include <Epetra_Comm.h>
#include "Peridigm_DataManager.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

std::vector<int> PeridigmNS::DataManager::allGlobalFieldIds;
Teuchos::RCP<const Epetra_BlockMap> PeridigmNS::DataManager::scalarGlobalMap;
Teuchos::RCP<const Epetra_BlockMap> PeridigmNS::DataManager::vectorGlobalMap;
std::map< std::pair<int, PeridigmNS::PeridigmField::Step>, Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::fieldIdAndStepToGlobalData;
std::vector< Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::scalarGlobalDataStateN;
std::vector< Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::scalarGlobalDataStateNP1;
std::vector< Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::scalarGlobalDataStateNONE;
std::vector< Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::vectorGlobalDataStateN;
std::vector< Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::vectorGlobalDataStateNP1;
std::vector< Teuchos::RCP<Epetra_Vector> > PeridigmNS::DataManager::vectorGlobalDataStateNONE;

void PeridigmNS::DataManager::allocateData(vector<int> fieldIds)
{
  // remove duplicates
  sort(fieldIds.begin(), fieldIds.end());
  vector<int>::iterator newEnd = unique(fieldIds.begin(), fieldIds.end());
  fieldIds.erase(newEnd, fieldIds.end());

  // temp vectors for global data bookkeeping
  vector<int> statelessScalarGlobalFieldIds;
  vector<int> statelessVectorGlobalFieldIds;
  vector<int> statefulScalarGlobalFieldIds;
  vector<int> statefulVectorGlobalFieldIds;

  for(unsigned int i=0; i<fieldIds.size() ; ++i){

    int fieldId = fieldIds[i];
    PeridigmNS::FieldSpec spec = fieldManager.getFieldSpec(fieldId);
    PeridigmField::Relation relation = spec.getRelation();
    PeridigmField::Length length = spec.getLength();
    PeridigmField::Temporal temporal = spec.getTemporal();

    // Global data
    if(relation == PeridigmField::GLOBAL){

      // Allocate global data only if it hasn't already been allocated
      if( find(allGlobalFieldIds.begin(), allGlobalFieldIds.end(), fieldId) == allGlobalFieldIds.end() ){

        if(length == PeridigmField::SCALAR && temporal == PeridigmField::CONSTANT)
          statelessScalarGlobalFieldIds.push_back(fieldId);

        else if(length == PeridigmField::SCALAR && temporal == PeridigmField::TWO_STEP)
          statefulScalarGlobalFieldIds.push_back(fieldId);

        else if(length == PeridigmField::VECTOR && temporal == PeridigmField::CONSTANT)
          statelessVectorGlobalFieldIds.push_back(fieldId);

        else if(length == PeridigmField::VECTOR && temporal == PeridigmField::TWO_STEP)
          statefulVectorGlobalFieldIds.push_back(fieldId);

        else
          TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::RangeError, 
                                     "PeridigmNS::DataManager::allocateData, invalid FieldSpec!");

        allGlobalFieldIds.push_back(fieldId);
      }
    }
    // Bond data
    else if(relation == PeridigmField::BOND){

      // Vector bond data is not supported
      TEUCHOS_TEST_FOR_EXCEPTION(length != PeridigmField::SCALAR, Teuchos::RangeError, 
                                 "PeridigmNS::DataManager::allocateData, invalid FieldSpec, BOND data must be SCALAR!");

      if(temporal == PeridigmField::CONSTANT)
        statelessBondFieldIds.push_back(fieldId);

      else if(temporal == PeridigmField::TWO_STEP)
        statefulBondFieldIds.push_back(fieldId);

      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::RangeError, 
                                   "PeridigmNS::DataManager::allocateData, invalid FieldSpec!");

      allFieldIds.push_back(fieldId);
    }
    // Element and node data
    else if(relation == PeridigmField::NODE || relation == PeridigmField::ELEMENT){

      if(temporal == PeridigmField::CONSTANT)
        statelessPointFieldIds[length].push_back(fieldId);

      else if(temporal == PeridigmField::TWO_STEP)
        statefulPointFieldIds[length].push_back(fieldId);

      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::RangeError, 
                                   "PeridigmNS::DataManager::allocateData, invalid FieldSpec!");

      allFieldIds.push_back(fieldId);
    }
  }

  // add global field ids into list of all field ids
  allFieldIds.insert(allFieldIds.end(), allGlobalFieldIds.begin(), allGlobalFieldIds.end());

  int numGlobalElements(0), numMyElements(0), *myGlobalElements(0), indexBase(0);

  // make sure maps exist before trying to create states
  if(statelessPointFieldIds.size() + statefulPointFieldIds.size() > 0){
    TEUCHOS_TEST_FOR_EXCEPTION(overlapScalarPointMap.is_null(), Teuchos::NullReferenceError, 
                               "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate point data with no map (forget setMaps()?).");
    numGlobalElements = overlapScalarPointMap->NumGlobalElements();
    numMyElements = overlapScalarPointMap->NumMyElements();
    myGlobalElements = overlapScalarPointMap->MyGlobalElements();
  }
  if(statelessBondFieldIds.size() + statefulBondFieldIds.size() > 0){
    TEUCHOS_TEST_FOR_EXCEPTION(ownedBondMap.is_null(), Teuchos::NullReferenceError, 
                               "Error in PeridigmNS::DataManager::allocateData(), attempting to allocate bond data with no map (forget setMaps()?).");
  }

  // create the global data (these are static members of DataManager)
  if(statelessScalarGlobalFieldIds.size() + statefulScalarGlobalFieldIds.size() + statelessVectorGlobalFieldIds.size() + statefulVectorGlobalFieldIds.size()> 0){

    // obtain a comm object for use in creating the Epetra_BlockMaps
    Teuchos::RCP<const Epetra_Comm> comm = getEpetraComm();
    
    // create the maps for global data, if needed
    if(scalarGlobalMap.is_null()){
      int tempNumGlobalElements = 1;
      int tempNumMyElements = 1;
      vector<int> tempMyGlobalElements(1);
      tempMyGlobalElements[0] = 0;
      int tempElementSize = 1;
      int tempIndexBase = 0;
      scalarGlobalMap = Teuchos::RCP<Epetra_BlockMap>(new Epetra_BlockMap(tempNumGlobalElements,
                                                                          tempNumMyElements,
                                                                          &tempMyGlobalElements[0],
                                                                          tempElementSize,
                                                                          tempIndexBase,
                                                                          *comm));
      tempElementSize = 3;
      vectorGlobalMap = Teuchos::RCP<Epetra_BlockMap>(new Epetra_BlockMap(tempNumGlobalElements,
                                                                          tempNumMyElements,
                                                                          &tempMyGlobalElements[0],
                                                                          tempElementSize,
                                                                          tempIndexBase,
                                                                          *comm));
    }

    // Allocate global data and set up the FieldId-to-Epetra_Vector map
    for(unsigned int i=0 ; i<statelessScalarGlobalFieldIds.size() ; ++i){
      int index = scalarGlobalDataStateNONE.size();
      scalarGlobalDataStateNONE.push_back( Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*scalarGlobalMap)) );
      fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(statelessScalarGlobalFieldIds[i], PeridigmField::STEP_NONE) ] = scalarGlobalDataStateNONE[index];
    }
    for(unsigned int i=0 ; i<statefulScalarGlobalFieldIds.size() ; ++i){
      int index = scalarGlobalDataStateN.size();
      scalarGlobalDataStateN.push_back( Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*scalarGlobalMap)) );
      fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(statefulScalarGlobalFieldIds[i], PeridigmField::STEP_N) ] = scalarGlobalDataStateN[index];
      scalarGlobalDataStateNP1.push_back( Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*scalarGlobalMap)) );
      fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(statefulScalarGlobalFieldIds[i], PeridigmField::STEP_NP1) ] = scalarGlobalDataStateNP1[index];
    }
    for(unsigned int i=0 ; i<statelessVectorGlobalFieldIds.size() ; ++i){
      int index = vectorGlobalDataStateNONE.size();
      vectorGlobalDataStateNONE.push_back( Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*vectorGlobalMap)) );
      fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(statelessVectorGlobalFieldIds[i], PeridigmField::STEP_NONE) ] = vectorGlobalDataStateNONE[index];
    }
    for(unsigned int i=0 ; i<statefulVectorGlobalFieldIds.size() ; ++i){
      int index = vectorGlobalDataStateN.size();
      vectorGlobalDataStateN.push_back( Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*vectorGlobalMap)) );
      fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(statefulVectorGlobalFieldIds[i], PeridigmField::STEP_N) ] = vectorGlobalDataStateN[index];
      vectorGlobalDataStateNP1.push_back( Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*vectorGlobalMap)) );
      fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(statefulVectorGlobalFieldIds[i], PeridigmField::STEP_NP1) ] = vectorGlobalDataStateNP1[index];
    }
  }

  map< PeridigmField::Length, vector<int> >::iterator it;
  Teuchos::RCP<const Epetra_BlockMap> map;

  // create the states
  if(statelessPointFieldIds.size() + statelessBondFieldIds.size() > 0){
    stateNONE = Teuchos::rcp(new State);
    for(it = statelessPointFieldIds.begin() ; it != statelessPointFieldIds.end() ; ++it){
      PeridigmField::Length length = it->first;
      vector<int>& fieldIds = it->second;
      if(length == PeridigmField::SCALAR)
        map = overlapScalarPointMap;
      else if(length == PeridigmField::VECTOR)
        map = overlapVectorPointMap;
      else
        map = Teuchos::RCP<Epetra_BlockMap>(new Epetra_BlockMap(numGlobalElements,
                                                                numMyElements,
                                                                myGlobalElements,
                                                                PeridigmField::variableDimension(length),
                                                                indexBase,
                                                                *getEpetraComm()));
      stateNONE->allocatePointData(length, fieldIds, map);
    }
    if(statelessBondFieldIds.size() > 0){
      stateNONE->allocateBondData(statelessBondFieldIds, ownedBondMap);
    }
  }
  if(statefulPointFieldIds.size() + statefulBondFieldIds.size() > 0){
    stateN = Teuchos::rcp(new State);
    stateNP1 = Teuchos::rcp(new State);
    for(it = statefulPointFieldIds.begin() ; it != statefulPointFieldIds.end() ; ++it){
      PeridigmField::Length length = it->first;
      vector<int>& fieldIds = it->second;
      if(length == PeridigmField::SCALAR)
        map = overlapScalarPointMap;
      else if(length == PeridigmField::VECTOR)
        map = overlapVectorPointMap;
      else
        map = Teuchos::RCP<Epetra_BlockMap>(new Epetra_BlockMap(numGlobalElements,
                                                                numMyElements,
                                                                myGlobalElements,
                                                                PeridigmField::variableDimension(length),
                                                                indexBase,
                                                                *getEpetraComm()));
      stateN->allocatePointData(length, fieldIds, map);
      stateNP1->allocatePointData(length, fieldIds, map);
    }
    if(statefulBondFieldIds.size() > 0){
      stateN->allocateBondData(statefulBondFieldIds, ownedBondMap);
      stateNP1->allocateBondData(statefulBondFieldIds, ownedBondMap);
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

  // Store information on owned map
  int numGlobalElements = ownedScalarPointMap->NumGlobalElements();
  int numMyElements = ownedScalarPointMap->NumMyElements();
  int* myGlobalElements = ownedScalarPointMap->MyGlobalElements();
  int indexBase(0);

  // Loop over the states and scatter to the ghosted points
  for(int iState=0 ; iState<3 ; ++iState){

    Teuchos::RCP<State> state = stateNONE;
    if(iState == 1)
      state = stateN;
    else if(iState == 2)
      state = stateNP1;

    for(int iMultiVector = 0 ; iMultiVector < state->getMaxPointDataElementSize() ; ++iMultiVector){

      PeridigmField::Length length;      
      switch(iMultiVector){
      case 0: length = PeridigmField::LENGTH_1;
        break;
      case 1: length = PeridigmField::LENGTH_2;
        break;
      case 2: length = PeridigmField::LENGTH_3;
        break;
      case 3: length = PeridigmField::LENGTH_4;
        break;
      case 4: length = PeridigmField::LENGTH_5;
        break;
      case 5: length = PeridigmField::LENGTH_6;
        break;
      case 6: length = PeridigmField::LENGTH_7;
        break;
      case 7: length = PeridigmField::LENGTH_8;
        break;
      case 8: length = PeridigmField::LENGTH_9;
        break;
      default: TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::RangeError, "\n****Error, invalid PeridigmField::Length.\n");
        break;
      }
      int elementSize = PeridigmField::variableDimension(length);
        
      Teuchos::RCP<Epetra_MultiVector> overlapPointMultiVector = state->getPointMultiVector(length);
      if(!overlapPointMultiVector.is_null()){

        const Epetra_BlockMap& overlapMap = overlapPointMultiVector->Map();

        Teuchos::RCP<const Epetra_BlockMap> ownedMap;
        if(length == PeridigmField::LENGTH_1)
          ownedMap = ownedScalarPointMap;
        else if(length == PeridigmField::LENGTH_3)
          ownedMap = ownedVectorPointMap;
        else{
          ownedMap = Teuchos::RCP<Epetra_BlockMap>(new Epetra_BlockMap(numGlobalElements,
                                                                       numMyElements,
                                                                       myGlobalElements,
                                                                       PeridigmField::variableDimension(length),
                                                                       indexBase,
                                                                       overlapMap.Comm()));
        }

        // create an owned (non-overlap) multivector
        int numVectors = overlapPointMultiVector->NumVectors();
        Teuchos::RCP<Epetra_MultiVector> ownedPointMultiVector = Teuchos::rcp(new Epetra_MultiVector(*ownedMap, numVectors));

        // copy data from the overlap vector into the owned (non-overlap) vector
        for(int iVec=0 ; iVec<numVectors ; ++iVec){
          double* overlapPointData = (*overlapPointMultiVector)[iVec];
          double* ownedPointData = (*ownedPointMultiVector)[iVec];
          for(int iLID=0 ; iLID<ownedMap->NumMyElements() ; ++iLID){
            int globalID = ownedMap->GID(iLID);
            int overlapMapLocalID = overlapMap.LID(globalID);
            for(int i=0 ; i<elementSize ; ++i)
              ownedPointData[iLID*elementSize+i] = overlapPointData[overlapMapLocalID*elementSize+i];
          }
        }

        // scatter the data back from the owned (non-overlap) multivector into the overlap multivector
        Teuchos::RCP<Epetra_Import> importer = Teuchos::rcp(new Epetra_Import(overlapMap, *ownedMap));
        overlapPointMultiVector->Import(*ownedPointMultiVector, *importer, Insert);
      }
    }

    // note: global data and bond data are not ghosted, so there's no need to scatter to ghosts.
  }
}

void PeridigmNS::DataManager::rebalance(Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedScalarPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOverlapScalarPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedVectorPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOverlapVectorPointMap,
                                        Teuchos::RCP<const Epetra_BlockMap> rebalancedOwnedBondMap)
{
  rebalanceCount++;

  // Rebalance involves importing from the original overlap multivectors to the new overlap multivectors.
  // Elements in the overlap (ghosted) portions of the original multivectors exist on multiple processors,
  // and there is no guarantee that the different processors hold the same values (they won't in general).
  // Therefore, prior to rebalancing, call scatterToGhosts() to get the same values on all processors.
  scatterToGhosts();

  // Store information used in the creation of rebalanced overlap maps
  int numGlobalElements = rebalancedOverlapScalarPointMap->NumGlobalElements();
  int numMyElements = rebalancedOverlapScalarPointMap->NumMyElements();
  int* myGlobalElements = rebalancedOverlapScalarPointMap->MyGlobalElements();
  int indexBase(0);
  Teuchos::RCP<const Epetra_Comm> comm = getEpetraComm();

  map< PeridigmField::Length, vector<int> >::iterator it;
  Teuchos::RCP<const Epetra_BlockMap> map;
    
  for(int iState=0 ; iState<3 ; ++iState){

    Teuchos::RCP<State> state;
    std::map< PeridigmField::Length, vector<int> > *pointFieldIds;
    vector<int> *bondFieldIds;
    if(iState == 0){
      state = stateNONE;
      pointFieldIds = &statelessPointFieldIds;
      bondFieldIds = &statelessBondFieldIds;
    }
    else if(iState == 1){
      state = stateN;
      pointFieldIds = &statefulPointFieldIds;
      bondFieldIds = &statefulBondFieldIds;
    }
    else if(iState == 2){
      state = stateNP1;
      pointFieldIds = &statefulPointFieldIds;
      bondFieldIds = &statefulBondFieldIds;
    }

    if(!state.is_null()){

      Teuchos::RCP<State> rebalancedState = Teuchos::rcp(new State);

      // Allocate point-wise data and import from the old State to the rebalanced State
      for(it = pointFieldIds->begin() ; it != pointFieldIds->end() ; ++it){
        PeridigmField::Length length = it->first;
        vector<int>& fieldIds = it->second;
        if(length == PeridigmField::SCALAR)
          map = rebalancedOverlapScalarPointMap;
        else if(length == PeridigmField::VECTOR)
          map = rebalancedOverlapVectorPointMap;
        else
          map = Teuchos::RCP<Epetra_BlockMap>(new Epetra_BlockMap(numGlobalElements,
                                                                  numMyElements,
                                                                  myGlobalElements,
                                                                  PeridigmField::variableDimension(length),
                                                                  indexBase,
                                                                  *getEpetraComm()));
        rebalancedState->allocatePointData(length, fieldIds, map);
        Epetra_Import importer(*map, state->getPointMultiVector(length)->Map());
        rebalancedState->getPointMultiVector(length)->Import(*state->getPointMultiVector(length), importer, Insert);
      }

      // Allocate bond data and import from the old State to the rebalanced State
      if(bondFieldIds->size() > 0){
        rebalancedState->allocateBondData(*bondFieldIds, rebalancedOwnedBondMap);      
        Epetra_Import importer(*rebalancedOwnedBondMap, *ownedBondMap);
        rebalancedState->getBondMultiVector()->Import(*state->getBondMultiVector(), importer, Insert);
      }

      // Set the State to the rebalanced State
      if(iState == 0)
        stateNONE = rebalancedState;
      else if(iState == 1)
        stateN = rebalancedState;
      else if(iState == 2)
        stateNP1 = rebalancedState;
    }
  }

  // Store the rebalanced maps
  ownedScalarPointMap = rebalancedOwnedScalarPointMap;
  overlapScalarPointMap = rebalancedOverlapScalarPointMap;
  ownedVectorPointMap = rebalancedOwnedVectorPointMap;
  overlapVectorPointMap = rebalancedOverlapVectorPointMap;
  ownedBondMap = rebalancedOwnedBondMap;
}

Teuchos::RCP<const Epetra_Comm> PeridigmNS::DataManager::getEpetraComm()
{
  Teuchos::RCP<const Epetra_Comm> comm;
  if(!ownedScalarPointMap.is_null())
    comm = Teuchos::rcpFromRef(ownedScalarPointMap->Comm());
  else if(!overlapScalarPointMap.is_null())
    comm = Teuchos::rcpFromRef(overlapScalarPointMap->Comm());
  else if(!ownedVectorPointMap.is_null())
    comm = Teuchos::rcpFromRef(ownedVectorPointMap->Comm());
  else if(!overlapVectorPointMap.is_null())
    comm = Teuchos::rcpFromRef(overlapVectorPointMap->Comm());
  else if(!ownedBondMap.is_null())
    comm = Teuchos::rcpFromRef(ownedBondMap->Comm());
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::NullReferenceError, 
                               "Error in PeridigmNS::DataManager::getEpetraComm(), no comm object available (forget setMaps()?).");
  return comm;
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

  // Check for data in State objects
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

  // Check for data in global data
  if(!hasData){
    if( fieldIdAndStepToGlobalData.find( pair<int, PeridigmField::Step>(fieldId, step) ) != fieldIdAndStepToGlobalData.end() )
      hasData = true;
  }

  return hasData;
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::DataManager::getData(int fieldId, PeridigmField::Step step)
{
  Teuchos::RCP<Epetra_Vector> data;

  if( fieldManager.isGlobalSpec(fieldId) ){
    data = fieldIdAndStepToGlobalData[ pair<int, PeridigmField::Step>(fieldId, step) ];
  }
  else{
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
  }

  if(data.is_null()){
    stringstream ss;
    ss << "**** Error, PeridigmNS::DataManager::getData(), fieldId and Step not found!\n";
    ss << "**** Spec: " << fieldManager.getFieldSpec(fieldId) << "\n";
    ss << "**** Step: " << step << "\n";
    TEUCHOS_TEST_FOR_EXCEPTION(data.is_null(), Teuchos::RangeError, ss.str());
  }

  return data;
}
