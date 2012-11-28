/*! \file Peridigm_State.cpp */

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

#include "Peridigm_State.hpp"
#include "Peridigm_Field.hpp"
#include <Epetra_Import.h>
#include <Teuchos_Assert.hpp>
#include <sstream>

using namespace std;

void PeridigmNS::State::allocateScalarPointData(vector<int> fieldIds, Teuchos::RCP<const Epetra_BlockMap> map)
{
  std::sort(fieldIds.begin(), fieldIds.end());
  scalarPointData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldIds.size()));
  for(unsigned int i=0 ; i<fieldIds.size() ; ++i)
    fieldIdToDataMap[fieldIds[i]] = Teuchos::rcp((*scalarPointData)(i), false);
}

void PeridigmNS::State::allocateVectorPointData(vector<int> fieldIds, Teuchos::RCP<const Epetra_BlockMap> map)
{
  std::sort(fieldIds.begin(), fieldIds.end());
  vectorPointData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldIds.size()));
  for(unsigned int i=0 ; i<fieldIds.size() ; ++i)
    fieldIdToDataMap[fieldIds[i]] = Teuchos::rcp((*vectorPointData)(i), false);
}

void PeridigmNS::State::allocateScalarBondData(vector<int> fieldIds, Teuchos::RCP<const Epetra_BlockMap> map)
{
  std::sort(fieldIds.begin(), fieldIds.end());
  scalarBondData = Teuchos::rcp(new Epetra_MultiVector(*map, fieldIds.size()));
  for(unsigned int i=0 ; i<fieldIds.size() ; ++i)
    fieldIdToDataMap[fieldIds[i]] = Teuchos::rcp((*scalarBondData)(i), false);
}

vector<int> PeridigmNS::State::getFieldIds(PeridigmField::Relation relation,
										   PeridigmField::Length length)
{
  vector<int> fieldIds;
  FieldManager& fieldManager = FieldManager::self();
  std::map< int, Teuchos::RCP<Epetra_Vector> >::const_iterator it;
  for(it = fieldIdToDataMap.begin() ; it != fieldIdToDataMap.end() ; ++it){
    PeridigmNS::FieldSpec spec = fieldManager.getFieldSpec(it->first);
    if(spec.getLength() == length && spec.getRelation() == relation)
      fieldIds.push_back(it->first);
  }
  return fieldIds;
}

bool PeridigmNS::State::hasData(int fieldId)
{
  std::map< int, Teuchos::RCP<Epetra_Vector> >::iterator lb = fieldIdToDataMap.lower_bound(fieldId);
  bool keyExists = ( lb != fieldIdToDataMap.end() && !(fieldIdToDataMap.key_comp()(fieldId, lb->first)) );
  return keyExists;
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::State::getData(int fieldId)
{
  // search for the data
  std::map< int, Teuchos::RCP<Epetra_Vector> >::iterator lb = fieldIdToDataMap.lower_bound(fieldId);
  // if the key does not exist, throw an exception
  bool keyExists = ( lb != fieldIdToDataMap.end() && !(fieldIdToDataMap.key_comp()(fieldId, lb->first)) );
  if(!keyExists){
    std::stringstream ss;
    ss << fieldId;
    TEUCHOS_TEST_FOR_EXCEPTION(!keyExists, Teuchos::RangeError, 
                       "****Error in PeridigmNS::State::getData(), field id does not exist: " + ss.str() + "\n");
  }
  return lb->second;
}

void PeridigmNS::State::copyLocallyOwnedDataFromState(Teuchos::RCP<PeridigmNS::State> source)
{
  // Make sure the source isn't a null ref-count pointer
  TEUCHOS_TEST_FOR_EXCEPTION(source.is_null(), Teuchos::NullReferenceError,
                     "PeridigmNS::State::copyLocallyOwnedDataFromState() called with null ref-count pointer.\n");

  if(!scalarPointData.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(source->getScalarPointMultiVector().is_null(), Teuchos::NullReferenceError,
                       "PeridigmNS::State::copyLocallyOwnedDataFromState() called with incompatible State.\n");
    copyLocallyOwnedMultiVectorData( *(source->getScalarPointMultiVector()), *scalarPointData );
  }
  if(!vectorPointData.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(source->getVectorPointMultiVector().is_null(), Teuchos::NullReferenceError,
                       "PeridigmNS::State::copyLocallyOwnedDataFromState() called with incompatible State.\n");
    copyLocallyOwnedMultiVectorData( *(source->getVectorPointMultiVector()), *vectorPointData );
  }
  if(!scalarBondData.is_null()){
    TEUCHOS_TEST_FOR_EXCEPTION(source->getScalarBondMultiVector().is_null(), Teuchos::NullReferenceError,
                       "PeridigmNS::State::copyLocallyOwnedDataFromState() called with incompatible State.\n");
    copyLocallyOwnedMultiVectorData( *(source->getScalarBondMultiVector()), *scalarBondData );
  }
}

void PeridigmNS::State::copyLocallyOwnedMultiVectorData(Epetra_MultiVector& source, Epetra_MultiVector& target)
{
  TEUCHOS_TEST_FOR_EXCEPTION(source.NumVectors() != target.NumVectors(), std::runtime_error,
                     "PeridigmNS::State::copyLocallyOwnedMultiVectorData() called with incompatible MultiVectors.\n");
  int numVectors = target.NumVectors();
  const Epetra_BlockMap& sourceMap = source.Map();
  const Epetra_BlockMap& targetMap = target.Map();
  for(int iVec=0 ; iVec<numVectors ; ++iVec){
    Epetra_Vector& sourceVector = *source(iVec);
    Epetra_Vector& targetVector = *target(iVec);
    for(int targetLID=0 ; targetLID<targetMap.NumMyElements() ; ++targetLID){
      int GID = targetMap.GID(targetLID);
      int sourceLID = sourceMap.LID(GID);
      TEUCHOS_TEST_FOR_EXCEPTION(sourceLID == -1, std::range_error,
                         "PeridigmNS::State::copyLocallyOwnedMultiVectorData() called with incompatible MultiVectors.\n");
      TEUCHOS_TEST_FOR_EXCEPTION(sourceMap.ElementSize(sourceLID) != targetMap.ElementSize(targetLID), std::range_error,
                         "PeridigmNS::State::copyLocallyOwnedMultiVectorData() called with incompatible MultiVectors.\n");
      int elementSize = targetMap.ElementSize(targetLID);
      int sourceFirstPointInElement = sourceMap.FirstPointInElement(sourceLID);
      int targetFirstPointInElement = targetMap.FirstPointInElement(targetLID);
      for(int i=0 ; i<elementSize ; ++i){
        targetVector[targetFirstPointInElement+i] = sourceVector[sourceFirstPointInElement+i];
      }
    }
  }
}
