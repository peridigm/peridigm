/*! \file Peridigm_FieldSpec.cpp */

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

#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>
#include <algorithm>

using namespace std;

PeridigmNS::FieldManager& PeridigmNS::FieldManager::self() {
  static FieldManager fieldManager;
  return fieldManager;
}

int PeridigmNS::FieldManager::getFieldId(PeridigmField::Relation relation_,
                                         PeridigmField::Length length_,
                                         PeridigmField::Temporal temporal_,
                                         std::string label_)
{
  // Begin with a check for an exact match
  for(vector<FieldSpec>::iterator it = fieldSpecs.begin() ; it != fieldSpecs.end() ; ++it){

    // If everything matches, return the field id
    if(relation_ == it->relation && length_ == it->length && temporal_ == it->temporal && label_ == it->label)
      return it->id;
  }

  // Check for fields that do not match but are suspiciously similar
  for(vector<FieldSpec>::iterator it = fieldSpecs.begin() ; it != fieldSpecs.end() ; ++it){

    // Check for labels that match or differ only by case
    string lcLabel_(label_);
    transform(lcLabel_.begin(), lcLabel_.end(), lcLabel_.begin(), ::tolower);
    string lcLabel(it->label);
    transform(lcLabel.begin(), lcLabel.end(), lcLabel.begin(), ::tolower);
    if(lcLabel_ == lcLabel){
      stringstream ss;
      ss << "\n**** Error:  getFieldId() found suspiciously similar field specifications.\n";
      ss << "\n****         Requested spec:\n";
      ss << "\n****           relation: " << relation_;
      ss << "\n****           length:   " << length_;
      ss << "\n****           temporal: " << temporal_;
      ss << "\n****           label:    " << label_;
      ss << "\n****         Close match:\n";
      ss << "\n****           relation: " << it->relation;
      ss << "\n****           length:   " << it->length;
      ss << "\n****           temporal: " << it->temporal;
      ss << "\n****           label:    " << it->label;
      TEUCHOS_TEST_FOR_EXCEPT_MSG(lcLabel_ == lcLabel, ss.str());
    }
  }

  // If the spec does not exist, create a new one
  int id = fieldSpecs.size();
  fieldSpecs.push_back( FieldSpec(relation_, length_, temporal_, label_, id) );
  labelToIdMap[label_] = id;
  if(relation_ == PeridigmField::GLOBAL)
    specIsGlobal.push_back(true);
  else
    specIsGlobal.push_back(false);

  return id;
}

bool PeridigmNS::FieldManager::hasField(std::string label)
{
  if(labelToIdMap.find(label) == labelToIdMap.end())
    return false;
  return true;
}

int PeridigmNS::FieldManager::getFieldId(std::string label)
{
  map<string, int>::iterator it = labelToIdMap.find(label);
  if(it == labelToIdMap.end()){
    string msg = "\n**** Error:  getFieldId(), label not found:  " + label + "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(it == labelToIdMap.end(), msg);
  }
  return it->second;
}

PeridigmNS::FieldSpec PeridigmNS::FieldManager::getFieldSpec(int fieldId)
{
  unsigned int id = static_cast<unsigned int>(fieldId);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(id < 0 || id >= fieldSpecs.size(), "\n**** Error:  getFieldSpec(), ID not found.\n");
  return fieldSpecs[id];
}

PeridigmNS::FieldSpec PeridigmNS::FieldManager::getFieldSpec(string label)
{
  return getFieldSpec( getFieldId(label) );
}

std::ostream& operator<<(std::ostream& os, const PeridigmNS::PeridigmField::Relation& relation){
  if(relation == PeridigmNS::PeridigmField::UNDEFINED_RELATION)
    os << "UNDEFINED_RELATION";
  else if(relation == PeridigmNS::PeridigmField::NODE)
    os << "NODE";
  else if(relation == PeridigmNS::PeridigmField::ELEMENT)
    os << "ELEMENT";
  else if(relation == PeridigmNS::PeridigmField::BOND)
    os << "BOND";
  else if(relation == PeridigmNS::PeridigmField::GLOBAL)
    os << "GLOBAL";
  return os;
}

std::ostream& operator<<(std::ostream& os, const PeridigmNS::PeridigmField::Length& length){
  if(length == PeridigmNS::PeridigmField::UNDEFINED_LENGTH)
    os << "UNDEFINED_LENGTH";
  else if(length == PeridigmNS::PeridigmField::SCALAR)
    os << "SCALAR";
  else if(length == PeridigmNS::PeridigmField::VECTOR)
    os << "VECTOR";
  return os;
}

std::ostream& operator<<(std::ostream& os, const PeridigmNS::PeridigmField::Temporal& temporal){
  if(temporal == PeridigmNS::PeridigmField::UNDEFINED_TEMPORAL)
    os << "UNDEFINED_TEMPORAL";
  else if(temporal == PeridigmNS::PeridigmField::CONSTANT)
    os << "CONSTANT";
  else if(temporal == PeridigmNS::PeridigmField::TWO_STEP)
    os << "TWO_STEP";
  return os;
}

std::ostream& operator<<(std::ostream& os, const PeridigmNS::PeridigmField::Step& step){
  if(step == PeridigmNS::PeridigmField::UNDEFINED_STEP)
    os << "UNDEFINED_STEP";
  else if(step == PeridigmNS::PeridigmField::STEP_NONE)
    os << "STEP_NONE";
  else if(step == PeridigmNS::PeridigmField::STEP_N)
    os << "STEP_N";
  else if(step == PeridigmNS::PeridigmField::STEP_NP1)
    os << "STEP_NP1";
  return os;
}

std::ostream& operator<<(std::ostream& os, const PeridigmNS::FieldSpec& fieldSpec){
  os << fieldSpec.getLabel() << " (" << fieldSpec.getRelation() << ", " << fieldSpec.getLength() << ", " << fieldSpec.getTemporal() << ")";
  return os;
}
