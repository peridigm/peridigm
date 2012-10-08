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

#include "Peridigm_FieldSpec.hpp"
#include <Teuchos_Assert.hpp>
#include <algorithm>

using namespace std;

const PeridigmNS::FieldSpecManager& PeridigmNS::FieldSpecManager::self() {
  static const FieldSpecManager fieldSpecManager;
  return fieldSpecManager;
}

int PeridigmNS::FieldSpecManager::getFieldSpecID(PeridigmField::Relation relation_,
                                                 PeridigmField::Length length_,
                                                 PeridigmField::Temporal temporal_,
                                                 std::string label_)
{
  // Begin with a check for an exact match
  for(vector<FieldSpec>::iterator it = fieldSpecs.begin() ; it != fieldSpecs.end() ; ++it){

    // If everything matches, return the spec id
    if(relation_ == it->relation && length_ == it->length && temporal_ == it->temporal && label_ == it->label)
      return it->id;
  }

  // Check for specs that do not match but are suspiciously similar
  for(vector<FieldSpec>::iterator it = fieldSpecs.begin() ; it != fieldSpecs.end() ; ++it){

    // Check for labels that match or differ only by case
    string lcLabel_(label_);
    transform(lcLabel_.begin(), lcLabel_.end(), lcLabel_.begin(), ::tolower);
    string lcLabel(it->label);
    transform(lcLabel.begin(), lcLabel.end(), lcLabel.begin(), ::tolower);
    if(lcLabel_ == lcLabel){
      stringstream ss;
      ss << "\n**** Error:  getFieldSpecId() found suspiciously similar field specs.\n";
      ss << "\n****         Requested spec:\n";
      ss << "\n****           relation: " << relation_ + "\n";
      ss << "\n****           length:   " << length_ + "\n";
      ss << "\n****           temporal: " << temporal_ + "\n";
      ss << "\n****           label:    " << label_ + "\n";
      ss << "\n****         Close match:\n";
      ss << "\n****           relation: " << it->relation + "\n";
      ss << "\n****           length:   " << it->length + "\n";
      ss << "\n****           temporal: " << it->temporal + "\n";
      ss << "\n****           label:    " << it->label + "\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(lcLabel_ == lcLabel, ss.str());
    }
  }

  // If the spec does not exist, create a new one
  int id = fieldSpecs.size();
  fieldSpecs.push_back( FieldSpec(relation_, length_, temporal_, label_, id) );
  labelToIdMap[label_] = id;
  
  return id;
}

bool PeridigmNS::FieldSpecManager::hasFieldSpec(std::string label)
{
  if(labelToIdMap.find(label) == labelToIdMap.end())
    return false;
  return true;
}

int PeridigmNS::FieldSpecManager::getFieldSpecID(std::string label)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(labelToIdMap.find(label) == labelToIdMap.end(), "\n**** Error:  getFieldSpecID(), label not found.\n");
  return labelToIdMap[label];
}
