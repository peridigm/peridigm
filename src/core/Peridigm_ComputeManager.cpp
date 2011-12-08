//*! \file Peridigm_ComputeManager.cpp */
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

#include <vector>
#include <string>
#include <iostream>

#include "Peridigm_ComputeManager.hpp"
#include "compute/compute_includes.hpp"
#include "Peridigm_DataManager.hpp"
#include "mesh_output/Field.h"

PeridigmNS::ComputeManager::ComputeManager( Teuchos::RCP<Teuchos::ParameterList>& params, PeridigmNS::Peridigm *peridigm_ ) {

  Teuchos::RCP<Compute> compute;

  // Hook up parent pointer
  peridigm = peridigm_;

  // No input to validate; no computes requested
  if (params == Teuchos::null) return;

  Teuchos::RCP<Teuchos::ParameterList> materialOutputFields; 
  // User-requested fields for output
  if (params->isSublist("Material Output Fields"))
    materialOutputFields = sublist(params, "Material Output Fields");
  else 
    // No input to validate; no computes requested
    return;

  Teuchos::ParameterList::ConstIterator i1;
  // Loop over the material types in the materialOutputFields parameterlist
  for (i1 = materialOutputFields->begin(); i1 != materialOutputFields->end(); ++i1) {
    const Teuchos::ParameterEntry& val1 = materialOutputFields->entry(i1);
    // const std::string& name1 = materialOutputFields->name(i1);
    // For each material type, loop over requested output fields and hook up pointers
    if (val1.isList()) { // each material type is a sublist
      const Teuchos::ParameterList& sublist = Teuchos::getValue<Teuchos::ParameterList>(val1);
      Teuchos::ParameterList::ConstIterator i2;
      for (i2=sublist.begin(); i2 != sublist.end(); ++i2) {
        const std::string& nm = sublist.name(i2);

        #define COMPUTE_CLASS
        #define ComputeClass(key, Class, Args) \
        if (nm == #key) { \
          compute = Teuchos::rcp( new PeridigmNS::Class(Args) ); \
          computeObjects.push_back( Teuchos::rcp_implicit_cast<Compute>(compute) ); \
        }
        #include "compute/compute_includes.hpp"
        #undef  COMPUTE_CLASS

      }
    }
  }

}

Teuchos::ParameterList PeridigmNS::ComputeManager::getValidParameterList() {
  Teuchos::ParameterList validParameterList("Output");
  return validParameterList;
}

std::vector<Field_NS::FieldSpec> PeridigmNS::ComputeManager::getFieldSpecs() {

  std::vector<Field_NS::FieldSpec> myFieldSpecs;

  // Loop over all compute objects, collect the field specs they compute
  for (unsigned int i=0; i < computeObjects.size(); i++) {
    Teuchos::RCP<const PeridigmNS::Compute> compute = computeObjects[i];
    std::vector<Field_NS::FieldSpec> computeFieldSpecs = compute->getFieldSpecs();
    myFieldSpecs.insert(myFieldSpecs.end(), computeFieldSpecs.begin(), computeFieldSpecs.end());
  }

  // remove duplicates
  std::sort(myFieldSpecs.begin(), myFieldSpecs.end());
  std::unique(myFieldSpecs.begin(), myFieldSpecs.end());

  return myFieldSpecs;
}

PeridigmNS::ComputeManager::~ComputeManager() {
}

void PeridigmNS::ComputeManager::compute(Teuchos::RCP<PeridigmNS::DataManager>& dataManager) {

  for (unsigned int i=0; i < computeObjects.size(); i++)
    computeObjects[i]->compute(dataManager);

}
