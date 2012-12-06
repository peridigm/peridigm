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

using namespace std;

PeridigmNS::ComputeManager::ComputeManager( Teuchos::RCP<Teuchos::ParameterList> params, Teuchos::RCP<const Epetra_Comm> epetraComm ) {

  Teuchos::RCP<Compute> compute;

  // No input to validate; no computes requested
  if (params == Teuchos::null) return;

  // Create a list of the compute classes that will be built based on user-supplied compute class parameters.
  // If the name of one of these compute classes appears in both the compute class parameters ParameterList and
  // the output variables ParameterList, we want to only instantiate one object and we want to use the compute
  // class parameters (as opposed to an empty parameter list, which is what is done when creating compute objects
  // based on the output variables ParameterList).
  vector<string> computeClassParameterNames;
  if (params->isSublist("Compute Class Parameters")) {
    Teuchos::RCP<Teuchos::ParameterList> computeClassParameters = sublist(params, "Compute Class Parameters");
    for (Teuchos::ParameterList::ConstIterator it = computeClassParameters->begin(); it != computeClassParameters->end(); ++it) {
      computeClassParameterNames.push_back(it->first);
    }
  }

  vector< pair<string, Teuchos::RCP<Teuchos::ParameterList> > > computeClassesToBuild;

  // Create compute classes based on requested output variables
  if (params->isSublist("Output Variables")) {
    Teuchos::RCP<Teuchos::ParameterList> outputVariables = sublist(params, "Output Variables");
    for (Teuchos::ParameterList::ConstIterator it = outputVariables->begin(); it != outputVariables->end(); ++it) {
      const string& name = it->first;
      if( find(computeClassParameterNames.begin(), computeClassParameterNames.end(), name) == computeClassParameterNames.end() ){
        Teuchos::RCP<Teuchos::ParameterList> nullRcp;
        pair<string, Teuchos::RCP<Teuchos::ParameterList> > nameAndParamsPair(name, nullRcp);
        computeClassesToBuild.push_back(nameAndParamsPair);
      }
    }
  }

  // Create compute classes based on user-supplied compute class parameters
  if (params->isSublist("Compute Class Parameters")) {
    Teuchos::RCP<Teuchos::ParameterList> computeClassParameters = sublist(params, "Compute Class Parameters");
    for (Teuchos::ParameterList::ConstIterator it = computeClassParameters->begin(); it != computeClassParameters->end(); ++it) {
      string parameterListName = it->first;
      Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcpFromRef( computeClassParameters->sublist(parameterListName) );
      string name = params->get<string>("Compute Class");
      pair<string, Teuchos::RCP<Teuchos::ParameterList> > nameAndParamsPair(name, params);
      computeClassesToBuild.push_back(nameAndParamsPair);
    }
  }

  // Instantiate the compute classes
  vector< pair<string, Teuchos::RCP<Teuchos::ParameterList> > >::iterator it;
  for (it = computeClassesToBuild.begin() ; it != computeClassesToBuild.end() ; it++) {
    string name = it->first;
    Teuchos::RCP<Teuchos::ParameterList> params = it->second;
    #define COMPUTE_CLASS
      #define ComputeClass(key, Class) \
      if (name == #key) { \
        compute = Teuchos::rcp( new PeridigmNS::Class(params, epetraComm) ); \
        computeObjects.push_back( Teuchos::rcp_implicit_cast<Compute>(compute) ); \
      }
      #include "compute/compute_includes.hpp"
    #undef  COMPUTE_CLASS
  }
}

Teuchos::ParameterList PeridigmNS::ComputeManager::getValidParameterList() {
  Teuchos::ParameterList validParameterList("Output");
  return validParameterList;
}

vector<int> PeridigmNS::ComputeManager::FieldIds() const {

  vector<int> myFieldIds;

  // Loop over all compute objects, collect the field ids they compute
  for (unsigned int i=0; i < computeObjects.size(); i++) {
    Teuchos::RCP<const PeridigmNS::Compute> compute = computeObjects[i];
    vector<int> computeFieldIds = compute->FieldIds();
    myFieldIds.insert(myFieldIds.end(), computeFieldIds.begin(), computeFieldIds.end());
  }

  // remove duplicates
  sort(myFieldIds.begin(), myFieldIds.end());
  vector<int>::iterator newEnd = unique(myFieldIds.begin(), myFieldIds.end());
  myFieldIds.erase(newEnd, myFieldIds.end());

  return myFieldIds;
}

PeridigmNS::ComputeManager::~ComputeManager() {
}

void PeridigmNS::ComputeManager::initialize(Teuchos::RCP< vector<PeridigmNS::Block> > blocks) {

  // \todo Identify what the desired behavior is for compute classes and multiple blocks!
  //       Calling initialize on each block individually may not make sense.

  for(unsigned int i=0 ; i<computeObjects.size() ; ++i){
     computeObjects[i]->initialize(blocks);
  }

}

void PeridigmNS::ComputeManager::compute(Teuchos::RCP< vector<PeridigmNS::Block> > blocks) {

  // \todo Identify what the desired behavior is for compute classes and multiple blocks!

  for(unsigned int i=0 ; i<computeObjects.size() ; ++i){
     computeObjects[i]->compute(blocks);
  }

}
