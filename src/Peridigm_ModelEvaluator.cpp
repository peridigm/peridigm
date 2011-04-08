/*! \file Peridigm_ModelEvaluator.hpp */

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

#include <Teuchos_RCPDecl.hpp>
#include <Teuchos_VerboseObject.hpp>
#include "PHAL_Dimension.hpp"
#include "PHAL_FactoryTraits.hpp"
#include "Peridigm_ModelEvaluator.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "PdGridData.h"
#include "PdQuickGrid.h"
#include "PdZoltan.h"
#include <sstream>

PeridigmNS::ModelEvaluator::ModelEvaluator(const Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > > materialModels_,
                                           const Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::ContactModel> > > contactModels_,
                                           const Teuchos::RCP<const Epetra_Comm>& comm)
  : materialModels(materialModels_),
    contactModels(contactModels_),
    analysisHasContact(false),
	numPID(comm->NumProc()),
	myPID(comm->MyPID()),
    verbose(false)
{
  if(contactModels->size() > 0)
    analysisHasContact = true;

  constructForceEvaluators();
  forceFieldManager->postRegistrationSetup(NULL);

  // \todo Call only if needed (implicit time integration)
  constructJacobianEvaluators();
  jacobianFieldManager->postRegistrationSetup(NULL);
}

PeridigmNS::ModelEvaluator::~ModelEvaluator(){
}

void 
PeridigmNS::ModelEvaluator::constructForceEvaluators()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::vector;
  using std::map;
  using PHX::DataLayout;
  using PHX::MDALayout;

  // Create a list of evaluators to build

  // Each field evaluator will be associated with one or more data
  // types and data layouts.

  map<string, RCP<ParameterList> > evaluatorsToBuild;

  // Create a dummy data layout that can be associated with an evaluator if needed
  //! \todo Determine how the field manager's data layout scheme works and set it up properly.
  RCP<DataLayout> dummy = rcp(new MDALayout<Dummy>(0));

  { // Update force state
  RCP<ParameterList> p = rcp(new ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_update_force_state;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", verbose);

  // Associate the ParameterLibrary with the evaluator
  // this grants outside access (e.g. Dakota) to the parameters in the evaluator
//   p->set<RCP<ParamLib> >("Parameter Library", paramLib);

  // set the coefficient to the value in the parameter list
  // the parameter list could have been read from file or set to a default value
//   p->set<double>("Coefficient", params->get("Coefficient",1.0)); 
//   p->set<Teuchos::ParameterList*>("Parameter List", &paramList);

  p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

  evaluatorsToBuild["UpdateForceState"] = p;
  }

  { // Evaluate force
  RCP<ParameterList> p = rcp(new ParameterList);
  int type = FactoryTraits<PHAL::PeridigmTraits>::id_evaluate_force;
  p->set<int>("Type", type); 
  p->set<bool>("Verbose", verbose);

  p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

  evaluatorsToBuild["EvaluateForce"] = p;
  }

  // Contact
  if(analysisHasContact){
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = FactoryTraits<PHAL::PeridigmTraits>::id_contact;
    p->set<int>("Type", type); 
    p->set<bool>("Verbose", verbose);

    p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

    evaluatorsToBuild["Contact"] = p;
  }

  // Build Field Evaluators for each evaluation type
  PHX::EvaluatorFactory<PHAL::PeridigmTraits, FactoryTraits<PHAL::PeridigmTraits> > factory;
  RCP< vector< RCP<PHX::Evaluator_TemplateManager<PHAL::PeridigmTraits> > > > evaluators;
  evaluators = factory.buildEvaluators(evaluatorsToBuild);

  // Create a FieldManager
  forceFieldManager = Teuchos::rcp(new PHX::FieldManager<PHAL::PeridigmTraits>);

  // Register all Evaluators with the field manager
  PHX::registerEvaluators(evaluators, *forceFieldManager);

  // List EvaluateForce as a required field for the field manager (other evaluators will be
  // called as perscribed by dependencies).
  PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> evaluate_force_tag("EvaluateForce", dummy);
  forceFieldManager->requireField<PHAL::PeridigmTraits::Residual>(evaluate_force_tag);

  // Require the contact force evaluation
  if(analysisHasContact){
    PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> contact_tag("Contact", dummy);
    forceFieldManager->requireField<PHAL::PeridigmTraits::Residual>(contact_tag);
  }
}

void 
PeridigmNS::ModelEvaluator::constructJacobianEvaluators()
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;
  using std::vector;
  using std::map;
  using PHX::DataLayout;
  using PHX::MDALayout;

  // Create a list of evaluators to build

  // Each field evaluator will be associated with one or more data
  // types and data layouts.
  map<string, RCP<ParameterList> > evaluatorsToBuild;

  // Create a dummy data layout that can be associated with an evaluator if needed
  //! \todo Determine how the field manager's data layout scheme works and set it up properly.
  RCP<DataLayout> dummy = rcp(new MDALayout<Dummy>(0));

  // Jacobian
  {
    RCP<ParameterList> p = rcp(new ParameterList);
    int type = FactoryTraits<PHAL::PeridigmTraits>::id_evaluate_jacobian;
    p->set<int>("Type", type); 
    p->set<bool>("Verbose", verbose);

    p->set< RCP<DataLayout> >("Dummy Data Layout", dummy);

    evaluatorsToBuild["EvaluateJacobian"] = p;
  }

  // Build Field Evaluators for each evaluation type
  PHX::EvaluatorFactory<PHAL::PeridigmTraits, FactoryTraits<PHAL::PeridigmTraits> > factory;
  RCP< vector< RCP<PHX::Evaluator_TemplateManager<PHAL::PeridigmTraits> > > > evaluators;
  evaluators = factory.buildEvaluators(evaluatorsToBuild);

  // Create a FieldManager
  jacobianFieldManager = Teuchos::rcp(new PHX::FieldManager<PHAL::PeridigmTraits>);

  // Register all Evaluators with the field manager
  PHX::registerEvaluators(evaluators, *jacobianFieldManager);

  // List EvaluateForce as a required field for the field manager (other evaluators will be
  // called as perscribed by dependencies).
  PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> evaluate_jacobian_tag("EvaluateJacobian", dummy);
  jacobianFieldManager->requireField<PHAL::PeridigmTraits::Residual>(evaluate_jacobian_tag);
}

Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::Material> > >
PeridigmNS::ModelEvaluator::getMaterialModels() const
{
  return materialModels;
}

Teuchos::RCP<std::vector<Teuchos::RCP<const PeridigmNS::ContactModel> > >
PeridigmNS::ModelEvaluator::getContactModels() const
{
  return contactModels;
}

void 
PeridigmNS::ModelEvaluator::evalModel(Teuchos::RCP<PHAL::Workset> workset) const
{
  // call field manager with workset
  forceFieldManager->evaluateFields<PHAL::PeridigmTraits::Residual>(*workset);
}

void 
PeridigmNS::ModelEvaluator::evalJacobian(Teuchos::RCP<PHAL::Workset> workset) const
{
  // call field manager with workset
  // \todo We should be using the template type to control force/jacobian calls.
  jacobianFieldManager->evaluateFields<PHAL::PeridigmTraits::Residual>(*workset);
}
