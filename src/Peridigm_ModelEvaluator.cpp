/*! \file Peridigm_ModelEvaluator.hpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

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

  constructEvaluators();
  fm->postRegistrationSetup(NULL);
}

PeridigmNS::ModelEvaluator::~ModelEvaluator(){
}

void 
PeridigmNS::ModelEvaluator::constructEvaluators()
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
  fm = Teuchos::rcp(new PHX::FieldManager<PHAL::PeridigmTraits>);

  // Register all Evaluators with the field manager
  PHX::registerEvaluators(evaluators, *fm);

  // List EvaluateForce as a required field for the field manager (other evaluators will be
  // called as perscribed by dependencies).
  PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> evaluate_force_tag("EvaluateForce", dummy);
  fm->requireField<PHAL::PeridigmTraits::Residual>(evaluate_force_tag);

  // Require the contact force evaluation
  if(analysisHasContact){
    PHX::Tag<PHAL::PeridigmTraits::Residual::ScalarT> contact_tag("Contact", dummy);
    fm->requireField<PHAL::PeridigmTraits::Residual>(contact_tag);
  }
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
  fm->evaluateFields<PHAL::PeridigmTraits::Residual>(*workset);
}
