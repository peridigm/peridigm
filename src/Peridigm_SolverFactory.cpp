/*! \file Peridigm_SolverFactory.cpp */

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

#include "Peridigm_SolverFactory.hpp"
#include "ENAT_RythmosSolver.hpp"
#include "Peridigm_VerletSolver.hpp"
#include "Peridigm_RythmosObserver.hpp"
#include "Peridigm_VerletObserver.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

Peridigm::SolverFactory::SolverFactory(const std::string inputFile, const MPI_Comm& appComm) 
{
  using Teuchos::RCP;
  using Teuchos::rcp;

#ifdef HAVE_MPI
    Comm = rcp(new Epetra_MpiComm(appComm));
#else
    Comm = rcp(new Epetra_SerialComm);
#endif

    // Set application parameters to default values
    appParams = rcp(new Teuchos::ParameterList());
    setProblemParamDefaults(appParams.ptr());
    setSolverParamDefaults(appParams.ptr());

    // Update parameters with data from xml file
    Teuchos::updateParametersFromXmlFile(inputFile, appParams.get());

    // create the model evaluator specific to the problem being solved
    model = rcp(new Peridigm::ModelEvaluator(Comm, appParams));

}

Teuchos::RCP<EpetraExt::ModelEvaluator> Peridigm::SolverFactory::create()
{
  Teuchos::RCP<Teuchos::ParameterList> solverPL = sublist(appParams,"Solver",true);
  Teuchos::RCP<EpetraExt::ModelEvaluator> retval;
  if (solverPL->isSublist("Rythmos")) { 
     // create the Rythmos Observer and Solver
    Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> > observer = Teuchos::rcp(new Peridigm::RythmosObserver(model,appParams));
    retval = Teuchos::rcp(new ENAT::RythmosSolver(appParams, model, observer));
  }
  else if (solverPL->isSublist("Verlet")) {
     // create the Verlet Observer and Solver
    Teuchos::RCP<Peridigm::VerletObserver> observer = Teuchos::rcp(new Peridigm::VerletObserver(model,appParams));
    retval = Teuchos::rcp(new Peridigm::VerletSolver(appParams, model, observer));
  }
  else { // Solver not recognized
    TEST_FOR_EXCEPTION( !solverPL->isSublist("Rythmos") && !solverPL->isSublist("Verlet"),
                        std::invalid_argument,
                        "Peridigm::SolverFactory: \"Solver\" list must include arguments for either \"Rythmos\" or \"Verlet\".");
  }
  return retval;
}

void Peridigm::SolverFactory::setProblemParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> appParams_)
{
  Teuchos::ParameterList& problemParams = appParams_->sublist("Problem");	

  // general settings
  problemParams.set("Verbose", false);
}

void Peridigm::SolverFactory::setSolverParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> appParams_)
{
  Teuchos::ParameterList& solverParams = appParams_->sublist("Solver");	

  // general settings
  solverParams.set("Verbose", false);

}
