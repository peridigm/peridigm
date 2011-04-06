/*! \file Peridigm_SolverFactory.cpp */

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
