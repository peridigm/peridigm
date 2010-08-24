/*! \file Peridigm_Factory.cpp */

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

#include <Teuchos_RCP.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Peridigm_Factory.hpp"
#include "Peridigm.hpp"

Peridigm::PeridigmFactory::PeridigmFactory()
{
}

Teuchos::RCP<Peridigm::Peridigm> Peridigm::PeridigmFactory::create(const std::string inputFile, const MPI_Comm& peridigmComm)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::RCP<Epetra_Comm> comm;
  #ifdef HAVE_MPI
    comm = rcp(new Epetra_MpiComm(peridigmComm));
  #else
    comm = rcp(new Epetra_SerialComm);
  #endif

  // Set application parameters to default values
  Teuchos::RCP<Teuchos::ParameterList> peridigmParams = rcp(new Teuchos::ParameterList());
  setProblemParamDefaults(peridigmParams.ptr());
  setSolverParamDefaults(peridigmParams.ptr());

  // Update parameters with data from xml file
  Teuchos::updateParametersFromXmlFile(inputFile, peridigmParams.get());

  // Create new Peridigm object
  return rcp(new Peridigm::Peridigm(comm, peridigmParams));

}

void Peridigm::PeridigmFactory::setProblemParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_)
{
  Teuchos::ParameterList& problemParams = peridigmParams_->sublist("Problem");

  // general settings
  problemParams.set("Verbose", false);
}

void Peridigm::PeridigmFactory::setSolverParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_)
{
  Teuchos::ParameterList& solverParams = peridigmParams_->sublist("Solver");

  // general settings
  solverParams.set("Verbose", false);
}
