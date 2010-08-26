/*! \file Peridigm_VerletObserver.hpp */
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

#ifndef PERIDIGM_VERLETOBSERVER_HPP
#define PERIDIGM_VERLETOBSERVER_HPP

#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>
#include <Rythmos_StepperBase.hpp>
#include <Rythmos_IntegrationObserverBase.hpp>
#include <Rythmos_TimeRange.hpp>
#include <EpetraExt_ModelEvaluator.h>
#include "Peridigm_ModelEvaluator.hpp"
#include "Peridigm_OutputManager.hpp"

namespace PeridigmNS {

  typedef double Scalar;

  class VerletObserver {
  public:
     //! Constructor
     VerletObserver (Teuchos::RCP<EpetraExt::ModelEvaluator>, const Teuchos::RCP<Teuchos::ParameterList>& params);

     ~VerletObserver () {};

     //! Called by stepper after time step completed.
    void observeCompletedTimeStep(Teuchos::RCP<const Epetra_Vector>, double);

  private:

   Teuchos::RCP<Peridigm::ModelEvaluator> model;
   Teuchos::RCP<Epetra_Map> solutionMap;
   Teuchos::RCP<Peridigm::OutputManager> outputManager;

   // Parameterlist to hold output parameters from input deck
   Teuchos::RCP<Teuchos::ParameterList> outputParams;
   
   // Parameterlist to hold description of force state data (from material model)
   Teuchos::RCP<Teuchos::ParameterList> forceStateDesc;
   
   // Flag to indicate if Observer is "turned on". Will update this later to be set from input script
   bool active;

};

}

#endif // PERIDIGM_VERLETOBSERVER_HPP
