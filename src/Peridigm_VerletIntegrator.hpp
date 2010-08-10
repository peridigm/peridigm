/*! \file Peridigm_VerletIntegrator.hpp */
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

#ifndef PERIDIGM_VERLETINTEGRATOR_HPP
#define PERIDIGM_VERLETINTEGRATOR_HPP

#include <Peridigm_VerletObserver.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace Peridigm {

  class VerletIntegrator {
  public:
     VerletIntegrator(Teuchos::RCP<EpetraExt::ModelEvaluator>, Teuchos::RCP<Teuchos::ParameterList>);

     ~VerletIntegrator() {};

     int setIntegrationObserver(Teuchos::RCP<Peridigm::VerletObserver>);

     int integrate(double);

  private:

     Teuchos::RCP<Peridigm::VerletObserver> observer;
     Teuchos::RCP<Teuchos::ParameterList> solverParams;
     Teuchos::RCP<EpetraExt::ModelEvaluator> model;
     double t_initial;
     double t_current;
     double dt;

};

}

#endif // PERIDIGM_VERLETINTEGRATOR_HPP
