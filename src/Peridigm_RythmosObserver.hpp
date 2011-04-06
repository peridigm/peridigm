/*! \file Peridigm_RythmosObserver.hpp */
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

#ifndef PERIDIGM_RYTHMOSOBSERVER_HPP
#define PERIDIGM_RYTHMOSOBSERVER_HPP

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

  /*!
   * \brief Class for strategy object that observes and reports on time integration by observing the stepper object
   */
  class RythmosObserver : public Rythmos::IntegrationObserverBase<Scalar> {
  public:
     //! Constructor
     RythmosObserver (Teuchos::RCP<EpetraExt::ModelEvaluator>, const Teuchos::RCP<Teuchos::ParameterList>& params);

     ~RythmosObserver () {};

     Teuchos::RCP<Rythmos::IntegrationObserverBase<Scalar> >  cloneIntegrationObserver() const {
       TEST_FOR_EXCEPT(true);
     };

     void resetIntegrationObserver( const Rythmos::TimeRange<Scalar> &integrationTimeDomain ) {};

     //! Called by stepper after time step completed.
     void observeCompletedTimeStep(
       const Rythmos::StepperBase<Scalar> &stepper,
       const Rythmos::StepControlInfo<Scalar> &stepCtrlInfo,
       const int timeStepIter
    );

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

#endif // PERIDIGM_RYTHMOSOBSERVER_HPP
