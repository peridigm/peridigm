/*! \file Peridigm_VerletObserver.cpp */
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

#include "Peridigm_VerletObserver.hpp"
#include <EpetraExt_ModelEvaluator.h>
#include <Epetra_Vector.h>
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include "Peridigm_OutputManager_VTK_XML.hpp"

PeridigmNS::VerletObserver::VerletObserver(Teuchos::RCP<EpetraExt::ModelEvaluator> model_,
                                         const Teuchos::RCP<Teuchos::ParameterList>& params) {

  active = false;

  model = Teuchos::rcp_dynamic_cast<PeridigmNS::ModelEvaluator>(model_);
  TEUCHOS_TEST_FOR_EXCEPT_MSG( model.get() == NULL, "PeridigmNS::RythmosObserver: PeridigmNS::ModelEvaluator not passed in.");

  if (params->isSublist("Output")) {
    active = true;
    outputParams =  Teuchos::rcp(&(params->sublist("Output")),false);
    outputParams->set("NumProc", (int)(model->get_x_map()->Comm()).NumProc());
    outputParams->set("MyPID", (int)(model->get_x_map()->Comm()).MyPID());
  }

  if (active) {
    // Make the default format "VTK_XML"
    string outputFormat = outputParams->get("Output File Type", "VTK_XML");
    TEUCHOS_TEST_FOR_EXCEPTION( outputFormat != "VTK_XML",
                        std::invalid_argument,
                        "PeridigmNS::RythmosObserver: \"Output File Type\" must be either \"VTK_XML\".");
    if (outputFormat == "VTK_XML")
       outputManager = Teuchos::rcp(new PeridigmNS::OutputManager_VTK_XML( outputParams ));
    else
      TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,"PeridigmNS::RythmosObserver: \"Output File Type\" must be \"VTK_XML\".");

    // Query material models for their force state data descriptions
    forceStateDesc = Teuchos::rcp( new Teuchos::ParameterList() );
    std::vector<Teuchos::RCP<const PeridigmNS::Material> > materials = *(model->getMaterials());
    for(unsigned int i=0; i<materials.size(); ++i){
      Teuchos::ParameterList& subList = forceStateDesc->sublist(materials[i]->Name());
      for(int j=0;j<materials[i]->NumScalarConstitutiveVariables(); ++j){
        subList.set( materials[i]->ScalarConstitutiveVariableName(j), j );
      }
    }
    // Initialize current time in this parameterlist
    forceStateDesc->set("Time", 0.0);
    // Set RCP to neighborlist
    forceStateDesc->set("Bond Family",model->getNeighborhoodData());
    // Ask OutputManager to write initial conditions to disk
    outputManager->write(model->get_x_init(),model->getScalarConstitutiveDataOverlap(),model->getNeighborhoodData(),forceStateDesc);
  }

  //  verbose = problemParams->get("Verbose", false);

}

void PeridigmNS::VerletObserver::observeCompletedTimeStep(Teuchos::RCP<const Epetra_Vector> currentSolution, double time) {

  // We have completed a time step, tell the Peridigm_ModelEvaluator to update its state information
  model->updateState();

  // Callback allowing the ModelEvaluator to update the contact configuration, if necessary
  model->updateContact(currentSolution);

  // Only report status if Observer is active
  if (!active) return;

  //cout << "PERIDIGM OBSERVER CALLED step=" <<  timeStepIter  << ",  time=" << stepper.getStepStatus().time << endl;

  // Set current time in this parameterlist
  forceStateDesc->set("Time", time);
  outputManager->write(currentSolution,model->getScalarConstitutiveDataOverlap(),model->getNeighborhoodData(),forceStateDesc);

}
