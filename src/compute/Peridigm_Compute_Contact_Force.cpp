/*! \file Peridigm_Compute_Contact_Force.cpp */

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

#include "Peridigm_Compute_Contact_Force.hpp"
#include "../Peridigm.hpp"

//! Standard constructor.
PeridigmNS::Compute_Contact_Force::Compute_Contact_Force(PeridigmNS::Peridigm *peridigm_ ){peridigm = peridigm_;}

//! Destructor.
PeridigmNS::Compute_Contact_Force::~Compute_Contact_Force(){}

//! Returns the fieldspecs computed by this class
std::vector<Field_NS::FieldSpec> PeridigmNS::Compute_Contact_Force::getFieldSpecs() const {
  std::vector<Field_NS::FieldSpec> myFieldSpecs;
  myFieldSpecs.push_back(Field_NS::CONTACT_FORCE3D);

  return myFieldSpecs;
}

//! Fill the contact force vector
int PeridigmNS::Compute_Contact_Force::compute(Teuchos::RCP<PeridigmNS::DataManager> dataManager) const {

  int retval;

  Teuchos::RCP<Epetra_Vector> contact_force, contact_force_density, volume;
  contact_force_density = dataManager->getData(Field_NS::CONTACT_FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1);
  contact_force         = dataManager->getData(Field_NS::CONTACT_FORCE3D, Field_NS::FieldSpec::STEP_NP1);
  volume                = dataManager->getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE);

  // Sanity check
  if ( (contact_force_density->Map().NumMyElements() != volume->Map().NumMyElements()) ||  (contact_force->Map().NumMyElements() != volume->Map().NumMyElements()) ) {
    retval = 1;
    return(retval);
  }

  *contact_force = *contact_force_density;

  double *volume_values = volume->Values();
  double *contact_force_values  = contact_force->Values();

  // volume is a scalar and force a vector, so maps are different; must do multiplication on per-element basis
  int numElements = volume->Map().NumMyElements();
  double vol;
  for (int i=0;i<numElements;i++) {
    vol = volume_values[i]; 
    contact_force_values[3*i] = vol*contact_force_values[3*i];
    contact_force_values[3*i+1] = vol*contact_force_values[3*i+1];
    contact_force_values[3*i+2] = vol*contact_force_values[3*i+2];
  }

  return(0);

}

