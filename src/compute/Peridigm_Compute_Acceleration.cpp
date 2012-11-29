/*! \file Peridigm_Compute_Acceleration.cpp */

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

#include "Peridigm_Compute_Acceleration.hpp"
#include "Peridigm_Field.hpp"

//! Standard constructor.
PeridigmNS::Compute_Acceleration::Compute_Acceleration(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                       Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(params, epetraComm_), m_forceDensityFieldId(-1), m_accelerationFieldId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_forceDensityFieldId = fieldManager.getFieldId("Force_Density");
  m_accelerationFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration");
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_accelerationFieldId);
}

//! Destructor.
PeridigmNS::Compute_Acceleration::~Compute_Acceleration(){}

//! Fill the acceleration vector
int PeridigmNS::Compute_Acceleration::compute( Teuchos::RCP< std::vector<Block> > blocks ) const {
  int retval(0);

  Teuchos::RCP<Epetra_Vector> force, acceleration;
  std::vector<Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    force        = blockIt->getData(m_forceDensityFieldId, PeridigmField::STEP_NP1);
    acceleration = blockIt->getData(m_accelerationFieldId, PeridigmField::STEP_NP1);
    *acceleration = *force;
    double density = blockIt->getMaterialModel()->Density();
    // Report if any calls to Scale() failed.
    retval = retval || acceleration->Scale(1.0/density);
  }

  return retval;

}

