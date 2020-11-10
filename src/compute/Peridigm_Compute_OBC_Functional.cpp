/*! \file Peridigm_Compute_OBC_Functional.cpp */

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

#include "Peridigm_Compute_OBC_Functional.hpp"
#include "Peridigm_Field.hpp"

PeridigmNS::Compute_OBC_Functional::Compute_OBC_Functional(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                           Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                           Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_), m_obcFunctionalFieldId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_obcFunctionalFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "OBC_Functional");
  m_fieldIds.push_back(m_obcFunctionalFieldId);
}

PeridigmNS::Compute_OBC_Functional::~Compute_OBC_Functional(){}

void PeridigmNS::Compute_OBC_Functional::initialize( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ){
  std::vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->getData(m_obcFunctionalFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  }
}

int PeridigmNS::Compute_OBC_Functional::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  // This class is a no-op, it's only purpose is to register the "OBC_Functional" field and to hook up
  // the plumbing so that this field can be output to the exodus file.

  // The filling of this field is handled by Albany in an optimization-based coupling (OBC) Albany-Peridigm problem.

  return 0;
}
