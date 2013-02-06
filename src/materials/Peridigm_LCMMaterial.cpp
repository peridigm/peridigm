/*! \file Peridigm_LCMMaterial.cpp */

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

#include "Peridigm_LCMMaterial.hpp"
#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>

#ifdef PERIDIGM_LCM

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
//#include "QCAD_MaterialDatabase.hpp"
#include "Phalanx.hpp"

#include "PHAL_AlbanyTraits.hpp"
#include "PHAL_SaveStateField.hpp"
#include "Albany_Utils.hpp"
#include "Albany_StateManager.hpp"
#include "Albany_TmplSTKMeshStruct.hpp"
#include "Albany_STKDiscretization.hpp"
#include "Albany_Layouts.hpp"

#include "LCM/evaluators/SetField.hpp"
#include "LCM/evaluators/Neohookean.hpp"
#include "LCM/evaluators/J2Stress.hpp"

#include "Tensor.h"

#endif

using namespace std;

PeridigmNS::LCMMaterial::LCMMaterial(const Teuchos::ParameterList& params)
  : Material(params)
{
#ifndef PERIDIGM_LCM
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  LCM material models are not enabled.  Recompile with -D USE_LCM:BOOL=ON.\n");
#endif

  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_shearModulus = params.get<double>("Shear Modulus");
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId           = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_damageFieldId           = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId     = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId       = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
}

PeridigmNS::LCMMaterial::~LCMMaterial()
{
}

void
PeridigmNS::LCMMaterial::initialize(const double dt,
                                      const int numOwnedPoints,
                                      const int* ownedIDs,
                                      const int* neighborhoodList,
                                      PeridigmNS::DataManager& dataManager) const
{
#ifdef PERIDIGM_LCM

  typedef PHX::MDField<PHAL::AlbanyTraits::Residual::ScalarT>::size_type size_type;
  typedef PHAL::AlbanyTraits::Residual Residual;
  typedef PHAL::AlbanyTraits::Residual::ScalarT ScalarT;
  typedef PHAL::AlbanyTraits Traits;

  // Set up the data layout
  const int worksetSize = 1;
  const int numQPts = 1;
  const int numDim = 3;
  const int numVertices = 1;
  const int numNodes = 1;
  const Teuchos::RCP<Albany::Layouts> dl =
    Teuchos::rcp(new Albany::Layouts(worksetSize,numVertices,numNodes,numQPts,numDim));

  // Instantiate the required evaluators with EvalT = PHAL::AlbanyTraits::Residual and Traits = PHAL::AlbanyTraits

  //---------------------------------------------------------------------------
  // Deformation gradient
  Teuchos::ArrayRCP<ScalarT> defgrad(9);

  defGrad[0] = 1.0;  defGrad[1] = 0.0;  defGrad[2] = 0.0;
  defGrad[0] = 0.0;  defGrad[1] = 1.0;  defGrad[2] = 0.0;
  defGrad[0] = 0.0;  defGrad[1] = 0.0;  defGrad[2] = 1.0;

  // SetField evaluator, which will be used to manually assign a value to the defgrad field
  Teuchos::ParameterList setDefGradParams("SetFieldDefGrad");
  setDefGradParams.set<string>("Evaluated Field Name", "F");
  setDefGradParams.set<Teuchos::RCP<PHX::DataLayout> >("Evaluated Field Data Layout", dl->qp_tensor);
  setDefGradParams.set< Teuchos::ArrayRCP<ScalarT> >("Field Values", defgrad);
  // Teuchos::RCP<LCM::SetField<Residual, Traits> > setFieldDefGrad =
  //   Teuchos::rcp(new LCM::SetField<Residual, Traits>(setDefGradP));


#endif
}

void
PeridigmNS::LCMMaterial::computeForce(const double dt,
                                        const int numOwnedPoints,
                                        const int* ownedIDs,
                                        const int* neighborhoodList,
                                        PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}
