/*! \file Peridigm_ElasticCorrespondenceMaterial2.cpp */

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

#include "Peridigm_ElasticCorrespondenceMaterial2.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticCorrespondenceMaterial2::ElasticCorrespondenceMaterial2(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial2(params),
    m_unrotatedRateOfDeformationXXFieldId(-1), m_unrotatedRateOfDeformationXYFieldId(-1), m_unrotatedRateOfDeformationXZFieldId(-1), 
    m_unrotatedRateOfDeformationYXFieldId(-1), m_unrotatedRateOfDeformationYYFieldId(-1), m_unrotatedRateOfDeformationYZFieldId(-1), 
    m_unrotatedRateOfDeformationZXFieldId(-1), m_unrotatedRateOfDeformationZYFieldId(-1), m_unrotatedRateOfDeformationZZFieldId(-1),
    m_stressXXFieldId(-1), m_stressXYFieldId(-1), m_stressXZFieldId(-1), 
    m_stressYXFieldId(-1), m_stressYYFieldId(-1), m_stressYZFieldId(-1), 
    m_stressZXFieldId(-1), m_stressZYFieldId(-1), m_stressZZFieldId(-1)
{

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_unrotatedRateOfDeformationXXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationXX");
  m_unrotatedRateOfDeformationXYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationXY");
  m_unrotatedRateOfDeformationXZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationXZ");
  m_unrotatedRateOfDeformationYXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationYX");
  m_unrotatedRateOfDeformationYYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationYY");
  m_unrotatedRateOfDeformationYZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationYZ");
  m_unrotatedRateOfDeformationZXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationZX");
  m_unrotatedRateOfDeformationZYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationZY");
  m_unrotatedRateOfDeformationZZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationZZ");
  m_stressXXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressXX");
  m_stressXYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressXY");
  m_stressXZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressXZ");
  m_stressYXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressYX");
  m_stressYYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressYY");
  m_stressYZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressYZ");
  m_stressZXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressZX");
  m_stressZYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressZY");
  m_stressZZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_StressZZ");

  m_fieldIds.push_back(m_unrotatedRateOfDeformationXXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZZFieldId);
  m_fieldIds.push_back(m_stressXXFieldId);
  m_fieldIds.push_back(m_stressXYFieldId);
  m_fieldIds.push_back(m_stressXZFieldId);
  m_fieldIds.push_back(m_stressYXFieldId);
  m_fieldIds.push_back(m_stressYYFieldId);
  m_fieldIds.push_back(m_stressYZFieldId);
  m_fieldIds.push_back(m_stressZXFieldId);
  m_fieldIds.push_back(m_stressZYFieldId);
  m_fieldIds.push_back(m_stressZZFieldId);
}

PeridigmNS::ElasticCorrespondenceMaterial2::~ElasticCorrespondenceMaterial2()
{
}

void
PeridigmNS::ElasticCorrespondenceMaterial2::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               PeridigmNS::DataManager& dataManager) const
{
  double *cauchyStressNXX, *cauchyStressNXY, *cauchyStressNXZ;
  double *cauchyStressNYX, *cauchyStressNYY, *cauchyStressNYZ;
  double *cauchyStressNZX, *cauchyStressNZY, *cauchyStressNZZ;
  dataManager.getData(m_stressXXFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNXX);
  dataManager.getData(m_stressXYFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNXY);
  dataManager.getData(m_stressXZFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNXZ);
  dataManager.getData(m_stressYXFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNYX);
  dataManager.getData(m_stressYYFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNYY);
  dataManager.getData(m_stressYZFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNYZ);
  dataManager.getData(m_stressZXFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNZX);
  dataManager.getData(m_stressZYFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNZY);
  dataManager.getData(m_stressZZFieldId, PeridigmField::STEP_N)->ExtractView(&cauchyStressNZZ);

  double *cauchyStressNP1XX, *cauchyStressNP1XY, *cauchyStressNP1XZ;
  double *cauchyStressNP1YX, *cauchyStressNP1YY, *cauchyStressNP1YZ;
  double *cauchyStressNP1ZX, *cauchyStressNP1ZY, *cauchyStressNP1ZZ;
  dataManager.getData(m_stressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1XX);
  dataManager.getData(m_stressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1XY);
  dataManager.getData(m_stressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1XZ);
  dataManager.getData(m_stressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1YX);
  dataManager.getData(m_stressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1YY);
  dataManager.getData(m_stressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1YZ);
  dataManager.getData(m_stressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1ZX);
  dataManager.getData(m_stressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1ZY);
  dataManager.getData(m_stressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1ZZ);

  double *unrotatedRateOfDeformationXX, *unrotatedRateOfDeformationXY, *unrotatedRateOfDeformationXZ;
  double *unrotatedRateOfDeformationYX, *unrotatedRateOfDeformationYY, *unrotatedRateOfDeformationYZ;
  double *unrotatedRateOfDeformationZX, *unrotatedRateOfDeformationZY, *unrotatedRateOfDeformationZZ;
  dataManager.getData(m_unrotatedRateOfDeformationXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationXX);
  dataManager.getData(m_unrotatedRateOfDeformationXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationXY);
  dataManager.getData(m_unrotatedRateOfDeformationXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationXZ);
  dataManager.getData(m_unrotatedRateOfDeformationYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationYX);
  dataManager.getData(m_unrotatedRateOfDeformationYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationYY);
  dataManager.getData(m_unrotatedRateOfDeformationYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationYZ);
  dataManager.getData(m_unrotatedRateOfDeformationZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationZX);
  dataManager.getData(m_unrotatedRateOfDeformationZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationZY);
  dataManager.getData(m_unrotatedRateOfDeformationZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformationZZ);

  CORRESPONDENCE::updateElasticCauchyStress(unrotatedRateOfDeformationXX, unrotatedRateOfDeformationXY, unrotatedRateOfDeformationXZ, 
                                            unrotatedRateOfDeformationYX, unrotatedRateOfDeformationYY, unrotatedRateOfDeformationYZ, 
                                            unrotatedRateOfDeformationZX, unrotatedRateOfDeformationZY, unrotatedRateOfDeformationZZ,
                                            cauchyStressNXX, cauchyStressNXY, cauchyStressNXZ, 
                                            cauchyStressNYX, cauchyStressNYY, cauchyStressNYZ, 
                                            cauchyStressNZX, cauchyStressNZY, cauchyStressNZZ,
                                            cauchyStressNP1XX, cauchyStressNP1XY, cauchyStressNP1XZ, 
                                            cauchyStressNP1YX, cauchyStressNP1YY, cauchyStressNP1YZ, 
                                            cauchyStressNP1ZX, cauchyStressNP1ZY, cauchyStressNP1ZZ,
                                            numOwnedPoints,
                                            m_bulkModulus,
                                            m_shearModulus,
                                            dt);
}
