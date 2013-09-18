/*! \file Peridigm_ElasticPlasticCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticPlasticCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_plastic_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticPlasticCorrespondenceMaterial::ElasticPlasticCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_yieldStress(0.0),
    m_unrotatedRateOfDeformationXXFieldId(-1), m_unrotatedRateOfDeformationXYFieldId(-1), m_unrotatedRateOfDeformationXZFieldId(-1), 
    m_unrotatedRateOfDeformationYXFieldId(-1), m_unrotatedRateOfDeformationYYFieldId(-1), m_unrotatedRateOfDeformationYZFieldId(-1), 
    m_unrotatedRateOfDeformationZXFieldId(-1), m_unrotatedRateOfDeformationZYFieldId(-1), m_unrotatedRateOfDeformationZZFieldId(-1),
    m_unrotatedCauchyStressXXFieldId(-1), m_unrotatedCauchyStressXYFieldId(-1), m_unrotatedCauchyStressXZFieldId(-1), 
    m_unrotatedCauchyStressYXFieldId(-1), m_unrotatedCauchyStressYYFieldId(-1), m_unrotatedCauchyStressYZFieldId(-1), 
    m_unrotatedCauchyStressZXFieldId(-1), m_unrotatedCauchyStressZYFieldId(-1), m_unrotatedCauchyStressZZFieldId(-1),
    m_vonMisesStressFieldId(-1), m_equivalentPlasticStrainFieldId(-1)
{
  m_yieldStress = params.get<double>("Yield Stress");

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
  m_unrotatedCauchyStressXXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressXX");
  m_unrotatedCauchyStressXYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressXY");
  m_unrotatedCauchyStressXZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressXZ");
  m_unrotatedCauchyStressYXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressYX");
  m_unrotatedCauchyStressYYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressYY");
  m_unrotatedCauchyStressYZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressYZ");
  m_unrotatedCauchyStressZXFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressZX");
  m_unrotatedCauchyStressZYFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressZY");
  m_unrotatedCauchyStressZZFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_StressZZ");

  m_vonMisesStressFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Von_Mises_Stress");

  m_equivalentPlasticStrainFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");

  m_fieldIds.push_back(m_unrotatedRateOfDeformationXXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZZFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressXXFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressXYFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressXZFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressYXFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressYYFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressYZFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressZXFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressZYFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressZZFieldId);
  m_fieldIds.push_back(m_vonMisesStressFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
}

PeridigmNS::ElasticPlasticCorrespondenceMaterial::~ElasticPlasticCorrespondenceMaterial()
{
}

void
PeridigmNS::ElasticPlasticCorrespondenceMaterial::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager) const
{

  PeridigmNS::CorrespondenceMaterial::initialize(dt,
                                                  numOwnedPoints,
                                                  ownedIDs,
                                                  neighborhoodList,
                                                  dataManager);

  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
}

void
PeridigmNS::ElasticPlasticCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               PeridigmNS::DataManager& dataManager) const
{
  double *unrotatedCauchyStressNXX, *unrotatedCauchyStressNXY, *unrotatedCauchyStressNXZ;
  double *unrotatedCauchyStressNYX, *unrotatedCauchyStressNYY, *unrotatedCauchyStressNYZ;
  double *unrotatedCauchyStressNZX, *unrotatedCauchyStressNZY, *unrotatedCauchyStressNZZ;
  dataManager.getData(m_unrotatedCauchyStressXXFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNXX);
  dataManager.getData(m_unrotatedCauchyStressXYFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNXY);
  dataManager.getData(m_unrotatedCauchyStressXZFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNXZ);
  dataManager.getData(m_unrotatedCauchyStressYXFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNYX);
  dataManager.getData(m_unrotatedCauchyStressYYFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNYY);
  dataManager.getData(m_unrotatedCauchyStressYZFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNYZ);
  dataManager.getData(m_unrotatedCauchyStressZXFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNZX);
  dataManager.getData(m_unrotatedCauchyStressZYFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNZY);
  dataManager.getData(m_unrotatedCauchyStressZZFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressNZZ);

  double *unrotatedCauchyStressNP1XX, *unrotatedCauchyStressNP1XY, *unrotatedCauchyStressNP1XZ;
  double *unrotatedCauchyStressNP1YX, *unrotatedCauchyStressNP1YY, *unrotatedCauchyStressNP1YZ;
  double *unrotatedCauchyStressNP1ZX, *unrotatedCauchyStressNP1ZY, *unrotatedCauchyStressNP1ZZ;
  dataManager.getData(m_unrotatedCauchyStressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1XX);
  dataManager.getData(m_unrotatedCauchyStressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1XY);
  dataManager.getData(m_unrotatedCauchyStressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1XZ);
  dataManager.getData(m_unrotatedCauchyStressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1YX);
  dataManager.getData(m_unrotatedCauchyStressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1YY);
  dataManager.getData(m_unrotatedCauchyStressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1YZ);
  dataManager.getData(m_unrotatedCauchyStressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1ZX);
  dataManager.getData(m_unrotatedCauchyStressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1ZY);
  dataManager.getData(m_unrotatedCauchyStressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1ZZ);

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

  double *vonMisesStress, *equivalentPlasticStrainN, *equivalentPlasticStrainNP1;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&equivalentPlasticStrainNP1);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&equivalentPlasticStrainN);

  CORRESPONDENCE::updateElasticPerfectlyPlasticCauchyStress(unrotatedRateOfDeformationXX, unrotatedRateOfDeformationXY, unrotatedRateOfDeformationXZ, 
                                                            unrotatedRateOfDeformationYX, unrotatedRateOfDeformationYY, unrotatedRateOfDeformationYZ, 
                                                            unrotatedRateOfDeformationZX, unrotatedRateOfDeformationZY, unrotatedRateOfDeformationZZ, 
                                                            unrotatedCauchyStressNXX, unrotatedCauchyStressNXY, unrotatedCauchyStressNXZ, 
                                                            unrotatedCauchyStressNYX, unrotatedCauchyStressNYY, unrotatedCauchyStressNYZ, 
                                                            unrotatedCauchyStressNZX, unrotatedCauchyStressNZY, unrotatedCauchyStressNZZ, 
                                                            unrotatedCauchyStressNP1XX, unrotatedCauchyStressNP1XY, unrotatedCauchyStressNP1XZ, 
                                                            unrotatedCauchyStressNP1YX, unrotatedCauchyStressNP1YY, unrotatedCauchyStressNP1YZ, 
                                                            unrotatedCauchyStressNP1ZX, unrotatedCauchyStressNP1ZY, unrotatedCauchyStressNP1ZZ, 
                                                            vonMisesStress,
                                                            equivalentPlasticStrainN, 
                                                            equivalentPlasticStrainNP1, 
                                                            numOwnedPoints, 
                                                            m_bulkModulus, 
                                                            m_shearModulus, 
                                                            m_yieldStress, 
                                                            dt);
}
