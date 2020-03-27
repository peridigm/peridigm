/*! \file Peridigm_ElasticPlasticHypoelasticCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticPlasticHypoelasticCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_plastic_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticPlasticHypoelasticCorrespondenceMaterial::ElasticPlasticHypoelasticCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : HypoelasticCorrespondenceMaterial(params),
    m_yieldStress(0.0),
    m_unrotatedRateOfDeformationFieldId(-1),
    m_unrotatedCauchyStressFieldId(-1),
    m_vonMisesStressFieldId(-1), 
    m_equivalentPlasticStrainFieldId(-1),
    m_stressTriaxialityFieldId(-1),
    m_flyingPointFlagFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationXXFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationXYFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationXZFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationYXFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationYYFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationYZFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationZXFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationZYFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationZZFieldId(-1),
    m_bondLevelUnrotatedCauchyStressXXFieldId(-1),
    m_bondLevelUnrotatedCauchyStressXYFieldId(-1),
    m_bondLevelUnrotatedCauchyStressXZFieldId(-1),
    m_bondLevelUnrotatedCauchyStressYXFieldId(-1),
    m_bondLevelUnrotatedCauchyStressYYFieldId(-1),
    m_bondLevelUnrotatedCauchyStressYZFieldId(-1),
    m_bondLevelUnrotatedCauchyStressZXFieldId(-1),
    m_bondLevelUnrotatedCauchyStressZYFieldId(-1),
    m_bondLevelUnrotatedCauchyStressZZFieldId(-1),
    m_bondLevelVonMisesStressFieldId(-1), 
    m_bondLevelEquivalentPlasticStrainFieldId(-1),
    m_bondLevelStressTriaxialityFieldId(-1)
{

  m_yieldStress = params.get<double>("Yield Stress");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_vonMisesStressFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Von_Mises_Stress");
  m_equivalentPlasticStrainFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_stressTriaxialityFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Stress_Triaxiality");
  m_flyingPointFlagFieldId            = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flying_Point_Flag");

  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_vonMisesStressFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_stressTriaxialityFieldId);
  m_fieldIds.push_back(m_flyingPointFlagFieldId);

  m_bondLevelUnrotatedRateOfDeformationXXFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_XX");
  m_bondLevelUnrotatedRateOfDeformationXYFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_XY");
  m_bondLevelUnrotatedRateOfDeformationXZFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_XZ");
  m_bondLevelUnrotatedRateOfDeformationYXFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_YX");
  m_bondLevelUnrotatedRateOfDeformationYYFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_YY");
  m_bondLevelUnrotatedRateOfDeformationYZFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_YZ");
  m_bondLevelUnrotatedRateOfDeformationZXFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_ZX");
  m_bondLevelUnrotatedRateOfDeformationZYFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_ZY");
  m_bondLevelUnrotatedRateOfDeformationZZFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_ZZ");
  m_bondLevelUnrotatedCauchyStressXXFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_XX");
  m_bondLevelUnrotatedCauchyStressXYFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_XY");
  m_bondLevelUnrotatedCauchyStressXZFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_XZ");
  m_bondLevelUnrotatedCauchyStressYXFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_YX");
  m_bondLevelUnrotatedCauchyStressYYFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_YY");
  m_bondLevelUnrotatedCauchyStressYZFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_YZ");
  m_bondLevelUnrotatedCauchyStressZXFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_ZX");
  m_bondLevelUnrotatedCauchyStressZYFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_ZY");
  m_bondLevelUnrotatedCauchyStressZZFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_ZZ");
  m_bondLevelVonMisesStressFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Von_Mises_Stress");
  m_bondLevelEquivalentPlasticStrainFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Equivalent_Plastic_Strain");
  m_bondLevelStressTriaxialityFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Stress_Triaxiality");

  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationXXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationXYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationXZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationYXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationYYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationYZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationZXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationZYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationZZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressXXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressXYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressXZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressYXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressYYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressYZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressZXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressZYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressZZFieldId);
  m_fieldIds.push_back(m_bondLevelVonMisesStressFieldId);
  m_fieldIds.push_back(m_bondLevelEquivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_bondLevelStressTriaxialityFieldId);
}

PeridigmNS::ElasticPlasticHypoelasticCorrespondenceMaterial::~ElasticPlasticHypoelasticCorrespondenceMaterial()
{
}

void
PeridigmNS::ElasticPlasticHypoelasticCorrespondenceMaterial::initialize(const double dt,
                                                                        const int numOwnedPoints,
                                                                        const int* ownedIDs,
                                                                        const int* neighborhoodList,
                                                                        PeridigmNS::DataManager& dataManager)
{

  PeridigmNS::HypoelasticCorrespondenceMaterial::initialize(dt,
                                                            numOwnedPoints,
                                                            ownedIDs,
                                                            neighborhoodList,
                                                            dataManager);

  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_stressTriaxialityFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelEquivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelEquivalentPlasticStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelStressTriaxialityFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
}

void
PeridigmNS::ElasticPlasticHypoelasticCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                                                 const int numOwnedPoints,
                                                                                 const int* neighborhoodList,
                                                                                 PeridigmNS::DataManager& dataManager) const
  {

  // Compute the node-level stress values
  // This is only done for output (visualization) purposes
  double *unrotatedCauchyStressN;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->ExtractView(&unrotatedCauchyStressN);

  double *unrotatedCauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);

  double *unrotatedRateOfDeformation;
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  double *vonMisesStress, *equivalentPlasticStrainN, *equivalentPlasticStrainNP1;
  dataManager.getData(m_vonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&vonMisesStress);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&equivalentPlasticStrainNP1);
  dataManager.getData(m_equivalentPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&equivalentPlasticStrainN);

  double *stressTriaxiality, *flyingPointFlag;
  dataManager.getData(m_stressTriaxialityFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressTriaxiality);
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_N)->ExtractView(&flyingPointFlag);

  CORRESPONDENCE::updateElasticPerfectlyPlasticCauchyStress(unrotatedRateOfDeformation, 
                                                            unrotatedCauchyStressN, 
                                                            unrotatedCauchyStressNP1,
                                                            vonMisesStress,
                                                            equivalentPlasticStrainN, 
                                                            equivalentPlasticStrainNP1, 
                                                            stressTriaxiality,
                                                            flyingPointFlag,
                                                            numOwnedPoints,
                                                            m_bulkModulus,
                                                            m_shearModulus,
                                                            m_yieldStress, 
                                                            dt);

  // Compute the bond-level stress values
  double *bondLevelUnrotatedCauchyStressXXN, *bondLevelUnrotatedCauchyStressXYN, *bondLevelUnrotatedCauchyStressXZN;
  double *bondLevelUnrotatedCauchyStressYXN, *bondLevelUnrotatedCauchyStressYYN, *bondLevelUnrotatedCauchyStressYZN;
  double *bondLevelUnrotatedCauchyStressZXN, *bondLevelUnrotatedCauchyStressZYN, *bondLevelUnrotatedCauchyStressZZN;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressXZN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressYZN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZXN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZYN);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelUnrotatedCauchyStressZZN);

  double *bondLevelUnrotatedCauchyStressXXNP1, *bondLevelUnrotatedCauchyStressXYNP1, *bondLevelUnrotatedCauchyStressXZNP1;
  double *bondLevelUnrotatedCauchyStressYXNP1, *bondLevelUnrotatedCauchyStressYYNP1, *bondLevelUnrotatedCauchyStressYZNP1;
  double *bondLevelUnrotatedCauchyStressZXNP1, *bondLevelUnrotatedCauchyStressZYNP1, *bondLevelUnrotatedCauchyStressZZNP1;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZZNP1);

  double *bondLevelUnrotatedRateOfDeformationXX, *bondLevelUnrotatedRateOfDeformationXY, *bondLevelUnrotatedRateOfDeformationXZ;
  double *bondLevelUnrotatedRateOfDeformationYX, *bondLevelUnrotatedRateOfDeformationYY, *bondLevelUnrotatedRateOfDeformationYZ;
  double *bondLevelUnrotatedRateOfDeformationZX, *bondLevelUnrotatedRateOfDeformationZY, *bondLevelUnrotatedRateOfDeformationZZ;
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZZ);

  double *bondLevelVonMisesStress, *bondLevelEquivalentPlasticStrainN, *bondLevelEquivalentPlasticStrainNP1, *bondLevelStressTriaxiality;
  dataManager.getData(m_bondLevelVonMisesStressFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVonMisesStress);
  dataManager.getData(m_bondLevelEquivalentPlasticStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelEquivalentPlasticStrainNP1);
  dataManager.getData(m_bondLevelEquivalentPlasticStrainFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelEquivalentPlasticStrainN);
  dataManager.getData(m_bondLevelStressTriaxialityFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelStressTriaxiality);

  CORRESPONDENCE::updateBondLevelElasticPerfectlyPlasticCauchyStress(bondLevelUnrotatedRateOfDeformationXX,
                                                                     bondLevelUnrotatedRateOfDeformationXY,
                                                                     bondLevelUnrotatedRateOfDeformationXZ,
                                                                     bondLevelUnrotatedRateOfDeformationYX,
                                                                     bondLevelUnrotatedRateOfDeformationYY,
                                                                     bondLevelUnrotatedRateOfDeformationYZ,
                                                                     bondLevelUnrotatedRateOfDeformationZX,
                                                                     bondLevelUnrotatedRateOfDeformationZY,
                                                                     bondLevelUnrotatedRateOfDeformationZZ,
                                                                     bondLevelUnrotatedCauchyStressXXN, 
                                                                     bondLevelUnrotatedCauchyStressXYN, 
                                                                     bondLevelUnrotatedCauchyStressXZN, 
                                                                     bondLevelUnrotatedCauchyStressYXN, 
                                                                     bondLevelUnrotatedCauchyStressYYN, 
                                                                     bondLevelUnrotatedCauchyStressYZN, 
                                                                     bondLevelUnrotatedCauchyStressZXN, 
                                                                     bondLevelUnrotatedCauchyStressZYN, 
                                                                     bondLevelUnrotatedCauchyStressZZN, 
                                                                     bondLevelUnrotatedCauchyStressXXNP1, 
                                                                     bondLevelUnrotatedCauchyStressXYNP1, 
                                                                     bondLevelUnrotatedCauchyStressXZNP1, 
                                                                     bondLevelUnrotatedCauchyStressYXNP1, 
                                                                     bondLevelUnrotatedCauchyStressYYNP1, 
                                                                     bondLevelUnrotatedCauchyStressYZNP1, 
                                                                     bondLevelUnrotatedCauchyStressZXNP1, 
                                                                     bondLevelUnrotatedCauchyStressZYNP1, 
                                                                     bondLevelUnrotatedCauchyStressZZNP1, 
                                                                     bondLevelVonMisesStress,
                                                                     bondLevelEquivalentPlasticStrainN,
                                                                     bondLevelEquivalentPlasticStrainNP1,
                                                                     bondLevelStressTriaxiality,
                                                                     flyingPointFlag,
                                                                     neighborhoodList,
                                                                     numOwnedPoints,
                                                                     m_bulkModulus,
                                                                     m_shearModulus,
                                                                     m_yieldStress,
                                                                     dt);
}
