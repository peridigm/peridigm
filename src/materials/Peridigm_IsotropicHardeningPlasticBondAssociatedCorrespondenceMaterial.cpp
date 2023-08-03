/*! \file Peridigm_IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial.cpp */

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

#include "Peridigm_IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "isotropic_hardening_correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

PeridigmNS::IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial::IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : BondAssociatedCorrespondenceMaterial(params),
    m_yieldStress(0.0), m_hardeningStrainConstant(0.0), m_hardeningExponent(0.0),
    m_initialYieldStress(0.0), m_saturatedYieldStress(0.0), m_exponentialConstant(0.0), m_linearConstant(0.0),
    m_vonMisesStressFieldId(-1), 
    m_equivalentPlasticStrainFieldId(-1),
    m_stressTriaxialityFieldId(-1),
    m_bondLevelVonMisesStressFieldId(-1), 
    m_bondLevelEquivalentPlasticStrainFieldId(-1),
    m_bondLevelStressTriaxialityFieldId(-1)
{

  m_hardeningRule = params.get<string>("Hardening Rule");
  Teuchos::ParameterList hardeningParams = params.sublist("Hardening Parameters");
  if(m_hardeningRule == "Power-law"){
    m_yieldStress = hardeningParams.get<double>("Yield Stress");
    m_hardeningStrainConstant = hardeningParams.get<double>("Hardening Characteristic Strain");
    m_hardeningExponent = hardeningParams.get<double>("Hardening Exponent");
  }
  else if(m_hardeningRule == "Saturation type exponential"){
    m_initialYieldStress = hardeningParams.get<double>("Initial Yield Stress");
    m_saturatedYieldStress = hardeningParams.get<double>("Saturated Yield Stress");
    m_exponentialConstant = hardeningParams.get<double>("Exponential Constant");
    m_linearConstant = hardeningParams.get<double>("Linear Constant");
  }
  else {
    string invalidHardeningRule("\n**** Unrecognized hardening rule type: ");
    invalidHardeningRule += m_hardeningRule;
    invalidHardeningRule += ", must be \"Power-law\", or \"Saturation type exponential.\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, invalidHardeningRule);
  }

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_vonMisesStressFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Von_Mises_Stress");
  m_equivalentPlasticStrainFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Equivalent_Plastic_Strain");
  m_stressTriaxialityFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Stress_Triaxiality");

  m_fieldIds.push_back(m_vonMisesStressFieldId);
  m_fieldIds.push_back(m_equivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_stressTriaxialityFieldId);

  m_bondLevelVonMisesStressFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Von_Mises_Stress");
  m_bondLevelEquivalentPlasticStrainFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Equivalent_Plastic_Strain");
  m_bondLevelStressTriaxialityFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Stress_Triaxiality");

  m_fieldIds.push_back(m_bondLevelVonMisesStressFieldId);
  m_fieldIds.push_back(m_bondLevelEquivalentPlasticStrainFieldId);
  m_fieldIds.push_back(m_bondLevelStressTriaxialityFieldId);
}

PeridigmNS::IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial::~IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial()
{
}

void
PeridigmNS::IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial::initialize(const double dt,
                                                                                      const int numOwnedPoints,
                                                                                      const int* ownedIDs,
                                                                                      const int* neighborhoodList,
                                                                                      PeridigmNS::DataManager& dataManager)
{

  PeridigmNS::BondAssociatedCorrespondenceMaterial::initialize(dt,
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
PeridigmNS::IsotropicHardeningPlasticBondAssociatedCorrespondenceMaterial::computeCauchyStress(const double dt,
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

  double *stressTriaxiality;
  dataManager.getData(m_stressTriaxialityFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressTriaxiality);

  if(m_hardeningRule == "Power-law"){
    CORRESPONDENCE::updateElasticIsotropicPowerlawHardeningPlasticCauchyStress(unrotatedRateOfDeformation, 
                                                                               unrotatedCauchyStressN, 
                                                                               unrotatedCauchyStressNP1,
                                                                               vonMisesStress,
                                                                               equivalentPlasticStrainN, 
                                                                               equivalentPlasticStrainNP1, 
                                                                               stressTriaxiality,
                                                                               numOwnedPoints,
                                                                               m_bulkModulus,
                                                                               m_shearModulus,
                                                                               m_yieldStress, 
                                                                               m_hardeningStrainConstant,
                                                                               m_hardeningExponent,
                                                                               dt);
  }
  else if(m_hardeningRule == "Saturation type exponential"){
    CORRESPONDENCE::updateElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress(unrotatedRateOfDeformation, 
                                                                                            unrotatedCauchyStressN, 
                                                                                            unrotatedCauchyStressNP1,
                                                                                            vonMisesStress,
                                                                                            equivalentPlasticStrainN, 
                                                                                            equivalentPlasticStrainNP1, 
                                                                                            stressTriaxiality,
                                                                                            numOwnedPoints,
                                                                                            m_bulkModulus,
                                                                                            m_shearModulus,
                                                                                            m_initialYieldStress,
                                                                                            m_saturatedYieldStress,
                                                                                            m_exponentialConstant,
                                                                                            m_linearConstant,
                                                                                            dt);
  }

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

  if(m_hardeningRule == "Power-law"){
    CORRESPONDENCE::updateBondLevelElasticIsotropicPowerlawHardeningPlasticCauchyStress(bondLevelUnrotatedRateOfDeformationXX,
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
                                                                                        neighborhoodList,
                                                                                        numOwnedPoints,
                                                                                        m_bulkModulus,
                                                                                        m_shearModulus,
                                                                                        m_yieldStress,
                                                                                        m_hardeningStrainConstant,
                                                                                        m_hardeningExponent,
                                                                                        dt);
  }
  else if(m_hardeningRule == "Saturation type exponential"){
    CORRESPONDENCE::updateBondLevelElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress(bondLevelUnrotatedRateOfDeformationXX,
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
                                                                                                     neighborhoodList,
                                                                                                     numOwnedPoints,
                                                                                                     m_bulkModulus,
                                                                                                     m_shearModulus,
                                                                                                     m_initialYieldStress,
                                                                                                     m_saturatedYieldStress,
                                                                                                     m_exponentialConstant,
                                                                                                     m_linearConstant,
                                                                                                     dt);
  }
}
