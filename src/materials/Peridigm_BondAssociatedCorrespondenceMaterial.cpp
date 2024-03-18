/*! \file Peridigm_BondAssociatedCorrespondenceMaterial.cpp */

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

#include "Peridigm_BondAssociatedCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include <Teuchos_Assert.hpp>
#include <Sacado.hpp> // for MPI_abort

PeridigmNS::BondAssociatedCorrespondenceMaterial::BondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0), 
    m_accuracyOrder(1), 
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_horizonFieldId(-1), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_displacementsFieldId(-1),
    m_coordinatesFieldId(-1), m_velocitiesFieldId(-1), 
    m_forceDensityFieldId(-1), m_damageFieldId(-1),
    m_bondDamageFieldId(-1), m_influenceStateFieldId(-1),
    m_gradientWeightXFieldId(-1),
    m_gradientWeightYFieldId(-1),
    m_gradientWeightZFieldId(-1),
    m_gradientWeightEvaluationFlagFieldId(-1),
    m_deformationGradientXFieldId(-1),
    m_deformationGradientYFieldId(-1),
    m_deformationGradientZFieldId(-1),
    m_deformationGradientDotXFieldId(-1),
    m_deformationGradientDotYFieldId(-1),
    m_deformationGradientDotZFieldId(-1),
    m_greenLagrangeStrainFieldId(-1),
    m_leftStretchTensorFieldId(-1),
    m_rotationTensorFieldId(-1), 
    m_unrotatedCauchyStressFieldId(-1),
    m_cauchyStressFieldId(-1), 
    m_unrotatedRateOfDeformationFieldId(-1),
    m_weightedVolumeFieldId(-1),
    m_bondLevelLeftStretchTensorXXFieldId(-1),
    m_bondLevelLeftStretchTensorXYFieldId(-1),
    m_bondLevelLeftStretchTensorXZFieldId(-1),
    m_bondLevelLeftStretchTensorYXFieldId(-1),
    m_bondLevelLeftStretchTensorYYFieldId(-1),
    m_bondLevelLeftStretchTensorYZFieldId(-1),
    m_bondLevelLeftStretchTensorZXFieldId(-1),
    m_bondLevelLeftStretchTensorZYFieldId(-1),
    m_bondLevelLeftStretchTensorZZFieldId(-1),
    m_bondLevelRotationTensorXXFieldId(-1), 
    m_bondLevelRotationTensorXYFieldId(-1), 
    m_bondLevelRotationTensorXZFieldId(-1), 
    m_bondLevelRotationTensorYXFieldId(-1), 
    m_bondLevelRotationTensorYYFieldId(-1), 
    m_bondLevelRotationTensorYZFieldId(-1), 
    m_bondLevelRotationTensorZXFieldId(-1), 
    m_bondLevelRotationTensorZYFieldId(-1), 
    m_bondLevelRotationTensorZZFieldId(-1), 
    m_bondLevelUnrotatedCauchyStressXXFieldId(-1),
    m_bondLevelUnrotatedCauchyStressXYFieldId(-1),
    m_bondLevelUnrotatedCauchyStressXZFieldId(-1),
    m_bondLevelUnrotatedCauchyStressYXFieldId(-1),
    m_bondLevelUnrotatedCauchyStressYYFieldId(-1),
    m_bondLevelUnrotatedCauchyStressYZFieldId(-1),
    m_bondLevelUnrotatedCauchyStressZXFieldId(-1),
    m_bondLevelUnrotatedCauchyStressZYFieldId(-1),
    m_bondLevelUnrotatedCauchyStressZZFieldId(-1),
    m_bondLevelCauchyStressXXFieldId(-1), 
    m_bondLevelCauchyStressXYFieldId(-1), 
    m_bondLevelCauchyStressXZFieldId(-1), 
    m_bondLevelCauchyStressYXFieldId(-1), 
    m_bondLevelCauchyStressYYFieldId(-1), 
    m_bondLevelCauchyStressYZFieldId(-1), 
    m_bondLevelCauchyStressZXFieldId(-1), 
    m_bondLevelCauchyStressZYFieldId(-1), 
    m_bondLevelCauchyStressZZFieldId(-1), 
    m_bondLevelUnrotatedRateOfDeformationXXFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationXYFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationXZFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationYXFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationYYFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationYZFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationZXFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationZYFieldId(-1),
    m_bondLevelUnrotatedRateOfDeformationZZFieldId(-1),
    m_bondLevelPiolaStressXXFieldId(-1), 
    m_bondLevelPiolaStressXYFieldId(-1), 
    m_bondLevelPiolaStressXZFieldId(-1), 
    m_bondLevelPiolaStressYXFieldId(-1), 
    m_bondLevelPiolaStressYYFieldId(-1), 
    m_bondLevelPiolaStressYZFieldId(-1), 
    m_bondLevelPiolaStressZXFieldId(-1), 
    m_bondLevelPiolaStressZYFieldId(-1), 
    m_bondLevelPiolaStressZZFieldId(-1), 
    m_bondLevelDeformationGradientInvXXFieldId(-1),
    m_bondLevelDeformationGradientInvXYFieldId(-1),
    m_bondLevelDeformationGradientInvXZFieldId(-1),
    m_bondLevelDeformationGradientInvYXFieldId(-1),
    m_bondLevelDeformationGradientInvYYFieldId(-1),
    m_bondLevelDeformationGradientInvYZFieldId(-1),
    m_bondLevelDeformationGradientInvZXFieldId(-1),
    m_bondLevelDeformationGradientInvZYFieldId(-1),
    m_bondLevelDeformationGradientInvZZFieldId(-1),
    m_bondLevelJacobianDeterminantFieldId(-1),
    m_stressIntegralFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");

  if(params.isParameter("Gradient Order Of Accuracy")){
    m_accuracyOrder = params.get<int>("Gradient Order Of Accuracy");
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticBondAssociatedCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticBondAssociatedCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticBondAssociatedCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_horizonFieldId                    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_displacementsFieldId              = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_damageFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_influenceStateFieldId             = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::CONSTANT, "Influence_State");
  m_gradientWeightXFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_X");
  m_gradientWeightYFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_Y");
  m_gradientWeightZFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_Z");
  m_gradientWeightEvaluationFlagFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_Evaluation_Flag");
  m_deformationGradientXFieldId       = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Deformation_Gradient_X");
  m_deformationGradientYFieldId       = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Deformation_Gradient_Y");
  m_deformationGradientZFieldId       = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Deformation_Gradient_Z");
  m_deformationGradientDotXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Deformation_Gradient_Dot_X");
  m_deformationGradientDotYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Deformation_Gradient_Dot_Y");
  m_deformationGradientDotZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Deformation_Gradient_Dot_Z");
  m_greenLagrangeStrainFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Green_Lagrange_Strain");
  m_leftStretchTensorFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_weightedVolumeFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
  
  // Currently BOND Field is limited to SCALAR. The reason is to avoid
  // allocating too much memory for the bond data. 
  // Perhaps later a Multi-Vector Bond Field can be added to avoid using 
  // multiple field managers here.
  m_bondLevelLeftStretchTensorXXFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_XX");
  m_bondLevelLeftStretchTensorXYFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_XY");
  m_bondLevelLeftStretchTensorXZFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_XZ");
  m_bondLevelLeftStretchTensorYXFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_YX");
  m_bondLevelLeftStretchTensorYYFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_YY");
  m_bondLevelLeftStretchTensorYZFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_YZ");
  m_bondLevelLeftStretchTensorZXFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_ZX");
  m_bondLevelLeftStretchTensorZYFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_ZY");
  m_bondLevelLeftStretchTensorZZFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor_ZZ");
  m_bondLevelRotationTensorXXFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_XX");
  m_bondLevelRotationTensorXYFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_XY");
  m_bondLevelRotationTensorXZFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_XZ");
  m_bondLevelRotationTensorYXFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_YX");
  m_bondLevelRotationTensorYYFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_YY");
  m_bondLevelRotationTensorYZFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_YZ");
  m_bondLevelRotationTensorZXFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_ZX");
  m_bondLevelRotationTensorZYFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_ZY");
  m_bondLevelRotationTensorZZFieldId             = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Rotation_Tensor_ZZ");
  m_bondLevelUnrotatedCauchyStressXXFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_XX");
  m_bondLevelUnrotatedCauchyStressXYFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_XY");
  m_bondLevelUnrotatedCauchyStressXZFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_XZ");
  m_bondLevelUnrotatedCauchyStressYXFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_YX");
  m_bondLevelUnrotatedCauchyStressYYFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_YY");
  m_bondLevelUnrotatedCauchyStressYZFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_YZ");
  m_bondLevelUnrotatedCauchyStressZXFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_ZX");
  m_bondLevelUnrotatedCauchyStressZYFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_ZY");
  m_bondLevelUnrotatedCauchyStressZZFieldId      = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress_ZZ");
  m_bondLevelCauchyStressXXFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_XX");
  m_bondLevelCauchyStressXYFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_XY");
  m_bondLevelCauchyStressXZFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_XZ");
  m_bondLevelCauchyStressYXFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_YX");
  m_bondLevelCauchyStressYYFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_YY");
  m_bondLevelCauchyStressYZFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_YZ");
  m_bondLevelCauchyStressZXFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_ZX");
  m_bondLevelCauchyStressZYFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_ZY");
  m_bondLevelCauchyStressZZFieldId               = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Cauchy_Stress_ZZ");
  m_bondLevelUnrotatedRateOfDeformationXXFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_XX");
  m_bondLevelUnrotatedRateOfDeformationXYFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_XY");
  m_bondLevelUnrotatedRateOfDeformationXZFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_XZ");
  m_bondLevelUnrotatedRateOfDeformationYXFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_YX");
  m_bondLevelUnrotatedRateOfDeformationYYFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_YY");
  m_bondLevelUnrotatedRateOfDeformationYZFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_YZ");
  m_bondLevelUnrotatedRateOfDeformationZXFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_ZX");
  m_bondLevelUnrotatedRateOfDeformationZYFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_ZY");
  m_bondLevelUnrotatedRateOfDeformationZZFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation_ZZ");
  m_bondLevelPiolaStressXXFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_XX");
  m_bondLevelPiolaStressXYFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_XY");
  m_bondLevelPiolaStressXZFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_XZ");
  m_bondLevelPiolaStressYXFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_YX");
  m_bondLevelPiolaStressYYFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_YY");
  m_bondLevelPiolaStressYZFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_YZ");
  m_bondLevelPiolaStressZXFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_ZX");
  m_bondLevelPiolaStressZYFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_ZY");
  m_bondLevelPiolaStressZZFieldId                = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Piola_Stress_ZZ");
  m_bondLevelDeformationGradientInvXXFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_XX");
  m_bondLevelDeformationGradientInvXYFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_XY");
  m_bondLevelDeformationGradientInvXZFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_XZ");
  m_bondLevelDeformationGradientInvYXFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_YX");
  m_bondLevelDeformationGradientInvYYFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_YY");
  m_bondLevelDeformationGradientInvYZFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_YZ");
  m_bondLevelDeformationGradientInvZXFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_ZX");
  m_bondLevelDeformationGradientInvZYFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_ZY");
  m_bondLevelDeformationGradientInvZZFieldId     = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_Gradient_Inv_ZZ");
  m_bondLevelJacobianDeterminantFieldId          = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Jacobian_Determinant");
  m_stressIntegralFieldId                        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Stress_Integral");

  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_displacementsFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_influenceStateFieldId);
  m_fieldIds.push_back(m_gradientWeightXFieldId);
  m_fieldIds.push_back(m_gradientWeightYFieldId);
  m_fieldIds.push_back(m_gradientWeightZFieldId);
  m_fieldIds.push_back(m_gradientWeightEvaluationFlagFieldId);
  m_fieldIds.push_back(m_deformationGradientXFieldId);
  m_fieldIds.push_back(m_deformationGradientYFieldId);
  m_fieldIds.push_back(m_deformationGradientZFieldId);
  m_fieldIds.push_back(m_deformationGradientDotXFieldId);
  m_fieldIds.push_back(m_deformationGradientDotYFieldId);
  m_fieldIds.push_back(m_deformationGradientDotZFieldId);
  m_fieldIds.push_back(m_greenLagrangeStrainFieldId);
  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);

  m_fieldIds.push_back(m_bondLevelLeftStretchTensorXXFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorXYFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorXZFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorYXFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorYYFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorYZFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorZXFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorZYFieldId);
  m_fieldIds.push_back(m_bondLevelLeftStretchTensorZZFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorXXFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorXYFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorXZFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorYXFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorYYFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorYZFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorZXFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorZYFieldId);
  m_fieldIds.push_back(m_bondLevelRotationTensorZZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressXXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressXYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressXZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressYXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressYYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressYZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressZXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressZYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedCauchyStressZZFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressXXFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressXYFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressXZFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressYXFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressYYFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressYZFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressZXFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressZYFieldId);
  m_fieldIds.push_back(m_bondLevelCauchyStressZZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationXXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationXYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationXZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationYXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationYYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationYZFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationZXFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationZYFieldId);
  m_fieldIds.push_back(m_bondLevelUnrotatedRateOfDeformationZZFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressXXFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressXYFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressXZFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressYXFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressYYFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressYZFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressZXFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressZYFieldId);
  m_fieldIds.push_back(m_bondLevelPiolaStressZZFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvXXFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvXYFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvXZFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvYXFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvYYFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvYZFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvZXFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvZYFieldId);
  m_fieldIds.push_back(m_bondLevelDeformationGradientInvZZFieldId);
  m_fieldIds.push_back(m_bondLevelJacobianDeterminantFieldId);
  m_fieldIds.push_back(m_stressIntegralFieldId);
}

PeridigmNS::BondAssociatedCorrespondenceMaterial::~BondAssociatedCorrespondenceMaterial()
{
}

void
PeridigmNS::BondAssociatedCorrespondenceMaterial::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager)
{
  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_gradientWeightEvaluationFlagFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientDotXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientDotYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientDotZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *leftStretchTensorN;
  double *leftStretchTensorNP1;
  double *rotationTensorN;
  double *rotationTensorNP1;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);

  //Initialize the left stretch and rotation tenor to the identity matrix
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(leftStretchTensorNP1, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(rotationTensorNP1, numOwnedPoints);

  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  dataManager.getData(m_bondLevelLeftStretchTensorXXFieldId, PeridigmField::STEP_N)->PutScalar(1.0);
  dataManager.getData(m_bondLevelLeftStretchTensorXYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorXZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorYXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorYYFieldId, PeridigmField::STEP_N)->PutScalar(1.0);
  dataManager.getData(m_bondLevelLeftStretchTensorYZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorZXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorZYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorZZFieldId, PeridigmField::STEP_N)->PutScalar(1.0);

  dataManager.getData(m_bondLevelLeftStretchTensorXXFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);
  dataManager.getData(m_bondLevelLeftStretchTensorXYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorXZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorYXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorYYFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);
  dataManager.getData(m_bondLevelLeftStretchTensorYZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorZXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorZYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelLeftStretchTensorZZFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);

  dataManager.getData(m_bondLevelRotationTensorXXFieldId, PeridigmField::STEP_N)->PutScalar(1.0);
  dataManager.getData(m_bondLevelRotationTensorXYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorXZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorYXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorYYFieldId, PeridigmField::STEP_N)->PutScalar(1.0);
  dataManager.getData(m_bondLevelRotationTensorYZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorZXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorZYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorZZFieldId, PeridigmField::STEP_N)->PutScalar(1.0);

  dataManager.getData(m_bondLevelRotationTensorXXFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);
  dataManager.getData(m_bondLevelRotationTensorXYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorXZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorYXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorYYFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);
  dataManager.getData(m_bondLevelRotationTensorYZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorZXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorZYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelRotationTensorZZFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);

  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);

  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(m_bondLevelCauchyStressXXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressXYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressXZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressYXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressYYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressYZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressZXFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressZYFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressZZFieldId, PeridigmField::STEP_N)->PutScalar(0.0);

  dataManager.getData(m_bondLevelCauchyStressXXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressXYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressXZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressYXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressYYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressYZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressZXFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressZYFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondLevelCauchyStressZZFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(m_bondLevelPiolaStressXXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressXYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressXZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressYXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressYYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressYZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressZXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressZYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelPiolaStressZZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  dataManager.getData(m_bondLevelDeformationGradientInvXXFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_bondLevelDeformationGradientInvXYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDeformationGradientInvXZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDeformationGradientInvYXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDeformationGradientInvYYFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_bondLevelDeformationGradientInvYZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDeformationGradientInvZXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDeformationGradientInvZYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_bondLevelDeformationGradientInvZZFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);

  dataManager.getData(m_bondLevelJacobianDeterminantFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);

  dataManager.getData(m_stressIntegralFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  dataManager.getData(m_damageFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *horizon, *volume; 
  double *weightedVolume;
  double *modelCoordinates;
  double *gradientWeightX, *gradientWeightY, *gradientWeightZ;
  double *gradientWeightEvaluationFlag;
  double *damage, *bondDamage;
  double *influenceState;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightX);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightY);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightZ);
  dataManager.getData(m_gradientWeightEvaluationFlagFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightEvaluationFlag);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_influenceStateFieldId, PeridigmField::STEP_NONE)->ExtractView(&influenceState);

  int gradientWeightReturnCode = CORRESPONDENCE::computeLagrangianGradientWeights(horizon,
                                                                                  modelCoordinates,
                                                                                  volume,
                                                                                  gradientWeightX,
                                                                                  gradientWeightY,
                                                                                  gradientWeightZ,
                                                                                  gradientWeightEvaluationFlag,
                                                                                  bondDamage,
                                                                                  influenceState,
                                                                                  m_accuracyOrder,
                                                                                  neighborhoodList,
                                                                                  numOwnedPoints);
          
  string gradientWeightErrorMessage =
    "**** Error:  BondAssociatedCorrespondenceMaterial::precompute() failed to compute gradient weights.\n";
  gradientWeightErrorMessage +=
    "****         Possible scenarios: 1) The horizon is too small, or 2) Too much damage around some points.\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(gradientWeightReturnCode != 0, gradientWeightErrorMessage);

  // Compute the weighted volume
  CORRESPONDENCE::computeWeightedVolume(volume,
                                        weightedVolume,
                                        influenceState,
                                        neighborhoodList,
                                        numOwnedPoints);
}

void
PeridigmNS::BondAssociatedCorrespondenceMaterial::computeForce(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* ownedIDs,
                                                               const int* neighborhoodList,
                                                               PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces 
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *horizon, *volume, *modelCoordinates, *coordinates, *velocities;
  double *gradientWeightX, *gradientWeightY, *gradientWeightZ;
  double *deformationGradientX, *deformationGradientY, *deformationGradientZ;
  double *deformationGradientDotX, *deformationGradientDotY, *deformationGradientDotZ;
  double *bondLevelDeformationGradientInvXX, *bondLevelDeformationGradientInvXY, *bondLevelDeformationGradientInvXZ;
  double *bondLevelDeformationGradientInvYX, *bondLevelDeformationGradientInvYY, *bondLevelDeformationGradientInvYZ;
  double *bondLevelDeformationGradientInvZX, *bondLevelDeformationGradientInvZY, *bondLevelDeformationGradientInvZZ;
  double *bondLevelJacobianDeterminant;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightX);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightY);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightZ);
  dataManager.getData(m_deformationGradientXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientX);
  dataManager.getData(m_deformationGradientYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientY);
  dataManager.getData(m_deformationGradientZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZ);
  dataManager.getData(m_deformationGradientDotXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientDotX);
  dataManager.getData(m_deformationGradientDotYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientDotY);
  dataManager.getData(m_deformationGradientDotZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientDotZ);
  dataManager.getData(m_bondLevelDeformationGradientInvXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvXX);
  dataManager.getData(m_bondLevelDeformationGradientInvXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvXY);
  dataManager.getData(m_bondLevelDeformationGradientInvXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvXZ);
  dataManager.getData(m_bondLevelDeformationGradientInvYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvYX);
  dataManager.getData(m_bondLevelDeformationGradientInvYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvYY);
  dataManager.getData(m_bondLevelDeformationGradientInvYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvYZ);
  dataManager.getData(m_bondLevelDeformationGradientInvZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvZX);
  dataManager.getData(m_bondLevelDeformationGradientInvZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvZY);
  dataManager.getData(m_bondLevelDeformationGradientInvZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelDeformationGradientInvZZ);
  dataManager.getData(m_bondLevelJacobianDeterminantFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelJacobianDeterminant);
    
  double *influenceState;
  dataManager.getData(m_influenceStateFieldId, PeridigmField::STEP_NONE)->ExtractView(&influenceState);

  double *greenLagrangeStrainN, *greenLagrangeStrainNP1;
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_N)->ExtractView(&greenLagrangeStrainN);
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_NP1)->ExtractView(&greenLagrangeStrainNP1);

  // This step is done only for visualization purposes.
  CORRESPONDENCE::updateGreenLagrangeStrain(deformationGradientX,
                                            deformationGradientY,
                                            deformationGradientZ,
                                            deformationGradientDotX,
                                            deformationGradientDotY,
                                            deformationGradientDotZ,
                                            greenLagrangeStrainN,
                                            greenLagrangeStrainNP1,
                                            numOwnedPoints,
                                            dt);

  double *leftStretchTensorN, *leftStretchTensorNP1, *rotationTensorN, *rotationTensorNP1, *unrotatedRateOfDeformation;
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->ExtractView(&unrotatedRateOfDeformation);

  double *bondLevelLeftStretchTensorXXN, *bondLevelLeftStretchTensorXXNP1, *bondLevelRotationTensorXXN, *bondLevelRotationTensorXXNP1, *bondLevelUnrotatedRateOfDeformationXX;
  double *bondLevelLeftStretchTensorXYN, *bondLevelLeftStretchTensorXYNP1, *bondLevelRotationTensorXYN, *bondLevelRotationTensorXYNP1, *bondLevelUnrotatedRateOfDeformationXY;
  double *bondLevelLeftStretchTensorXZN, *bondLevelLeftStretchTensorXZNP1, *bondLevelRotationTensorXZN, *bondLevelRotationTensorXZNP1, *bondLevelUnrotatedRateOfDeformationXZ;
  double *bondLevelLeftStretchTensorYXN, *bondLevelLeftStretchTensorYXNP1, *bondLevelRotationTensorYXN, *bondLevelRotationTensorYXNP1, *bondLevelUnrotatedRateOfDeformationYX;
  double *bondLevelLeftStretchTensorYYN, *bondLevelLeftStretchTensorYYNP1, *bondLevelRotationTensorYYN, *bondLevelRotationTensorYYNP1, *bondLevelUnrotatedRateOfDeformationYY;
  double *bondLevelLeftStretchTensorYZN, *bondLevelLeftStretchTensorYZNP1, *bondLevelRotationTensorYZN, *bondLevelRotationTensorYZNP1, *bondLevelUnrotatedRateOfDeformationYZ;
  double *bondLevelLeftStretchTensorZXN, *bondLevelLeftStretchTensorZXNP1, *bondLevelRotationTensorZXN, *bondLevelRotationTensorZXNP1, *bondLevelUnrotatedRateOfDeformationZX;
  double *bondLevelLeftStretchTensorZYN, *bondLevelLeftStretchTensorZYNP1, *bondLevelRotationTensorZYN, *bondLevelRotationTensorZYNP1, *bondLevelUnrotatedRateOfDeformationZY;
  double *bondLevelLeftStretchTensorZZN, *bondLevelLeftStretchTensorZZNP1, *bondLevelRotationTensorZZN, *bondLevelRotationTensorZZNP1, *bondLevelUnrotatedRateOfDeformationZZ;
  dataManager.getData(m_bondLevelLeftStretchTensorXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorXXN);
  dataManager.getData(m_bondLevelLeftStretchTensorXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorXYN);
  dataManager.getData(m_bondLevelLeftStretchTensorXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorXZN);
  dataManager.getData(m_bondLevelLeftStretchTensorYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorYXN);
  dataManager.getData(m_bondLevelLeftStretchTensorYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorYYN);
  dataManager.getData(m_bondLevelLeftStretchTensorYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorYZN);
  dataManager.getData(m_bondLevelLeftStretchTensorZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorZXN);
  dataManager.getData(m_bondLevelLeftStretchTensorZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorZYN);
  dataManager.getData(m_bondLevelLeftStretchTensorZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelLeftStretchTensorZZN);
  dataManager.getData(m_bondLevelRotationTensorXXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorXXN);
  dataManager.getData(m_bondLevelRotationTensorXYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorXYN);
  dataManager.getData(m_bondLevelRotationTensorXZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorXZN);
  dataManager.getData(m_bondLevelRotationTensorYXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorYXN);
  dataManager.getData(m_bondLevelRotationTensorYYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorYYN);
  dataManager.getData(m_bondLevelRotationTensorYZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorYZN);
  dataManager.getData(m_bondLevelRotationTensorZXFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorZXN);
  dataManager.getData(m_bondLevelRotationTensorZYFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorZYN);
  dataManager.getData(m_bondLevelRotationTensorZZFieldId, PeridigmField::STEP_N)->ExtractView(&bondLevelRotationTensorZZN);
  dataManager.getData(m_bondLevelLeftStretchTensorXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorXXNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorXYNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorXZNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorYXNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorYYNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorYZNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorZXNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorZYNP1);
  dataManager.getData(m_bondLevelLeftStretchTensorZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelLeftStretchTensorZZNP1);
  dataManager.getData(m_bondLevelRotationTensorXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorXXNP1);
  dataManager.getData(m_bondLevelRotationTensorXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorXYNP1);
  dataManager.getData(m_bondLevelRotationTensorXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorXZNP1);
  dataManager.getData(m_bondLevelRotationTensorYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorYXNP1);
  dataManager.getData(m_bondLevelRotationTensorYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorYYNP1);
  dataManager.getData(m_bondLevelRotationTensorYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorYZNP1);
  dataManager.getData(m_bondLevelRotationTensorZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorZXNP1);
  dataManager.getData(m_bondLevelRotationTensorZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorZYNP1);
  dataManager.getData(m_bondLevelRotationTensorZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelRotationTensorZZNP1);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationXZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationYZ);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZX);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZY);
  dataManager.getData(m_bondLevelUnrotatedRateOfDeformationZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelUnrotatedRateOfDeformationZZ);

  // Compute left stretch tensor, rotation tensor, and unrotated rate-of-deformation.
  // Performs a polar decomposition via Flanagan & Taylor (1987) algorithm.
  
  // Compute node-level values
  // this is only done for output (visualization) purposes
  int nodeLevelRotationTensorReturnCode = CORRESPONDENCE::computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor(deformationGradientX,
                                                                                                                      deformationGradientY,
                                                                                                                      deformationGradientZ,
                                                                                                                      deformationGradientDotX,
                                                                                                                      deformationGradientDotY,
                                                                                                                      deformationGradientDotZ,
                                                                                                                      leftStretchTensorN,
                                                                                                                      rotationTensorN,
                                                                                                                      leftStretchTensorNP1,
                                                                                                                      rotationTensorNP1,
                                                                                                                      unrotatedRateOfDeformation,
                                                                                                                      numOwnedPoints, 
                                                                                                                      dt);
  string nodeLevelRotationTensorErrorMessage =
    "**** Error:  BondAssociatedCorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
  nodeLevelRotationTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_TERMINATION(nodeLevelRotationTensorReturnCode != 0, nodeLevelRotationTensorErrorMessage);

  // Compute bond-level values
  int bondLevelRotationTensorReturnCode = CORRESPONDENCE::computeBondLevelUnrotatedRateOfDeformationAndRotationTensor(modelCoordinates,
                                                                                                                      coordinates,
                                                                                                                      velocities,
                                                                                                                      deformationGradientX,
                                                                                                                      deformationGradientY,
                                                                                                                      deformationGradientZ,
                                                                                                                      deformationGradientDotX,
                                                                                                                      deformationGradientDotY,
                                                                                                                      deformationGradientDotZ,
                                                                                                                      bondLevelDeformationGradientInvXX,
                                                                                                                      bondLevelDeformationGradientInvXY,
                                                                                                                      bondLevelDeformationGradientInvXZ,
                                                                                                                      bondLevelDeformationGradientInvYX,
                                                                                                                      bondLevelDeformationGradientInvYY,
                                                                                                                      bondLevelDeformationGradientInvYZ,
                                                                                                                      bondLevelDeformationGradientInvZX,
                                                                                                                      bondLevelDeformationGradientInvZY,
                                                                                                                      bondLevelDeformationGradientInvZZ,
                                                                                                                      bondLevelJacobianDeterminant,
                                                                                                                      bondLevelLeftStretchTensorXXN,
                                                                                                                      bondLevelLeftStretchTensorXYN,
                                                                                                                      bondLevelLeftStretchTensorXZN,
                                                                                                                      bondLevelLeftStretchTensorYXN,
                                                                                                                      bondLevelLeftStretchTensorYYN,
                                                                                                                      bondLevelLeftStretchTensorYZN,
                                                                                                                      bondLevelLeftStretchTensorZXN,
                                                                                                                      bondLevelLeftStretchTensorZYN,
                                                                                                                      bondLevelLeftStretchTensorZZN,
                                                                                                                      bondLevelRotationTensorXXN,
                                                                                                                      bondLevelRotationTensorXYN,
                                                                                                                      bondLevelRotationTensorXZN,
                                                                                                                      bondLevelRotationTensorYXN,
                                                                                                                      bondLevelRotationTensorYYN,
                                                                                                                      bondLevelRotationTensorYZN,
                                                                                                                      bondLevelRotationTensorZXN,
                                                                                                                      bondLevelRotationTensorZYN,
                                                                                                                      bondLevelRotationTensorZZN,
                                                                                                                      bondLevelLeftStretchTensorXXNP1,
                                                                                                                      bondLevelLeftStretchTensorXYNP1,
                                                                                                                      bondLevelLeftStretchTensorXZNP1,
                                                                                                                      bondLevelLeftStretchTensorYXNP1,
                                                                                                                      bondLevelLeftStretchTensorYYNP1,
                                                                                                                      bondLevelLeftStretchTensorYZNP1,
                                                                                                                      bondLevelLeftStretchTensorZXNP1,
                                                                                                                      bondLevelLeftStretchTensorZYNP1,
                                                                                                                      bondLevelLeftStretchTensorZZNP1,
                                                                                                                      bondLevelRotationTensorXXNP1,
                                                                                                                      bondLevelRotationTensorXYNP1,
                                                                                                                      bondLevelRotationTensorXZNP1,
                                                                                                                      bondLevelRotationTensorYXNP1,
                                                                                                                      bondLevelRotationTensorYYNP1,
                                                                                                                      bondLevelRotationTensorYZNP1,
                                                                                                                      bondLevelRotationTensorZXNP1,
                                                                                                                      bondLevelRotationTensorZYNP1,
                                                                                                                      bondLevelRotationTensorZZNP1,
                                                                                                                      bondLevelUnrotatedRateOfDeformationXX,
                                                                                                                      bondLevelUnrotatedRateOfDeformationXY,
                                                                                                                      bondLevelUnrotatedRateOfDeformationXZ,
                                                                                                                      bondLevelUnrotatedRateOfDeformationYX,
                                                                                                                      bondLevelUnrotatedRateOfDeformationYY,
                                                                                                                      bondLevelUnrotatedRateOfDeformationYZ,
                                                                                                                      bondLevelUnrotatedRateOfDeformationZX,
                                                                                                                      bondLevelUnrotatedRateOfDeformationZY,
                                                                                                                      bondLevelUnrotatedRateOfDeformationZZ,
                                                                                                                      influenceState,
                                                                                                                      neighborhoodList, 
                                                                                                                      numOwnedPoints, 
                                                                                                                      dt);
  string bondLevelRotationTensorErrorMessage =
    "**** Error:  BondAssociatedCorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
  bondLevelRotationTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_TERMINATION(bondLevelRotationTensorReturnCode != 0, bondLevelRotationTensorErrorMessage);

  // Evaluate the Cauchy stress using the routine implemented in the derived class (specific correspondence material model)
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  // computeCauchyStress() typically uses the following fields which are accessed via the DataManager:
  //   Input:  unrotated rate-of-deformation tensor
  //   Input:  unrotated Cauchy stress at step N
  //   Input:  internal state data (managed in the derived class)
  //   Output: unrotated Cauchy stress at step N+1
  computeCauchyStress(dt, numOwnedPoints, neighborhoodList, dataManager);

  // rotate back to the Eulerian frame
  // perform on the node level 
  // this is only done for output (visualization) purposes
  double *unrotatedCauchyStressNP1, *cauchyStressNP1;
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&unrotatedCauchyStressNP1);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->ExtractView(&cauchyStressNP1);

  CORRESPONDENCE::rotateCauchyStress(rotationTensorNP1,
                                     unrotatedCauchyStressNP1,
                                     cauchyStressNP1,
                                     numOwnedPoints);
  
  // perform on the bond level 
  double *bondLevelUnrotatedCauchyStressXXNP1, *bondLevelCauchyStressXXNP1;
  double *bondLevelUnrotatedCauchyStressXYNP1, *bondLevelCauchyStressXYNP1;
  double *bondLevelUnrotatedCauchyStressXZNP1, *bondLevelCauchyStressXZNP1;
  double *bondLevelUnrotatedCauchyStressYXNP1, *bondLevelCauchyStressYXNP1;
  double *bondLevelUnrotatedCauchyStressYYNP1, *bondLevelCauchyStressYYNP1;
  double *bondLevelUnrotatedCauchyStressYZNP1, *bondLevelCauchyStressYZNP1;
  double *bondLevelUnrotatedCauchyStressZXNP1, *bondLevelCauchyStressZXNP1;
  double *bondLevelUnrotatedCauchyStressZYNP1, *bondLevelCauchyStressZYNP1;
  double *bondLevelUnrotatedCauchyStressZZNP1, *bondLevelCauchyStressZZNP1;
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressXZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressYZNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZXNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZYNP1);
  dataManager.getData(m_bondLevelUnrotatedCauchyStressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelUnrotatedCauchyStressZZNP1);
  dataManager.getData(m_bondLevelCauchyStressXXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressXXNP1);
  dataManager.getData(m_bondLevelCauchyStressXYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressXYNP1);
  dataManager.getData(m_bondLevelCauchyStressXZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressXZNP1);
  dataManager.getData(m_bondLevelCauchyStressYXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressYXNP1);
  dataManager.getData(m_bondLevelCauchyStressYYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressYYNP1);
  dataManager.getData(m_bondLevelCauchyStressYZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressYZNP1);
  dataManager.getData(m_bondLevelCauchyStressZXFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressZXNP1);
  dataManager.getData(m_bondLevelCauchyStressZYFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressZYNP1);
  dataManager.getData(m_bondLevelCauchyStressZZFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondLevelCauchyStressZZNP1);

  CORRESPONDENCE::rotateBondLevelCauchyStress(bondLevelRotationTensorXXNP1,
                                              bondLevelRotationTensorXYNP1,
                                              bondLevelRotationTensorXZNP1,
                                              bondLevelRotationTensorYXNP1,
                                              bondLevelRotationTensorYYNP1,
                                              bondLevelRotationTensorYZNP1,
                                              bondLevelRotationTensorZXNP1,
                                              bondLevelRotationTensorZYNP1,
                                              bondLevelRotationTensorZZNP1,
                                              bondLevelUnrotatedCauchyStressXXNP1,
                                              bondLevelUnrotatedCauchyStressXYNP1,
                                              bondLevelUnrotatedCauchyStressXZNP1,
                                              bondLevelUnrotatedCauchyStressYXNP1,
                                              bondLevelUnrotatedCauchyStressYYNP1,
                                              bondLevelUnrotatedCauchyStressYZNP1,
                                              bondLevelUnrotatedCauchyStressZXNP1,
                                              bondLevelUnrotatedCauchyStressZYNP1,
                                              bondLevelUnrotatedCauchyStressZZNP1,
                                              bondLevelCauchyStressXXNP1,
                                              bondLevelCauchyStressXYNP1,
                                              bondLevelCauchyStressXZNP1,
                                              bondLevelCauchyStressYXNP1,
                                              bondLevelCauchyStressYYNP1,
                                              bondLevelCauchyStressYZNP1,
                                              bondLevelCauchyStressZXNP1,
                                              bondLevelCauchyStressZYNP1,
                                              bondLevelCauchyStressZZNP1,
                                              neighborhoodList,
                                              numOwnedPoints);

  double *bondLevelPiolaStressXX;
  double *bondLevelPiolaStressXY;
  double *bondLevelPiolaStressXZ;
  double *bondLevelPiolaStressYX;
  double *bondLevelPiolaStressYY;
  double *bondLevelPiolaStressYZ;
  double *bondLevelPiolaStressZX;
  double *bondLevelPiolaStressZY;
  double *bondLevelPiolaStressZZ;
  dataManager.getData(m_bondLevelPiolaStressXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressXX);
  dataManager.getData(m_bondLevelPiolaStressXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressXY);
  dataManager.getData(m_bondLevelPiolaStressXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressXZ);
  dataManager.getData(m_bondLevelPiolaStressYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressYX);
  dataManager.getData(m_bondLevelPiolaStressYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressYY);
  dataManager.getData(m_bondLevelPiolaStressYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressYZ);
  dataManager.getData(m_bondLevelPiolaStressZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressZX);
  dataManager.getData(m_bondLevelPiolaStressZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressZY);
  dataManager.getData(m_bondLevelPiolaStressZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelPiolaStressZZ);

  CORRESPONDENCE::computeBondLevelPiolaStress(bondLevelJacobianDeterminant,
                                              bondLevelCauchyStressXXNP1,
                                              bondLevelCauchyStressXYNP1,
                                              bondLevelCauchyStressXZNP1,
                                              bondLevelCauchyStressYXNP1,
                                              bondLevelCauchyStressYYNP1,
                                              bondLevelCauchyStressYZNP1,
                                              bondLevelCauchyStressZXNP1,
                                              bondLevelCauchyStressZYNP1,
                                              bondLevelCauchyStressZZNP1,
                                              bondLevelDeformationGradientInvXX,
                                              bondLevelDeformationGradientInvXY,
                                              bondLevelDeformationGradientInvXZ,
                                              bondLevelDeformationGradientInvYX,
                                              bondLevelDeformationGradientInvYY,
                                              bondLevelDeformationGradientInvYZ,
                                              bondLevelDeformationGradientInvZX,
                                              bondLevelDeformationGradientInvZY,
                                              bondLevelDeformationGradientInvZZ,
                                              bondLevelPiolaStressXX,
                                              bondLevelPiolaStressXY,
                                              bondLevelPiolaStressXZ,
                                              bondLevelPiolaStressYX,
                                              bondLevelPiolaStressYY,
                                              bondLevelPiolaStressYZ,
                                              bondLevelPiolaStressZX,
                                              bondLevelPiolaStressZY,
                                              bondLevelPiolaStressZZ,
                                              neighborhoodList,
                                              numOwnedPoints);

  double *stressIntegral;
  dataManager.getData(m_stressIntegralFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressIntegral);

  double *weightedVolume;
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);

  CORRESPONDENCE::computeStressIntegral(volume,
                                        weightedVolume,
                                        modelCoordinates,
                                        bondLevelPiolaStressXX,
                                        bondLevelPiolaStressXY,
                                        bondLevelPiolaStressXZ,
                                        bondLevelPiolaStressYX,
                                        bondLevelPiolaStressYY,
                                        bondLevelPiolaStressYZ,
                                        bondLevelPiolaStressZX,
                                        bondLevelPiolaStressZY,
                                        bondLevelPiolaStressZZ,
                                        stressIntegral,
                                        influenceState,
                                        neighborhoodList,
                                        numOwnedPoints);

  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double *w0 = weightedVolume;

  double* stressXX = bondLevelPiolaStressXX;
  double* stressXY = bondLevelPiolaStressXY;
  double* stressXZ = bondLevelPiolaStressXZ;
  double* stressYX = bondLevelPiolaStressYX;
  double* stressYY = bondLevelPiolaStressYY;
  double* stressYZ = bondLevelPiolaStressYZ;
  double* stressZX = bondLevelPiolaStressZX;
  double* stressZY = bondLevelPiolaStressZY;
  double* stressZZ = bondLevelPiolaStressZZ;

  double* stressInt = stressIntegral;

  double* phiX = gradientWeightX;
  double* phiY = gradientWeightY;
  double* phiZ = gradientWeightZ;

  double *omega = influenceState;

  double *modelCoordinatesPtr, *neighborModelmodelCoordinatesPtr, *forceDensityPtr, *neighborForceDensityPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLengthSq;
  double TX, TY, TZ, vol, neighborVol;
  int numNeighbors, neighborIndex;

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID, ++w0, stressInt+=9){
          
    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
            omega++, phiX++, phiY++, phiZ++, 
            stressXX++, stressXY++, stressXZ++, 
            stressYX++, stressYY++, stressYZ++, 
            stressZX++, stressZY++, stressZZ++){

      if(*omega > 0.0){
        neighborIndex = *neighborListPtr;
        neighborModelmodelCoordinatesPtr = modelCoordinates + 3*neighborIndex;

        undeformedBondX = *(neighborModelmodelCoordinatesPtr)   - *(modelCoordinatesPtr);
        undeformedBondY = *(neighborModelmodelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
        undeformedBondZ = *(neighborModelmodelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);
        undeformedBondLengthSq = undeformedBondX*undeformedBondX +
                                 undeformedBondY*undeformedBondY +
                                 undeformedBondZ*undeformedBondZ;

        TX = *omega / *w0 * ( *stressXX * undeformedBondX + *stressXY * undeformedBondY + *stressXZ * undeformedBondZ ) / undeformedBondLengthSq;
        TY = *omega / *w0 * ( *stressYX * undeformedBondX + *stressYY * undeformedBondY + *stressYZ * undeformedBondZ ) / undeformedBondLengthSq;
        TZ = *omega / *w0 * ( *stressZX * undeformedBondX + *stressZY * undeformedBondY + *stressZZ * undeformedBondZ ) / undeformedBondLengthSq;

        TX += ( *(stressInt+0) * *phiX + *(stressInt+1) * *phiY + *(stressInt+2) * *phiZ ); 
        TY += ( *(stressInt+3) * *phiX + *(stressInt+4) * *phiY + *(stressInt+5) * *phiZ ); 
        TZ += ( *(stressInt+6) * *phiX + *(stressInt+7) * *phiY + *(stressInt+8) * *phiZ ); 

        vol = volume[iID];
        neighborVol = volume[neighborIndex];

        forceDensityPtr = forceDensity + 3*iID;
        neighborForceDensityPtr = forceDensity + 3*neighborIndex;

        *(forceDensityPtr)   += TX * neighborVol;
        *(forceDensityPtr+1) += TY * neighborVol;
        *(forceDensityPtr+2) += TZ * neighborVol;
        *(neighborForceDensityPtr)   -= TX * vol;
        *(neighborForceDensityPtr+1) -= TY * vol;
        *(neighborForceDensityPtr+2) -= TZ * vol;
      }
    }
  }
}

void PeridigmNS::BondAssociatedCorrespondenceMaterial::precompute(const double dt,
                                                                  const int numOwnedPoints,
                                                                  const int* ownedIDs,
                                                                  const int* neighborhoodList,
                                                                  PeridigmNS::DataManager& dataManager) const
{
  // Zero out the data  
  dataManager.getData(m_deformationGradientXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientDotXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientDotYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientDotZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  double *horizon, *volume; 
  double *weightedVolume;
  double *modelCoordinates, *displacements, *velocities;
  double *gradientWeightX, *gradientWeightY, *gradientWeightZ;
  double *deformationGradientX, *deformationGradientY, *deformationGradientZ;
  double *deformationGradientDotX, *deformationGradientDotY, *deformationGradientDotZ;
  double *gradientWeightEvaluationFlag;
  double *damageN, *damageNP1, *bondDamage, *influenceState;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_displacementsFieldId, PeridigmField::STEP_NP1)->ExtractView(&displacements);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightX);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightY);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightZ);
  dataManager.getData(m_gradientWeightEvaluationFlagFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightEvaluationFlag);
  dataManager.getData(m_deformationGradientXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientX);
  dataManager.getData(m_deformationGradientYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientY);
  dataManager.getData(m_deformationGradientZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZ);
  dataManager.getData(m_deformationGradientDotXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientDotX);
  dataManager.getData(m_deformationGradientDotYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientDotY);
  dataManager.getData(m_deformationGradientDotZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientDotZ);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_N)->ExtractView(&damageN);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageNP1);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_influenceStateFieldId, PeridigmField::STEP_NONE)->ExtractView(&influenceState);

  // Identify points with growing damage
  CORRESPONDENCE::updateGradientWeightEvaluationFlag(damageN,
                                                     damageNP1,
                                                     gradientWeightEvaluationFlag,
                                                     numOwnedPoints);

  // The weights will be updated only when damage grows
  int gradientWeightReturnCode = CORRESPONDENCE::computeLagrangianGradientWeights(horizon,
                                                                                  modelCoordinates,
                                                                                  volume,
                                                                                  gradientWeightX,
                                                                                  gradientWeightY,
                                                                                  gradientWeightZ,
                                                                                  gradientWeightEvaluationFlag,
                                                                                  bondDamage,
                                                                                  influenceState,
                                                                                  m_accuracyOrder,
                                                                                  neighborhoodList,
                                                                                  numOwnedPoints);
          
  string gradientWeightErrorMessage =
    "**** Error:  BondAssociatedCorrespondenceMaterial::precompute() failed to compute gradient weights.\n";
  gradientWeightErrorMessage +=
    "****         Possible scenarios: 1) The horizon is too small, or 2) Too much damage around some points.\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(gradientWeightReturnCode != 0, gradientWeightErrorMessage);

  CORRESPONDENCE::computeDeformationGradient(volume,
                                             displacements,
                                             velocities,
                                             gradientWeightX,
                                             gradientWeightY,
                                             gradientWeightZ,
                                             deformationGradientX,
                                             deformationGradientY,
                                             deformationGradientZ,
                                             deformationGradientDotX,
                                             deformationGradientDotY,
                                             deformationGradientDotZ,
                                             neighborhoodList,
                                             numOwnedPoints);
}
