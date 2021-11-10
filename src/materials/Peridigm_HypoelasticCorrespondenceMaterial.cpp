/*! \file Peridigm_HypoelasticCorrespondenceMaterial.cpp */

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

#include "Peridigm_HypoelasticCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include <Teuchos_Assert.hpp>
#include <Sacado.hpp> // for MPI_abort

using namespace std;

PeridigmNS::HypoelasticCorrespondenceMaterial::HypoelasticCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0), m_actualHorizon(0.0),
    m_accuracyOrder(1), 
    m_OMEGA(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_horizonFieldId(-1), m_volumeFieldId(-1),
    m_coordinatesFieldId(-1), m_velocitiesFieldId(-1), 
    m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
    m_gradientWeightXFieldId(-1),
    m_gradientWeightYFieldId(-1),
    m_gradientWeightZFieldId(-1),
    m_velocityGradientFieldId(-1),
    m_velocityGradientXFieldId(-1),
    m_velocityGradientYFieldId(-1),
    m_velocityGradientZFieldId(-1),
    m_deformationGradientFieldId(-1),
    m_greenLagrangeStrainFieldId(-1),
    m_leftStretchTensorFieldId(-1),
    m_rotationTensorFieldId(-1), 
    m_unrotatedCauchyStressFieldId(-1),
    m_cauchyStressFieldId(-1), 
    m_unrotatedRateOfDeformationFieldId(-1),
    m_jacobianDeterminantFieldId(-1),
    m_weightedVolumeFieldId(-1),
    m_undamagedWeightedVolumeFieldId(-1),
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
    m_bondLevelVelocityGradientXXFieldId(-1),
    m_bondLevelVelocityGradientXYFieldId(-1),
    m_bondLevelVelocityGradientXZFieldId(-1),
    m_bondLevelVelocityGradientYXFieldId(-1),
    m_bondLevelVelocityGradientYYFieldId(-1),
    m_bondLevelVelocityGradientYZFieldId(-1),
    m_bondLevelVelocityGradientZXFieldId(-1),
    m_bondLevelVelocityGradientZYFieldId(-1),
    m_bondLevelVelocityGradientZZFieldId(-1),
    m_nonhomogeneousIntegralFieldId(-1),
    m_flyingPointFlagFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");

  m_actualHorizon = params.get<double>("Actual Horizon");

  if(params.isParameter("Gradient Order Of Accuracy")){
    m_accuracyOrder = params.get<int>("Gradient Order Of Accuracy");
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticHypoelasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticHypoelasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticHypoelasticCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_horizonFieldId                    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_volumeFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_coordinatesFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  m_forceDensityFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId                 = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_gradientWeightXFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_X");
  m_gradientWeightYFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_Y");
  m_gradientWeightZFieldId            = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Gradient_Weight_Z");
  m_velocityGradientFieldId           = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Velocity_Gradient");
  m_velocityGradientXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Velocity_Gradient_X");
  m_velocityGradientYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Velocity_Gradient_Y");
  m_velocityGradientZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Velocity_Gradient_Z");
  m_deformationGradientFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Deformation_Gradient");
  m_greenLagrangeStrainFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Green_Lagrange_Strain");
  m_leftStretchTensorFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Left_Stretch_Tensor");
  m_rotationTensorFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Rotation_Tensor");
  m_unrotatedCauchyStressFieldId      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Unrotated_Cauchy_Stress");
  m_cauchyStressFieldId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Cauchy_Stress");
  m_unrotatedRateOfDeformationFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_Deformation");
  m_jacobianDeterminantFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Jacobian_Determinant");
  m_weightedVolumeFieldId             = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
  m_undamagedWeightedVolumeFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Undamaged_Weighted_Volume");
  
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
  m_bondLevelVelocityGradientXXFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_XX");
  m_bondLevelVelocityGradientXYFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_XY");
  m_bondLevelVelocityGradientXZFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_XZ");
  m_bondLevelVelocityGradientYXFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_YX");
  m_bondLevelVelocityGradientYYFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_YY");
  m_bondLevelVelocityGradientYZFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_YZ");
  m_bondLevelVelocityGradientZXFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_ZX");
  m_bondLevelVelocityGradientZYFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_ZY");
  m_bondLevelVelocityGradientZZFieldId           = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Velocity_Gradient_ZZ");
  m_nonhomogeneousIntegralFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Nonhomogeneous_Integral");
  m_flyingPointFlagFieldId                       = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flying_Point_Flag");

  m_fieldIds.push_back(m_horizonFieldId);
  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_gradientWeightXFieldId);
  m_fieldIds.push_back(m_gradientWeightYFieldId);
  m_fieldIds.push_back(m_gradientWeightZFieldId);
  m_fieldIds.push_back(m_velocityGradientFieldId);
  m_fieldIds.push_back(m_velocityGradientXFieldId);
  m_fieldIds.push_back(m_velocityGradientYFieldId);
  m_fieldIds.push_back(m_velocityGradientZFieldId);
  m_fieldIds.push_back(m_deformationGradientFieldId);
  m_fieldIds.push_back(m_greenLagrangeStrainFieldId);
  m_fieldIds.push_back(m_leftStretchTensorFieldId);
  m_fieldIds.push_back(m_rotationTensorFieldId);
  m_fieldIds.push_back(m_unrotatedCauchyStressFieldId);
  m_fieldIds.push_back(m_cauchyStressFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationFieldId);
  m_fieldIds.push_back(m_jacobianDeterminantFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);
  m_fieldIds.push_back(m_undamagedWeightedVolumeFieldId);

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
  m_fieldIds.push_back(m_bondLevelVelocityGradientXXFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientXYFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientXZFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientYXFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientYYFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientYZFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientZXFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientZYFieldId);
  m_fieldIds.push_back(m_bondLevelVelocityGradientZZFieldId);
  m_fieldIds.push_back(m_nonhomogeneousIntegralFieldId);
  m_fieldIds.push_back(m_flyingPointFlagFieldId);
}

PeridigmNS::HypoelasticCorrespondenceMaterial::~HypoelasticCorrespondenceMaterial()
{
}

void
PeridigmNS::HypoelasticCorrespondenceMaterial::initialize(const double dt,
                                                          const int numOwnedPoints,
                                                          const int* ownedIDs,
                                                          const int* neighborhoodList,
                                                          PeridigmNS::DataManager& dataManager)
{
  // The horizon size needs to be corrected here. 
  // Horizon was specified before, only to form a bigger neighborhood lists.
  // Here an Eulerian frame is used, so we need to keep track of potential
  // points currently out of the scope of a material point but become neighbor
  // later.
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->PutScalar(m_actualHorizon);

  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_unrotatedRateOfDeformationFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_unrotatedCauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_cauchyStressFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *deformationGradientN;
  double *deformationGradientNP1;
  double *leftStretchTensorN;
  double *leftStretchTensorNP1;
  double *rotationTensorN;
  double *rotationTensorNP1;
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_N)->ExtractView(&deformationGradientN);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientNP1);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_N)->ExtractView(&leftStretchTensorN);
  dataManager.getData(m_leftStretchTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&leftStretchTensorNP1);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_N)->ExtractView(&rotationTensorN);
  dataManager.getData(m_rotationTensorFieldId, PeridigmField::STEP_NP1)->ExtractView(&rotationTensorNP1);

  //Initialize the left stretch and rotation tenor to the identity matrix
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(deformationGradientN, numOwnedPoints);
  CORRESPONDENCE::setOnesOnDiagonalFullTensor(deformationGradientNP1, numOwnedPoints);
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

  dataManager.getData(m_nonhomogeneousIntegralFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_undamagedWeightedVolumeFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_N)->PutScalar(1.0);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_NP1)->PutScalar(1.0);

  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_N)->PutScalar(-1.0);
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_NP1)->PutScalar(-1.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
}

void
PeridigmNS::HypoelasticCorrespondenceMaterial::computeForce(const double dt,
                                                            const int numOwnedPoints,
                                                            const int* ownedIDs,
                                                            const int* neighborhoodList,
                                                            PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces 
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *horizon, *volume, *coordinates, *velocities;
  double *weightedVolume, *jacobianDeterminantN, *jacobianDeterminantNP1;
  double *gradientWeightX, *gradientWeightY, *gradientWeightZ;
  double *velocityGradient, *velocityGradientX, *velocityGradientY, *velocityGradientZ;
  double *bondLevelVelocityGradientXX, *bondLevelVelocityGradientXY, *bondLevelVelocityGradientXZ;
  double *bondLevelVelocityGradientYX, *bondLevelVelocityGradientYY, *bondLevelVelocityGradientYZ;
  double *bondLevelVelocityGradientZX, *bondLevelVelocityGradientZY, *bondLevelVelocityGradientZZ;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_N)->ExtractView(&jacobianDeterminantN);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_NP1)->ExtractView(&jacobianDeterminantNP1);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightX);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightY);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightZ);
  dataManager.getData(m_velocityGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradient);
  dataManager.getData(m_velocityGradientXFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradientX);
  dataManager.getData(m_velocityGradientYFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradientY);
  dataManager.getData(m_velocityGradientZFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradientZ);
  dataManager.getData(m_bondLevelVelocityGradientXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientXX);
  dataManager.getData(m_bondLevelVelocityGradientXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientXY);
  dataManager.getData(m_bondLevelVelocityGradientXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientXZ);
  dataManager.getData(m_bondLevelVelocityGradientYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientYX);
  dataManager.getData(m_bondLevelVelocityGradientYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientYY);
  dataManager.getData(m_bondLevelVelocityGradientYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientYZ);
  dataManager.getData(m_bondLevelVelocityGradientZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientZX);
  dataManager.getData(m_bondLevelVelocityGradientZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientZY);
  dataManager.getData(m_bondLevelVelocityGradientZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&bondLevelVelocityGradientZZ);
  
  double *flyingPointFlag, *bondDamage;
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_N)->ExtractView(&flyingPointFlag);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamage);

  double *deformationGradientN, *deformationGradientNP1;
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_N)->ExtractView(&deformationGradientN);
  dataManager.getData(m_deformationGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientNP1);

  CORRESPONDENCE::computeBondLevelVelocityGradient(coordinates,
                                                   velocities,
                                                   velocityGradientX,
                                                   velocityGradientY,
                                                   velocityGradientZ,
                                                   bondLevelVelocityGradientXX,
                                                   bondLevelVelocityGradientXY,
                                                   bondLevelVelocityGradientXZ,
                                                   bondLevelVelocityGradientYX,
                                                   bondLevelVelocityGradientYY,
                                                   bondLevelVelocityGradientYZ,
                                                   bondLevelVelocityGradientZX,
                                                   bondLevelVelocityGradientZY,
                                                   bondLevelVelocityGradientZZ,
                                                   flyingPointFlag,
                                                   neighborhoodList,
                                                   numOwnedPoints);

  // Update the deformation gradient so that later Green Strain tensor can be
  // computed. This step is done only for visualization purposes.
  CORRESPONDENCE::updateDeformationGradient(velocityGradient,
                                            deformationGradientN,
                                            deformationGradientNP1,
                                            flyingPointFlag,
                                            numOwnedPoints,
                                            dt);

  double *greenLagrangeStrain;
  dataManager.getData(m_greenLagrangeStrainFieldId, PeridigmField::STEP_NONE)->ExtractView(&greenLagrangeStrain);

  // This step is done only for visualization purposes.
  CORRESPONDENCE::computeGreenLagrangeStrain(deformationGradientNP1,
                                             greenLagrangeStrain,
                                             flyingPointFlag,
                                             numOwnedPoints);

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
  int nodeLevelRotationTensorReturnCode = CORRESPONDENCE::computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor(velocityGradient,
                                                                                                                      leftStretchTensorN,
                                                                                                                      rotationTensorN,
                                                                                                                      leftStretchTensorNP1,
                                                                                                                      rotationTensorNP1,
                                                                                                                      unrotatedRateOfDeformation,
                                                                                                                      flyingPointFlag,
                                                                                                                      numOwnedPoints, 
                                                                                                                      dt);
  string nodeLevelRotationTensorErrorMessage =
    "**** Error:  HypoelasticCorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
  nodeLevelRotationTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";
  TEUCHOS_TEST_FOR_TERMINATION(nodeLevelRotationTensorReturnCode != 0, nodeLevelRotationTensorErrorMessage);

  // Compute bond-level values
  int bondLevelRotationTensorReturnCode = CORRESPONDENCE::computeBondLevelUnrotatedRateOfDeformationAndRotationTensor(bondLevelVelocityGradientXX,
                                                                                                                      bondLevelVelocityGradientXY,
                                                                                                                      bondLevelVelocityGradientXZ,
                                                                                                                      bondLevelVelocityGradientYX,
                                                                                                                      bondLevelVelocityGradientYY,
                                                                                                                      bondLevelVelocityGradientYZ,
                                                                                                                      bondLevelVelocityGradientZX,
                                                                                                                      bondLevelVelocityGradientZY,
                                                                                                                      bondLevelVelocityGradientZZ,
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
                                                                                                                      flyingPointFlag,
                                                                                                                      neighborhoodList, 
                                                                                                                      numOwnedPoints, 
                                                                                                                      dt);
  string bondLevelRotationTensorErrorMessage =
    "**** Error:  HypoelasticCorrespondenceMaterial::computeForce() failed to compute rotation tensor.\n";
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
                                     flyingPointFlag,
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
                                              flyingPointFlag,
                                              neighborhoodList,
                                              numOwnedPoints);

  // Cauchy stress is now updated and in the rotated state. Proceed with the
  // evaluation of the integral (related to the non-homogeneous part of the
  // deformation)

  double *nonhomogeneousIntegral;
  dataManager.getData(m_nonhomogeneousIntegralFieldId, PeridigmField::STEP_NONE)->ExtractView(&nonhomogeneousIntegral);

  double *undamagedWeightedVolume;
  dataManager.getData(m_undamagedWeightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&undamagedWeightedVolume);

  CORRESPONDENCE::computeNonhomogeneityIntegral(volume,
                                                undamagedWeightedVolume,
                                                jacobianDeterminantN,
                                                horizon,
                                                coordinates,
                                                bondLevelCauchyStressXXNP1,
                                                bondLevelCauchyStressXYNP1,
                                                bondLevelCauchyStressXZNP1,
                                                bondLevelCauchyStressYXNP1,
                                                bondLevelCauchyStressYYNP1,
                                                bondLevelCauchyStressYZNP1,
                                                bondLevelCauchyStressZXNP1,
                                                bondLevelCauchyStressZYNP1,
                                                bondLevelCauchyStressZZNP1,
                                                nonhomogeneousIntegral,
                                                flyingPointFlag,
                                                bondDamage,
                                                neighborhoodList,
                                                numOwnedPoints);

  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double *delta = horizon;
  double *w0 = undamagedWeightedVolume;
  double* stressXX = bondLevelCauchyStressXXNP1;
  double* stressXY = bondLevelCauchyStressXYNP1;
  double* stressXZ = bondLevelCauchyStressXZNP1;
  double* stressYX = bondLevelCauchyStressYXNP1;
  double* stressYY = bondLevelCauchyStressYYNP1;
  double* stressYZ = bondLevelCauchyStressYZNP1;
  double* stressZX = bondLevelCauchyStressZXNP1;
  double* stressZY = bondLevelCauchyStressZYNP1;
  double* stressZZ = bondLevelCauchyStressZZNP1;
  double* nonhomoIntegral = nonhomogeneousIntegral;

  double* phiX = gradientWeightX;
  double* phiY = gradientWeightY;
  double* phiZ = gradientWeightZ;

  double* flyingPointFlg = flyingPointFlag;
  double *bondDamagePtr = bondDamage;

  double *coordinatesPtr, *neighborCoordinatesPtr, *forceDensityPtr, *neighborForceDensityPtr;
  double deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  double TX, TY, TZ, omega, vol, neighborVol;
  int numNeighbors, neighborIndex;

  string matrixInversionErrorMessage =
    "**** Error:  HypoelasticCorrespondenceMaterial::computeForce() failed to invert deformation gradient.\n";
  matrixInversionErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID, ++delta, ++w0, 
        nonhomoIntegral+=9, flyingPointFlg++){
          
    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Loop over the neighbors and compute contribution to force densities
      coordinatesPtr = coordinates + 3*iID;
      numNeighbors = *neighborListPtr; neighborListPtr++;

      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++, 
              phiX++, phiY++, phiZ++, 
              stressXX++, stressXY++, stressXZ++, 
              stressYX++, stressYY++, stressYZ++, 
              stressZX++, stressZY++, stressZZ++){

        neighborIndex = *neighborListPtr;
        neighborCoordinatesPtr = coordinates + 3*neighborIndex;

        deformedBondX = *(neighborCoordinatesPtr)   - *(coordinatesPtr);
        deformedBondY = *(neighborCoordinatesPtr+1) - *(coordinatesPtr+1);
        deformedBondZ = *(neighborCoordinatesPtr+2) - *(coordinatesPtr+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);
        deformedBondLengthSq = deformedBondX*deformedBondX +
                               deformedBondY*deformedBondY +
                               deformedBondZ*deformedBondZ;

        omega = (1.0 - *bondDamagePtr) * m_OMEGA(deformedBondLength, *delta);

        if(omega > 0.0){
          TX = omega / *w0 * ( *stressXX * deformedBondX + *stressXY * deformedBondY + *stressXZ * deformedBondZ ) / deformedBondLengthSq;
          TY = omega / *w0 * ( *stressYX * deformedBondX + *stressYY * deformedBondY + *stressYZ * deformedBondZ ) / deformedBondLengthSq;
          TZ = omega / *w0 * ( *stressZX * deformedBondX + *stressZY * deformedBondY + *stressZZ * deformedBondZ ) / deformedBondLengthSq;

          TX += ( *(nonhomoIntegral+0) * *phiX + *(nonhomoIntegral+1) * *phiY + *(nonhomoIntegral+2) * *phiZ ); 
          TY += ( *(nonhomoIntegral+3) * *phiX + *(nonhomoIntegral+4) * *phiY + *(nonhomoIntegral+5) * *phiZ ); 
          TZ += ( *(nonhomoIntegral+6) * *phiX + *(nonhomoIntegral+7) * *phiY + *(nonhomoIntegral+8) * *phiZ ); 

          vol = jacobianDeterminantN[iID] * volume[iID];
          neighborVol = jacobianDeterminantN[neighborIndex] * volume[neighborIndex];

          forceDensityPtr = forceDensity + 3*iID;
          neighborForceDensityPtr = forceDensity + 3*neighborIndex;

          *(forceDensityPtr)   += jacobianDeterminantN[iID] * TX * neighborVol;
          *(forceDensityPtr+1) += jacobianDeterminantN[iID] * TY * neighborVol;
          *(forceDensityPtr+2) += jacobianDeterminantN[iID] * TZ * neighborVol;
          *(neighborForceDensityPtr)   -= jacobianDeterminantN[neighborIndex] * TX * vol;
          *(neighborForceDensityPtr+1) -= jacobianDeterminantN[neighborIndex] * TY * vol;
          *(neighborForceDensityPtr+2) -= jacobianDeterminantN[neighborIndex] * TZ * vol;
        }
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors;
      phiX += numNeighbors; phiY += numNeighbors; phiZ += numNeighbors; 
      stressXX += numNeighbors; stressXY += numNeighbors; stressXZ += numNeighbors; 
      stressYX += numNeighbors; stressYY += numNeighbors; stressYZ += numNeighbors; 
      stressZX += numNeighbors; stressZY += numNeighbors; stressZZ += numNeighbors;
    }
  }
}

void PeridigmNS::HypoelasticCorrespondenceMaterial::precompute(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* ownedIDs,
                                                               const int* neighborhoodList,
                                                               PeridigmNS::DataManager& dataManager) const
{
  // Zero out the data  
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_undamagedWeightedVolumeFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_velocityGradientZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);

  double *horizon, *volume; 
  double *weightedVolume, *undamagedWeightedVolume, *jacobianDeterminantN, *jacobianDeterminantNP1;
  double *coordinates, *velocities;
  double *gradientWeightX, *gradientWeightY, *gradientWeightZ;
  double *velocityGradient, *velocityGradientX, *velocityGradientY, *velocityGradientZ;
  double *flyingPointFlag, *bondDamage;
  dataManager.getData(m_horizonFieldId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_undamagedWeightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&undamagedWeightedVolume);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_N)->ExtractView(&jacobianDeterminantN);
  dataManager.getData(m_jacobianDeterminantFieldId, PeridigmField::STEP_NP1)->ExtractView(&jacobianDeterminantNP1);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);
  dataManager.getData(m_gradientWeightXFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightX);
  dataManager.getData(m_gradientWeightYFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightY);
  dataManager.getData(m_gradientWeightZFieldId, PeridigmField::STEP_NONE)->ExtractView(&gradientWeightZ);
  dataManager.getData(m_velocityGradientFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradient);
  dataManager.getData(m_velocityGradientXFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradientX);
  dataManager.getData(m_velocityGradientYFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradientY);
  dataManager.getData(m_velocityGradientZFieldId, PeridigmField::STEP_NONE)->ExtractView(&velocityGradientZ);
  dataManager.getData(m_flyingPointFlagFieldId, PeridigmField::STEP_N)->ExtractView(&flyingPointFlag);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamage);

  // Compute the current damaged weighted volume
  CORRESPONDENCE::computeWeightedVolume(volume,
                                        weightedVolume,
                                        jacobianDeterminantN,
                                        horizon,
                                        coordinates,
                                        flyingPointFlag,
                                        bondDamage,
                                        neighborhoodList,
                                        numOwnedPoints);

  int gradientWeightReturnCode = CORRESPONDENCE::computeGradientWeights(horizon,
                                                                        coordinates,
                                                                        volume,
                                                                        jacobianDeterminantN,
                                                                        gradientWeightX,
                                                                        gradientWeightY,
                                                                        gradientWeightZ,
                                                                        m_accuracyOrder,
                                                                        flyingPointFlag,
                                                                        bondDamage,
                                                                        neighborhoodList,
                                                                        numOwnedPoints);

  string gradientWeightErrorMessage =
    "**** Error:  HypoelasticCorrespondenceMaterial::precompute() failed to compute gradient weights.\n";
  gradientWeightErrorMessage +=
    "****         Possible scenarios: 1) The horizon is too small, or 2) Too much damage around some points.\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(gradientWeightReturnCode != 0, gradientWeightErrorMessage);

  CORRESPONDENCE::computeVelocityGradient(volume,
                                          jacobianDeterminantN,
                                          jacobianDeterminantNP1,
                                          velocities,
                                          gradientWeightX,
                                          gradientWeightY,
                                          gradientWeightZ,
                                          velocityGradient,
                                          velocityGradientX,
                                          velocityGradientY,
                                          velocityGradientZ,
                                          flyingPointFlag,
                                          neighborhoodList,
                                          numOwnedPoints,
                                          dt);

  // Compute the current undamaged weighted volume
  CORRESPONDENCE::computeUndamagedWeightedVolume(volume,
                                                 undamagedWeightedVolume,
                                                 jacobianDeterminantN,
                                                 horizon,
                                                 coordinates,
                                                 neighborhoodList,
                                                 numOwnedPoints);
}
