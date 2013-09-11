/*! \file Peridigm_CorrespondenceMaterial2.cpp */

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

#include "Peridigm_CorrespondenceMaterial2.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::CorrespondenceMaterial2::CorrespondenceMaterial2(const Teuchos::ParameterList& params)
  : Material(params),
    m_density(0.0), m_horizon(0.0), m_hourglassCoefficient(0.0), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_velocitiesFieldId(-1), 
    m_hourglassForceDensityFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
    m_deformationGradientXXFieldId(-1), m_deformationGradientXYFieldId(-1), m_deformationGradientXZFieldId(-1), 
    m_deformationGradientYXFieldId(-1), m_deformationGradientYYFieldId(-1), m_deformationGradientYZFieldId(-1), 
    m_deformationGradientZXFieldId(-1), m_deformationGradientZYFieldId(-1), m_deformationGradientZZFieldId(-1),
    m_shapeTensorInverseXXFieldId(-1), m_shapeTensorInverseXYFieldId(-1), m_shapeTensorInverseXZFieldId(-1), 
    m_shapeTensorInverseYXFieldId(-1), m_shapeTensorInverseYYFieldId(-1), m_shapeTensorInverseYZFieldId(-1), 
    m_shapeTensorInverseZXFieldId(-1), m_shapeTensorInverseZYFieldId(-1), m_shapeTensorInverseZZFieldId(-1),
    m_leftStretchTensorXXFieldId(-1), m_leftStretchTensorXYFieldId(-1), m_leftStretchTensorXZFieldId(-1), 
    m_leftStretchTensorYXFieldId(-1), m_leftStretchTensorYYFieldId(-1), m_leftStretchTensorYZFieldId(-1), 
    m_leftStretchTensorZXFieldId(-1), m_leftStretchTensorZYFieldId(-1), m_leftStretchTensorZZFieldId(-1),
    m_rotationTensorXXFieldId(-1), m_rotationTensorXYFieldId(-1), m_rotationTensorXZFieldId(-1), 
    m_rotationTensorYXFieldId(-1), m_rotationTensorYYFieldId(-1), m_rotationTensorYZFieldId(-1), 
    m_rotationTensorZXFieldId(-1), m_rotationTensorZYFieldId(-1), m_rotationTensorZZFieldId(-1),
    m_cauchyStressXXFieldId(-1), m_cauchyStressXYFieldId(-1), m_cauchyStressXZFieldId(-1), 
    m_cauchyStressYXFieldId(-1), m_cauchyStressYYFieldId(-1), m_cauchyStressYZFieldId(-1), 
    m_cauchyStressZXFieldId(-1), m_cauchyStressZYFieldId(-1), m_cauchyStressZZFieldId(-1),
    m_unrotatedRateOfDeformationXXFieldId(-1), m_unrotatedRateOfDeformationXYFieldId(-1), m_unrotatedRateOfDeformationXZFieldId(-1), 
    m_unrotatedRateOfDeformationYXFieldId(-1), m_unrotatedRateOfDeformationYYFieldId(-1), m_unrotatedRateOfDeformationYZFieldId(-1), 
    m_unrotatedRateOfDeformationZXFieldId(-1), m_unrotatedRateOfDeformationZYFieldId(-1), m_unrotatedRateOfDeformationZZFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence2 material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence2 material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticCorrespondence2 material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_velocitiesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP,  "Velocity");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_hourglassForceDensityFieldId   = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Hourglass_Force_Density");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_deformationGradientXXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXX");
  m_deformationGradientXYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXY");
  m_deformationGradientXZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXZ");
  m_deformationGradientYXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYX");
  m_deformationGradientYYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYY");
  m_deformationGradientYZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYZ");
  m_deformationGradientZXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZX");
  m_deformationGradientZYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZY");
  m_deformationGradientZZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZZ");
  m_leftStretchTensorXXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorXX");
  m_leftStretchTensorXYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorXY");
  m_leftStretchTensorXZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorXZ");
  m_leftStretchTensorYXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorYX");
  m_leftStretchTensorYYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorYY");
  m_leftStretchTensorYZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorYZ");
  m_leftStretchTensorZXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorZX");
  m_leftStretchTensorZYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorZY");
  m_leftStretchTensorZZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Left_Stretch_TensorZZ");
  m_rotationTensorXXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorXX");
  m_rotationTensorXYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorXY");
  m_rotationTensorXZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorXZ");
  m_rotationTensorYXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorYX");
  m_rotationTensorYYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorYY");
  m_rotationTensorYZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorYZ");
  m_rotationTensorZXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorZX");
  m_rotationTensorZYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorZY");
  m_rotationTensorZZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Rotation_TensorZZ");
  m_shapeTensorInverseXXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXX");
  m_shapeTensorInverseXYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXY");
  m_shapeTensorInverseXZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXZ");
  m_shapeTensorInverseYXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYX");
  m_shapeTensorInverseYYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYY");
  m_shapeTensorInverseYZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYZ");
  m_shapeTensorInverseZXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZX");
  m_shapeTensorInverseZYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZY");
  m_shapeTensorInverseZZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZZ");
  m_cauchyStressXXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressXX");
  m_cauchyStressXYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressXY");
  m_cauchyStressXZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressXZ");
  m_cauchyStressYXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressYX");
  m_cauchyStressYYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressYY");
  m_cauchyStressYZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressYZ");
  m_cauchyStressZXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressZX");
  m_cauchyStressZYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressZY");
  m_cauchyStressZZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Cauchy_StressZZ");
  m_unrotatedRateOfDeformationXXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationXX");
  m_unrotatedRateOfDeformationXYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationXY");
  m_unrotatedRateOfDeformationXZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationXZ");
  m_unrotatedRateOfDeformationYXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationYX");
  m_unrotatedRateOfDeformationYYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationYY");
  m_unrotatedRateOfDeformationYZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationYZ");
  m_unrotatedRateOfDeformationZXFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationZX");
  m_unrotatedRateOfDeformationZYFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationZY");
  m_unrotatedRateOfDeformationZZFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Unrotated_Rate_Of_DeformationZZ");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocitiesFieldId);
  m_fieldIds.push_back(m_hourglassForceDensityFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_deformationGradientXXFieldId);
  m_fieldIds.push_back(m_deformationGradientXYFieldId);
  m_fieldIds.push_back(m_deformationGradientXZFieldId);
  m_fieldIds.push_back(m_deformationGradientYXFieldId);
  m_fieldIds.push_back(m_deformationGradientYYFieldId);
  m_fieldIds.push_back(m_deformationGradientYZFieldId);
  m_fieldIds.push_back(m_deformationGradientZXFieldId);
  m_fieldIds.push_back(m_deformationGradientZYFieldId);
  m_fieldIds.push_back(m_deformationGradientZZFieldId);
  m_fieldIds.push_back(m_leftStretchTensorXXFieldId);
  m_fieldIds.push_back(m_leftStretchTensorXYFieldId);
  m_fieldIds.push_back(m_leftStretchTensorXZFieldId);
  m_fieldIds.push_back(m_leftStretchTensorYXFieldId);
  m_fieldIds.push_back(m_leftStretchTensorYYFieldId);
  m_fieldIds.push_back(m_leftStretchTensorYZFieldId);
  m_fieldIds.push_back(m_leftStretchTensorZXFieldId);
  m_fieldIds.push_back(m_leftStretchTensorZYFieldId);
  m_fieldIds.push_back(m_leftStretchTensorZZFieldId);
  m_fieldIds.push_back(m_rotationTensorXXFieldId);
  m_fieldIds.push_back(m_rotationTensorXYFieldId);
  m_fieldIds.push_back(m_rotationTensorXZFieldId);
  m_fieldIds.push_back(m_rotationTensorYXFieldId);
  m_fieldIds.push_back(m_rotationTensorYYFieldId);
  m_fieldIds.push_back(m_rotationTensorYZFieldId);
  m_fieldIds.push_back(m_rotationTensorZXFieldId);
  m_fieldIds.push_back(m_rotationTensorZYFieldId);
  m_fieldIds.push_back(m_rotationTensorZZFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseXXFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseXYFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseXZFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseYXFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseYYFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseYZFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseZXFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseZYFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseZZFieldId);
  m_fieldIds.push_back(m_cauchyStressXXFieldId);
  m_fieldIds.push_back(m_cauchyStressXYFieldId);
  m_fieldIds.push_back(m_cauchyStressXZFieldId);
  m_fieldIds.push_back(m_cauchyStressYXFieldId);
  m_fieldIds.push_back(m_cauchyStressYYFieldId);
  m_fieldIds.push_back(m_cauchyStressYZFieldId);
  m_fieldIds.push_back(m_cauchyStressZXFieldId);
  m_fieldIds.push_back(m_cauchyStressZYFieldId);
  m_fieldIds.push_back(m_cauchyStressZZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationXZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationYZFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZXFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZYFieldId);
  m_fieldIds.push_back(m_unrotatedRateOfDeformationZZFieldId);
}

PeridigmNS::CorrespondenceMaterial2::~CorrespondenceMaterial2()
{
}

void
PeridigmNS::CorrespondenceMaterial2::initialize(const double dt,
                                               const int numOwnedPoints,
                                               const int* ownedIDs,
                                               const int* neighborhoodList,
                                               PeridigmNS::DataManager& dataManager) const
{
  dataManager.getData(m_deformationGradientXXFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_deformationGradientXYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientXZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientYXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientYYFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_deformationGradientYZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientZXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientZYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_deformationGradientZZFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);

  dataManager.getData(m_leftStretchTensorXYFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_leftStretchTensorXYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorXZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorYXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorYYFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_leftStretchTensorYZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorZXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorZYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_leftStretchTensorZZFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  
  dataManager.getData(m_rotationTensorXXFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_rotationTensorXYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorXZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorYXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorYYFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
  dataManager.getData(m_rotationTensorYZFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorZXFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorZYFieldId, PeridigmField::STEP_NONE)->PutScalar(0.0);
  dataManager.getData(m_rotationTensorZZFieldId, PeridigmField::STEP_NONE)->PutScalar(1.0);
}

void
PeridigmNS::CorrespondenceMaterial2::computeForce(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *volume, *modelCoordinates, *coordinates, *velocities;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
  dataManager.getData(m_velocitiesFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocities);

  double *shapeTensorInverseXX, *shapeTensorInverseXY, *shapeTensorInverseXZ;
  double *shapeTensorInverseYX, *shapeTensorInverseYY, *shapeTensorInverseYZ;
  double *shapeTensorInverseZX, *shapeTensorInverseZY, *shapeTensorInverseZZ;
  dataManager.getData(m_shapeTensorInverseXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseXX);
  dataManager.getData(m_shapeTensorInverseXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseXY);
  dataManager.getData(m_shapeTensorInverseXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseXZ);
  dataManager.getData(m_shapeTensorInverseYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseYX);
  dataManager.getData(m_shapeTensorInverseYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseYY);
  dataManager.getData(m_shapeTensorInverseYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseYZ);
  dataManager.getData(m_shapeTensorInverseZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseZX);
  dataManager.getData(m_shapeTensorInverseZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseZY);
  dataManager.getData(m_shapeTensorInverseZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseZZ);

  double *deformationGradientXX, *deformationGradientXY, *deformationGradientXZ;
  double *deformationGradientYX, *deformationGradientYY, *deformationGradientYZ;
  double *deformationGradientZX, *deformationGradientZY, *deformationGradientZZ;
  dataManager.getData(m_deformationGradientXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientXX);
  dataManager.getData(m_deformationGradientXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientXY);
  dataManager.getData(m_deformationGradientXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientXZ);
  dataManager.getData(m_deformationGradientYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientYX);
  dataManager.getData(m_deformationGradientYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientYY);
  dataManager.getData(m_deformationGradientYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientYZ);
  dataManager.getData(m_deformationGradientZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZX);
  dataManager.getData(m_deformationGradientZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZY);
  dataManager.getData(m_deformationGradientZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZZ);

  // Compute the inverse of the shape tensor and the approximate deformation gradient
  // The approximate deformation gradient will be used by the derived class (specific correspondence material model)
  // to compute the Cauchy stress.
  // The inverse of the shape tensor is stored for later use after the Cauchy stress calculation

  int shapeTensorReturnCode = 
    CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                               modelCoordinates,
                                                                               coordinates,
                                                                               shapeTensorInverseXX, shapeTensorInverseXY, shapeTensorInverseXZ,
                                                                               shapeTensorInverseYX, shapeTensorInverseYY, shapeTensorInverseYZ,
                                                                               shapeTensorInverseZX, shapeTensorInverseZY, shapeTensorInverseZZ,
                                                                               deformationGradientXX, deformationGradientXY, deformationGradientXZ,
                                                                               deformationGradientYX, deformationGradientYY, deformationGradientYZ,
                                                                               deformationGradientZX, deformationGradientZY, deformationGradientZZ,
                                                                               neighborhoodList,
                                                                               numOwnedPoints,
                                                                               m_horizon);

  string shapeTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial2::computeForce() failed to compute shape tensor.\n";
  shapeTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";

  TEUCHOS_TEST_FOR_EXCEPT_MSG(shapeTensorReturnCode != 0, shapeTensorErrorMessage);

  double *leftStretchTensorXX, *leftStretchTensorXY, *leftStretchTensorXZ;
  double *leftStretchTensorYX, *leftStretchTensorYY, *leftStretchTensorYZ;
  double *leftStretchTensorZX, *leftStretchTensorZY, *leftStretchTensorZZ;
  dataManager.getData(m_leftStretchTensorXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorXX);
  dataManager.getData(m_leftStretchTensorXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorXY);
  dataManager.getData(m_leftStretchTensorXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorXZ);
  dataManager.getData(m_leftStretchTensorYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorYX);
  dataManager.getData(m_leftStretchTensorYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorYY);
  dataManager.getData(m_leftStretchTensorYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorYZ);
  dataManager.getData(m_leftStretchTensorZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorZX);
  dataManager.getData(m_leftStretchTensorZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorZY);
  dataManager.getData(m_leftStretchTensorZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&leftStretchTensorZZ);

  double *rotationTensorXX, *rotationTensorXY, *rotationTensorXZ;
  double *rotationTensorYX, *rotationTensorYY, *rotationTensorYZ;
  double *rotationTensorZX, *rotationTensorZY, *rotationTensorZZ;
  dataManager.getData(m_rotationTensorXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorXX);
  dataManager.getData(m_rotationTensorXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorXY);
  dataManager.getData(m_rotationTensorXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorXZ);
  dataManager.getData(m_rotationTensorYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorYX);
  dataManager.getData(m_rotationTensorYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorYY);
  dataManager.getData(m_rotationTensorYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorYZ);
  dataManager.getData(m_rotationTensorZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorZX);
  dataManager.getData(m_rotationTensorZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorZY);
  dataManager.getData(m_rotationTensorZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&rotationTensorZZ);

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

  // Compute left stretch tensor and rotation tensor, returns these updated
  // quantities along with the unrotated rate-of-deformation.  Performs a polar
  // decomposition via Flanagan & Taylor (1987) algorithm.
  int rotationTensorReturnCode = 
      CORRESPONDENCE::computeUnrotatedRateOfDeformationAndRotationTensor(volume,
                                                                         modelCoordinates, 
                                                                         coordinates, 
                                                                         velocities, 
                                                                         deformationGradientXX, deformationGradientXY, deformationGradientXZ,
                                                                         deformationGradientYX, deformationGradientYY, deformationGradientYZ,
                                                                         deformationGradientZX, deformationGradientZY, deformationGradientZZ,  
                                                                         shapeTensorInverseXX, shapeTensorInverseXY, shapeTensorInverseXZ,
                                                                         shapeTensorInverseYX, shapeTensorInverseYY, shapeTensorInverseYZ, 
                                                                         shapeTensorInverseZX, shapeTensorInverseZY, shapeTensorInverseZZ,
                                                                         leftStretchTensorXX, leftStretchTensorXY, leftStretchTensorXZ, 
                                                                         leftStretchTensorYX, leftStretchTensorYY, leftStretchTensorYZ, 
                                                                         leftStretchTensorZX, leftStretchTensorZY, leftStretchTensorZZ, 
                                                                         rotationTensorXX, rotationTensorXY, rotationTensorXZ, 
                                                                         rotationTensorYX, rotationTensorYY, rotationTensorYZ, 
                                                                         rotationTensorZX, rotationTensorZY, rotationTensorZZ, 
                                                                         unrotatedRateOfDeformationXX, unrotatedRateOfDeformationXY, unrotatedRateOfDeformationXZ, 
                                                                         unrotatedRateOfDeformationYX, unrotatedRateOfDeformationYY, unrotatedRateOfDeformationYZ, 
                                                                         unrotatedRateOfDeformationZX, unrotatedRateOfDeformationZY, unrotatedRateOfDeformationZZ,
                                                                         neighborhoodList, numOwnedPoints, m_horizon, dt);
  
  string rotationTensorErrorMessage =
    "**** Error:  CorrespondenceMaterial2::computeForce() failed to compute rotation tensor.\n";
  rotationTensorErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";

  TEUCHOS_TEST_FOR_EXCEPT_MSG(rotationTensorReturnCode != 0, shapeTensorErrorMessage);

  // Evaluate the Cauchy stress using the routine implemented in the derived class (specific correspondence material model)
  // The general idea is to compute the stress based on:
  //   1) The unrotated rate-of-deformation tensor
  //   2) The time step
  //   3) Whatever state variables are managed by the derived class
  //
  double *cauchyStressXX, *cauchyStressXY, *cauchyStressXZ;
  double *cauchyStressYX, *cauchyStressYY, *cauchyStressYZ;
  double *cauchyStressZX, *cauchyStressZY, *cauchyStressZZ;
  dataManager.getData(m_cauchyStressXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressXX);
  dataManager.getData(m_cauchyStressXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressXY);
  dataManager.getData(m_cauchyStressXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressXZ);
  dataManager.getData(m_cauchyStressYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressYX);
  dataManager.getData(m_cauchyStressYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressYY);
  dataManager.getData(m_cauchyStressYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressYZ);
  dataManager.getData(m_cauchyStressZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressZX);
  dataManager.getData(m_cauchyStressZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressZY);
  dataManager.getData(m_cauchyStressZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressZZ);

  // unrotate the Cauchy stress into the corotational frame
  CORRESPONDENCE::unrotateCauchyStress(rotationTensorXX, rotationTensorXY, rotationTensorXZ, 
                                       rotationTensorYX, rotationTensorYY, rotationTensorYZ, 
                                       rotationTensorZX, rotationTensorZY, rotationTensorZZ,
                                       cauchyStressXX, cauchyStressXY, cauchyStressXZ,
                                       cauchyStressYX, cauchyStressYY, cauchyStressYZ,
                                       cauchyStressZX, cauchyStressZY, cauchyStressZZ,
                                       numOwnedPoints);

  // updateCauchyStress implemented in the derived class.  
  //   Input: unrotated rate-of-deformation tensor
  //   Input: unrotated Cauchy stress at step N
  //   Output unrotated Cauchy stress at step N+1
  //
  //   Internal state variables are managed in the derived class
  computeCauchyStress(dt, numOwnedPoints, dataManager);

  // rotate back to the Eulerian frame
  CORRESPONDENCE::rotateCauchyStress(rotationTensorXX, rotationTensorXY, rotationTensorXZ, 
                                     rotationTensorYX, rotationTensorYY, rotationTensorYZ, 
                                     rotationTensorZX, rotationTensorZY, rotationTensorZZ,
                                     cauchyStressXX, cauchyStressXY, cauchyStressXZ,
                                     cauchyStressYX, cauchyStressYY, cauchyStressYZ,
                                     cauchyStressZX, cauchyStressZY, cauchyStressZZ,
                                     numOwnedPoints);

  //Cauchy stress is now updated and in the rotated state.  Proceed with
  //conversion to Piola-Kirchoff and force-vector states.
  double* stressXX = cauchyStressXX;
  double* stressXY = cauchyStressXY;
  double* stressXZ = cauchyStressXZ;
  double* stressYX = cauchyStressYX;
  double* stressYY = cauchyStressYY;
  double* stressYZ = cauchyStressYZ;
  double* stressZX = cauchyStressZX;
  double* stressZY = cauchyStressZY;
  double* stressZZ = cauchyStressZZ;

  double* defGradXX = deformationGradientXX;
  double* defGradXY = deformationGradientXY;
  double* defGradXZ = deformationGradientXZ;
  double* defGradYX = deformationGradientYX;
  double* defGradYY = deformationGradientYY;
  double* defGradYZ = deformationGradientYZ;
  double* defGradZX = deformationGradientZX;
  double* defGradZY = deformationGradientZY;
  double* defGradZZ = deformationGradientZZ;

  double* shapeTensorInvXX = shapeTensorInverseXX;
  double* shapeTensorInvXY = shapeTensorInverseXY;
  double* shapeTensorInvXZ = shapeTensorInverseXZ;
  double* shapeTensorInvYX = shapeTensorInverseYX;
  double* shapeTensorInvYY = shapeTensorInverseYY;
  double* shapeTensorInvYZ = shapeTensorInverseYZ;
  double* shapeTensorInvZX = shapeTensorInverseZX;
  double* shapeTensorInvZY = shapeTensorInverseZY;
  double* shapeTensorInvZZ = shapeTensorInverseZZ;

  double defGradInvXX, defGradInvXY, defGradInvXZ;
  double defGradInvYX, defGradInvYY, defGradInvYZ;
  double defGradInvZX, defGradInvZY, defGradInvZZ;

  double piolaStressXX, piolaStressXY, piolaStressXZ;
  double piolaStressYX, piolaStressYY, piolaStressYZ;
  double piolaStressZX, piolaStressZY, piolaStressZZ;

  double tempXX, tempXY, tempXZ;
  double tempYX, tempYY, tempYZ;
  double tempZX, tempZY, tempZZ;

  double *forceDensity;
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double *modelCoordinatesPtr, *neighborModelCoordinatesPtr, *forceDensityPtr, *neighborForceDensityPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double TX, TY, TZ;
  double omega, vol, neighborVol;
  int numNeighbors, neighborIndex;

  string matrixInversionErrorMessage =
    "**** Error:  CorrespondenceMaterial2::computeForce() failed to invert deformation gradient.\n";
  matrixInversionErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";

  // Loop over the material points and convert the Cauchy stress into pairwise peridynamic force densities

  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID,
        ++defGradXX, ++defGradXY, ++defGradXZ,
        ++defGradYX, ++defGradYY, ++defGradYZ,
        ++defGradZX, ++defGradZY, ++defGradZZ,
        ++stressXX, ++stressXY, ++stressXZ,
        ++stressYX, ++stressYY, ++stressYZ,
        ++stressZX, ++stressZY, ++stressZZ,
        ++shapeTensorInvXX, ++shapeTensorInvXY, ++shapeTensorInvXZ,
        ++shapeTensorInvYX, ++shapeTensorInvYY, ++shapeTensorInvYZ,
        ++shapeTensorInvZX, ++shapeTensorInvZY, ++shapeTensorInvZZ){

    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient
    int matrixInversionReturnCode =
      CORRESPONDENCE::invert3by3Matrix(*defGradXX, *defGradXY, *defGradXZ,
                                       *defGradYX, *defGradYY, *defGradYZ,
                                       *defGradZX, *defGradZY, *defGradZZ,
                                       defGradInvXX, defGradInvXY, defGradInvXZ,
                                       defGradInvYX, defGradInvYY, defGradInvYZ,
                                       defGradInvZX, defGradInvZY, defGradInvZZ);

    TEUCHOS_TEST_FOR_EXCEPT_MSG(matrixInversionReturnCode != 0, matrixInversionErrorMessage);

    // Compute the 1st Piola stress as the inner product of the Cauchy stress 
    // and the tranpose of the inverse of the deformation gradient
    CORRESPONDENCE::MatrixMultiply(*stressXX, *stressXY, *stressXZ,
                                   *stressYX, *stressYY, *stressYZ,
                                   *stressZX, *stressZY, *stressZZ,
                                   defGradInvXX, defGradInvYX, defGradInvZX, /* note: TRANSPOSE of the inverse def grad */
                                   defGradInvXY, defGradInvYY, defGradInvZY,
                                   defGradInvXZ, defGradInvYZ, defGradInvZZ,
                                   piolaStressXX, piolaStressXY, piolaStressXZ,
                                   piolaStressYX, piolaStressYY, piolaStressYZ,
                                   piolaStressZX, piolaStressZY, piolaStressZZ);

    // Inner product of Piola stress and the inverse of the shape tensor
    CORRESPONDENCE::MatrixMultiply(piolaStressXX, piolaStressXY, piolaStressXZ,
                                   piolaStressYX, piolaStressYY, piolaStressYZ,
                                   piolaStressZX, piolaStressZY, piolaStressZZ,
                                   *shapeTensorInvXX, *shapeTensorInvXY, *shapeTensorInvXZ,
                                   *shapeTensorInvYX, *shapeTensorInvYY, *shapeTensorInvYZ,
                                   *shapeTensorInvZX, *shapeTensorInvZY, *shapeTensorInvZZ,
                                   tempXX, tempXY, tempXZ,
                                   tempYX, tempYY, tempYZ,
                                   tempZX, tempZY, tempZZ);

    // Loop over the neighbors and compute contribution to force densities
    modelCoordinatesPtr = modelCoordinates + 3*iID;
    numNeighbors = *neighborListPtr; neighborListPtr++;

    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborModelCoordinatesPtr = modelCoordinates + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoordinatesPtr) - *(modelCoordinatesPtr);
      undeformedBondY = *(neighborModelCoordinatesPtr+1) - *(modelCoordinatesPtr+1);
      undeformedBondZ = *(neighborModelCoordinatesPtr+2) - *(modelCoordinatesPtr+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, m_horizon);

      TX = omega * (tempXX*undeformedBondX + tempXY*undeformedBondY + tempXZ*undeformedBondZ);
      TY = omega * (tempYX*undeformedBondX + tempYY*undeformedBondY + tempYZ*undeformedBondZ);
      TZ = omega * (tempZX*undeformedBondX + tempZY*undeformedBondY + tempZZ*undeformedBondZ);

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

  // Compute hourglass forces for stabilization of low-enery and/or zero-energy modes
  dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *hourglassForceDensity;
  dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&hourglassForceDensity);

  CORRESPONDENCE::computeHourglassForce(volume,
                                        modelCoordinates,
                                        coordinates,
                                        deformationGradientXX, deformationGradientXY, deformationGradientXZ,
                                        deformationGradientYX, deformationGradientYY, deformationGradientYZ,
                                        deformationGradientZX, deformationGradientZY, deformationGradientZZ,
                                        hourglassForceDensity,
                                        neighborhoodList,
                                        numOwnedPoints,
                                        m_horizon,
                                        m_bulkModulus,
                                        m_hourglassCoefficient);

  // Sum the hourglass force densities into the force densities
  Teuchos::RCP<Epetra_Vector> forceDensityVector = dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1);
  Teuchos::RCP<Epetra_Vector> hourglassForceDensityVector = dataManager.getData(m_hourglassForceDensityFieldId, PeridigmField::STEP_NP1);
  forceDensityVector->Update(1.0, *hourglassForceDensityVector, 1.0);
}
