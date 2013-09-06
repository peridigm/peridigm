/*! \file Peridigm_ElasticCorrespondenceMaterial.cpp */

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

#include "Peridigm_ElasticCorrespondenceMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::ElasticCorrespondenceMaterial::ElasticCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_youngsModulus(0.0), m_poissonsRatio(0.0), m_density(0.0), m_horizon(0.0), m_hourglassCoefficient(0.0), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_hourglassForceDensityFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
    m_deformationGradientXXFieldId(-1), m_deformationGradientXYFieldId(-1), m_deformationGradientXZFieldId(-1), 
    m_deformationGradientYXFieldId(-1), m_deformationGradientYYFieldId(-1), m_deformationGradientYZFieldId(-1), 
    m_deformationGradientZXFieldId(-1), m_deformationGradientZYFieldId(-1), m_deformationGradientZZFieldId(-1),
    m_strainXXFieldId(-1), m_strainXYFieldId(-1), m_strainXZFieldId(-1), 
    m_strainYXFieldId(-1), m_strainYYFieldId(-1), m_strainYZFieldId(-1), 
    m_strainZXFieldId(-1), m_strainZYFieldId(-1), m_strainZZFieldId(-1),
    m_stressXXFieldId(-1), m_stressXYFieldId(-1), m_stressXZFieldId(-1), 
    m_stressYXFieldId(-1), m_stressYYFieldId(-1), m_stressYZFieldId(-1), 
    m_stressZXFieldId(-1), m_stressZYFieldId(-1), m_stressZZFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_youngsModulus = (9.0*m_bulkModulus*m_shearModulus)/(3.0*m_bulkModulus + m_shearModulus);
  m_poissonsRatio = (3.0*m_bulkModulus - 2.0*m_shearModulus)/(6.0*m_bulkModulus + 2.0*m_shearModulus);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  m_hourglassCoefficient = params.get<double>("Hourglass Coefficient");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
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
  m_shapeTensorInverseXXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXX");
  m_shapeTensorInverseXYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXY");
  m_shapeTensorInverseXZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXZ");
  m_shapeTensorInverseYXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYX");
  m_shapeTensorInverseYYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYY");
  m_shapeTensorInverseYZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYZ");
  m_shapeTensorInverseZXFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZX");
  m_shapeTensorInverseZYFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZY");
  m_shapeTensorInverseZZFieldId    = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZZ");
  m_strainXXFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainXX");
  m_strainXYFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainXY");
  m_strainXZFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainXZ");
  m_strainYXFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainYX");
  m_strainYYFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainYY");
  m_strainYZFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainYZ");
  m_strainZXFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainZX");
  m_strainZYFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainZY");
  m_strainZZFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StrainZZ");
  m_stressXXFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressXX");
  m_stressXYFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressXY");
  m_stressXZFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressXZ");
  m_stressYXFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressYX");
  m_stressYYFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressYY");
  m_stressYZFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressYZ");
  m_stressZXFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressZX");
  m_stressZYFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressZY");
  m_stressZZFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "StressZZ");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
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
  m_fieldIds.push_back(m_shapeTensorInverseXXFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseXYFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseXZFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseYXFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseYYFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseYZFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseZXFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseZYFieldId);
  m_fieldIds.push_back(m_shapeTensorInverseZZFieldId);
  m_fieldIds.push_back(m_strainXXFieldId);
  m_fieldIds.push_back(m_strainXYFieldId);
  m_fieldIds.push_back(m_strainXZFieldId);
  m_fieldIds.push_back(m_strainYXFieldId);
  m_fieldIds.push_back(m_strainYYFieldId);
  m_fieldIds.push_back(m_strainYZFieldId);
  m_fieldIds.push_back(m_strainZXFieldId);
  m_fieldIds.push_back(m_strainZYFieldId);
  m_fieldIds.push_back(m_strainZZFieldId);
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

PeridigmNS::ElasticCorrespondenceMaterial::~ElasticCorrespondenceMaterial()
{
}

void
PeridigmNS::ElasticCorrespondenceMaterial::computeForce(const double dt,
                                                        const int numOwnedPoints,
                                                        const int* ownedIDs,
                                                        const int* neighborhoodList,
                                                        PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // \todo Optimize this function and move it to a cxx file; get rid of base-class implementation
  // \todo Separate def grad calculation from shape tensor inverse calculation and move to material_utilities.cxx (create one with damage and one without damage, the one with damage needs heursitics to keep running, check for possible use in Sierra/SM)

  computeApproximateDeformationGradient(numOwnedPoints, ownedIDs, neighborhoodList, dataManager);

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

  double *strainXX, *strainXY, *strainXZ;
  double *strainYX, *strainYY, *strainYZ;
  double *strainZX, *strainZY, *strainZZ;
  dataManager.getData(m_strainXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainXX);
  dataManager.getData(m_strainXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainXY);
  dataManager.getData(m_strainXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainXZ);
  dataManager.getData(m_strainYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainYX);
  dataManager.getData(m_strainYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainYY);
  dataManager.getData(m_strainYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainYZ);
  dataManager.getData(m_strainZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainZX);
  dataManager.getData(m_strainZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainZY);
  dataManager.getData(m_strainZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&strainZZ);

  // Green-Lagrange Strain E = 0.5*(F^T F - I)
  CORRESPONDENCE::computeGreenLagrangeStrain(deformationGradientXX, deformationGradientXY, deformationGradientXZ,
                                             deformationGradientYX, deformationGradientYY, deformationGradientYZ,
                                             deformationGradientZX, deformationGradientZY, deformationGradientZZ,
                                             strainXX, strainXY, strainXZ,
                                             strainYX, strainYY, strainYZ,
                                             strainZX, strainZY, strainZZ,
                                             numOwnedPoints);

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

  double *cauchyStressXX, *cauchyStressXY, *cauchyStressXZ;
  double *cauchyStressYX, *cauchyStressYY, *cauchyStressYZ;
  double *cauchyStressZX, *cauchyStressZY, *cauchyStressZZ;
  dataManager.getData(m_stressXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressXX);
  dataManager.getData(m_stressXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressXY);
  dataManager.getData(m_stressXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressXZ);
  dataManager.getData(m_stressYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressYX);
  dataManager.getData(m_stressYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressYY);
  dataManager.getData(m_stressYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressYZ);
  dataManager.getData(m_stressZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressZX);
  dataManager.getData(m_stressZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressZY);
  dataManager.getData(m_stressZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&cauchyStressZZ);

  double *volume, *modelCoordinates, *forceDensity;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  CORRESPONDENCE::computeClassicalElasticStress(strainXX, strainXY, strainXZ,
                                                strainYX, strainYY, strainYZ,
                                                strainZX, strainZY, strainZZ,
                                                cauchyStressXX, cauchyStressXY, cauchyStressXZ,
                                                cauchyStressYX, cauchyStressYY, cauchyStressYZ,
                                                cauchyStressZX, cauchyStressZY, cauchyStressZZ,
                                                numOwnedPoints,
                                                m_youngsModulus,
                                                m_poissonsRatio);

  double* defGradXX = deformationGradientXX;
  double* defGradXY = deformationGradientXY;
  double* defGradXZ = deformationGradientXZ;
  double* defGradYX = deformationGradientYX;
  double* defGradYY = deformationGradientYY;
  double* defGradYZ = deformationGradientYZ;
  double* defGradZX = deformationGradientZX;
  double* defGradZY = deformationGradientZY;
  double* defGradZZ = deformationGradientZZ;

  double* stressXX = cauchyStressXX;
  double* stressXY = cauchyStressXY;
  double* stressXZ = cauchyStressXZ;
  double* stressYX = cauchyStressYX;
  double* stressYY = cauchyStressYY;
  double* stressYZ = cauchyStressYZ;
  double* stressZX = cauchyStressZX;
  double* stressZY = cauchyStressZY;
  double* stressZZ = cauchyStressZZ;

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

  double *modelCoordinatesPtr, *neighborModelCoordinatesPtr, *forceDensityPtr, *neighborForceDensityPtr;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double TX, TY, TZ;
  double omega, vol, neighborVol;
  int numNeighbors, neighborIndex;

  string matrixInversionErrorMessage =
    "**** Error:  ElasticCorrespondenceMaterial::computeForce() failed to invert deformation gradient.\n";
  matrixInversionErrorMessage +=
    "****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n";

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

    TEUCHOS_TEST_FOR_EXCEPT_MSG(matrixInversionReturnCode != 0, matrixInversionErrorMessage;)

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

  double *coordinates, *hourglassForceDensity;
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
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
