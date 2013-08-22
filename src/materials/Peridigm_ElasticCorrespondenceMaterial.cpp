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
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::ElasticCorrespondenceMaterial::ElasticCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_youngsModulus(0.0), m_poissonsRatio(0.0), m_density(0.0), m_horizon(0.0), m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
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

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the ElasticCorrespondence material model.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Thermal Expansion Coefficient"), "**** Error:  Thermal expansion is not currently supported for the ElasticCorrespondence material model.\n");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
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

// void
// PeridigmNS::ElasticCorrespondenceMaterial::initialize(const double dt,
//                                                       const int numOwnedPoints,
//                                                       const int* ownedIDs,
//                                                       const int* neighborhoodList,
//                                                       PeridigmNS::DataManager& dataManager) const
// {
//   PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
//   if(fieldManager.hasFieldId("Bond_Damage"){
//       int bondDamageFieId = fieldManager.getFieldId("Bond_Damage");
//       fieldManager.getData(bondDamangeFieldId, PeridigmField::STEP_N)->PutScalar(0.0);
//       fieldManager.getData(bondDamangeFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);
//     }
// }

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
  // \todo Create stand-alone matrix inversion routine.
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
  int nodeId(0);
  double defGrad[3][3], greenLagrangeStrain[3][3];
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    nodeId = ownedIDs[iID];

    defGrad[0][0] = deformationGradientXX[nodeId] ; defGrad[0][1] = deformationGradientXY[nodeId] ; defGrad[0][2] = deformationGradientXZ[nodeId] ; 
    defGrad[1][0] = deformationGradientYX[nodeId] ; defGrad[1][1] = deformationGradientYY[nodeId] ; defGrad[1][2] = deformationGradientYZ[nodeId] ; 
    defGrad[2][0] = deformationGradientZX[nodeId] ; defGrad[2][1] = deformationGradientZY[nodeId] ; defGrad[2][2] = deformationGradientZZ[nodeId] ; 
    for(int i=0 ; i<3 ; ++i){
      for(int j=0 ; j<3 ; ++j){
        greenLagrangeStrain[i][j] = 0.5 * ( defGrad[0][i]*defGrad[0][j] + defGrad[1][i]*defGrad[1][j] + defGrad[2][i]*defGrad[2][j] );
        if(i == j)
          greenLagrangeStrain[i][j] -= 0.5;
      }
    }
    strainXX[nodeId] = greenLagrangeStrain[0][0] ; strainXY[nodeId] = greenLagrangeStrain[0][1] ; strainXZ[nodeId] = greenLagrangeStrain[0][2] ; 
    strainYX[nodeId] = greenLagrangeStrain[1][0] ; strainYY[nodeId] = greenLagrangeStrain[1][1] ; strainYZ[nodeId] = greenLagrangeStrain[1][2] ; 
    strainZX[nodeId] = greenLagrangeStrain[2][0] ; strainZY[nodeId] = greenLagrangeStrain[2][1] ; strainZZ[nodeId] = greenLagrangeStrain[2][2] ; 
  }

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

  double *stressXX, *stressXY, *stressXZ;
  double *stressYX, *stressYY, *stressYZ;
  double *stressZX, *stressZY, *stressZZ;
  dataManager.getData(m_stressXXFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressXX);
  dataManager.getData(m_stressXYFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressXY);
  dataManager.getData(m_stressXZFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressXZ);
  dataManager.getData(m_stressYXFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressYX);
  dataManager.getData(m_stressYYFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressYY);
  dataManager.getData(m_stressYZFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressYZ);
  dataManager.getData(m_stressZXFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressZX);
  dataManager.getData(m_stressZYFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressZY);
  dataManager.getData(m_stressZZFieldId, PeridigmField::STEP_NONE)->ExtractView(&stressZZ);

  double *volume, *modelCoordinates, *forceDensity;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&forceDensity);

  double cauchyStress[3][3], piolaStress[3][3], defGradInverse[3][3], minor[9], det, temp[3][3];
  const int *neighborListPtr = neighborhoodList;

  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    nodeId = ownedIDs[iID];

    // Hooke's Law    
    // \todo Take advantange of symmetry, don't store full tensor
    double constant = m_youngsModulus/((1.0 + m_poissonsRatio)*(1.0 - 2.0*m_poissonsRatio));
    cauchyStress[0][0] = constant*( (1.0 - m_poissonsRatio)*strainXX[nodeId] +         m_poissonsRatio*strainYY[nodeId] +         m_poissonsRatio*strainZZ[nodeId] ) ;
    cauchyStress[1][1] = constant*(        m_poissonsRatio*strainXX[nodeId]  + (1.0 - m_poissonsRatio)*strainYY[nodeId] +         m_poissonsRatio*strainZZ[nodeId] ) ;
    cauchyStress[2][2] = constant*(        m_poissonsRatio*strainXX[nodeId]  +         m_poissonsRatio*strainYY[nodeId] + (1.0 - m_poissonsRatio)*strainZZ[nodeId] ) ;
    cauchyStress[0][1] = constant*(1.0 - 2.0*m_poissonsRatio)*strainXY[nodeId];
    cauchyStress[1][2] = constant*(1.0 - 2.0*m_poissonsRatio)*strainYZ[nodeId];
    cauchyStress[2][0] = constant*(1.0 - 2.0*m_poissonsRatio)*strainZX[nodeId];
    cauchyStress[1][0] = cauchyStress[0][1];
    cauchyStress[2][1] = cauchyStress[1][2];
    cauchyStress[0][2] = cauchyStress[2][0];

    // Store the Cauchy stress for output
    stressXX[nodeId] = cauchyStress[0][0] ; stressXY[nodeId] = cauchyStress[0][1] ; stressXZ[nodeId] = cauchyStress[0][2] ; 
    stressYX[nodeId] = cauchyStress[1][0] ; stressYY[nodeId] = cauchyStress[1][1] ; stressYZ[nodeId] = cauchyStress[1][2] ; 
    stressZX[nodeId] = cauchyStress[2][0] ; stressZY[nodeId] = cauchyStress[2][1] ; stressZZ[nodeId] = cauchyStress[2][2] ; 

    // first Piola-Kirchhoff stress = J * cauchyStress * defGrad^-T

    // Invert the deformation gradient

    defGrad[0][0] = deformationGradientXX[nodeId] ; defGrad[0][1] = deformationGradientXY[nodeId] ; defGrad[0][2] = deformationGradientXZ[nodeId] ; 
    defGrad[1][0] = deformationGradientYX[nodeId] ; defGrad[1][1] = deformationGradientYY[nodeId] ; defGrad[1][2] = deformationGradientYZ[nodeId] ; 
    defGrad[2][0] = deformationGradientZX[nodeId] ; defGrad[2][1] = deformationGradientZY[nodeId] ; defGrad[2][2] = deformationGradientZZ[nodeId] ; 
    minor[0] = defGrad[1][1]*defGrad[2][2] - defGrad[1][2]*defGrad[2][1] ;
    minor[1] = defGrad[1][0]*defGrad[2][2] - defGrad[1][2]*defGrad[2][0] ;
    minor[2] = defGrad[1][0]*defGrad[2][1] - defGrad[1][1]*defGrad[2][0] ;
    minor[3] = defGrad[0][1]*defGrad[2][2] - defGrad[0][2]*defGrad[2][1] ;
    minor[4] = defGrad[0][0]*defGrad[2][2] - defGrad[2][0]*defGrad[0][2] ;
    minor[5] = defGrad[0][0]*defGrad[2][1] - defGrad[0][1]*defGrad[2][0] ;
    minor[6] = defGrad[0][1]*defGrad[1][2] - defGrad[0][2]*defGrad[1][1] ;
    minor[7] = defGrad[0][0]*defGrad[1][2] - defGrad[0][2]*defGrad[1][0] ;
    minor[8] = defGrad[0][0]*defGrad[1][1] - defGrad[0][1]*defGrad[1][0] ;
    det = defGrad[0][0]*minor[0] - defGrad[0][1]*minor[1] + defGrad[0][2]*minor[2] ;

    TEUCHOS_TEST_FOR_EXCEPT_MSG(det == 0.0, "**** Error:  ElasticCorrespondenceMaterial::computeForce() failed to invert deformation gradient.\n****         Note that all nodes must have a minimum of three neighbors.  Is the horizon too small?\n");

    defGradInverse[0][0] = minor[0]/det;
    defGradInverse[0][1] = -1.0*minor[3]/det;
    defGradInverse[0][2] = minor[6]/det;
    defGradInverse[1][0] = -1.0*minor[1]/det;
    defGradInverse[1][1] = minor[4]/det;
    defGradInverse[1][2] = -1.0*minor[7]/det;
    defGradInverse[2][0] = minor[2]/det;
    defGradInverse[2][1] = -1.0*minor[5]/det;
    defGradInverse[2][2] = minor[8]/det;

    for(int i=0 ; i<3 ; ++i)
      for(int j=0 ; j<3 ; ++j)
        piolaStress[i][j] = ( cauchyStress[i][0]*defGradInverse[j][0] + cauchyStress[i][1]*defGradInverse[j][1] +  cauchyStress[i][2]*defGradInverse[j][2] ) * det;

    // Inner product of Piola stress and the inverse of the shape tensor
    temp[0][0] =  piolaStress[0][0]*shapeTensorInverseXX[iID] + piolaStress[0][1]*shapeTensorInverseYX[iID] + piolaStress[0][2]*shapeTensorInverseZX[iID];
    temp[0][1] =  piolaStress[0][0]*shapeTensorInverseXY[iID] + piolaStress[0][1]*shapeTensorInverseYY[iID] + piolaStress[0][2]*shapeTensorInverseZY[iID];
    temp[0][2] =  piolaStress[0][0]*shapeTensorInverseXZ[iID] + piolaStress[0][1]*shapeTensorInverseYZ[iID] + piolaStress[0][2]*shapeTensorInverseZZ[iID];
    temp[1][0] =  piolaStress[1][0]*shapeTensorInverseXX[iID] + piolaStress[1][1]*shapeTensorInverseYX[iID] + piolaStress[1][2]*shapeTensorInverseZX[iID];
    temp[1][1] =  piolaStress[1][0]*shapeTensorInverseXY[iID] + piolaStress[1][1]*shapeTensorInverseYY[iID] + piolaStress[1][2]*shapeTensorInverseZY[iID];
    temp[1][2] =  piolaStress[1][0]*shapeTensorInverseXZ[iID] + piolaStress[1][1]*shapeTensorInverseYZ[iID] + piolaStress[1][2]*shapeTensorInverseZZ[iID];
    temp[2][0] =  piolaStress[2][0]*shapeTensorInverseXX[iID] + piolaStress[2][1]*shapeTensorInverseYX[iID] + piolaStress[2][2]*shapeTensorInverseZX[iID];
    temp[2][1] =  piolaStress[2][0]*shapeTensorInverseXY[iID] + piolaStress[2][1]*shapeTensorInverseYY[iID] + piolaStress[2][2]*shapeTensorInverseZY[iID];
    temp[2][2] =  piolaStress[2][0]*shapeTensorInverseXZ[iID] + piolaStress[2][1]*shapeTensorInverseYZ[iID] + piolaStress[2][2]*shapeTensorInverseZZ[iID];

    // Loop over the neighbors and compute contribution to force densities
    double undeformedBond[3], undeformedBondLength, omega, vol, neighborVol, T[3];
    int neighborIndex;
    int numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){
      neighborIndex = *neighborListPtr;
      undeformedBond[0] = modelCoordinates[3*neighborIndex]   - modelCoordinates[3*iID];
      undeformedBond[1] = modelCoordinates[3*neighborIndex+1] - modelCoordinates[3*iID+1];
      undeformedBond[2] = modelCoordinates[3*neighborIndex+2] - modelCoordinates[3*iID+2];
      undeformedBondLength = sqrt(undeformedBond[0]*undeformedBond[0] +
                                  undeformedBond[1]*undeformedBond[1] +
                                  undeformedBond[2]*undeformedBond[2]);
      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, m_horizon);

      T[0] = omega * (temp[0][0]*undeformedBond[0] + temp[0][1]*undeformedBond[1] + temp[0][2]*undeformedBond[2]);
      T[1] = omega * (temp[1][0]*undeformedBond[0] + temp[1][1]*undeformedBond[1] + temp[1][2]*undeformedBond[2]);
      T[2] = omega * (temp[2][0]*undeformedBond[0] + temp[2][1]*undeformedBond[1] + temp[2][2]*undeformedBond[2]);

      vol = volume[iID];
      neighborVol = volume[neighborIndex];

      forceDensity[3*iID]   += T[0] * neighborVol;
      forceDensity[3*iID+1] += T[1] * neighborVol;
      forceDensity[3*iID+2] += T[2] * neighborVol;
      forceDensity[3*neighborIndex]   -= T[0] * vol;
      forceDensity[3*neighborIndex+1] -= T[1] * vol;
      forceDensity[3*neighborIndex+2] -= T[2] * vol;
    }
  }
}
