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
#include "correspondence.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>

using namespace std;

PeridigmNS::ElasticCorrespondenceMaterial::ElasticCorrespondenceMaterial(const Teuchos::ParameterList& params)
  : CorrespondenceMaterial(params),
    m_youngsModulus(0.0), m_poissonsRatio(0.0),
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
  // m_bulkModulus and m_shearModulus are set within the CorrespondenceMaterial base class,
  // use them to compute Young's modulus and Poisson's ratio
  m_youngsModulus = (9.0*m_bulkModulus*m_shearModulus)/(3.0*m_bulkModulus + m_shearModulus);
  m_poissonsRatio = (3.0*m_bulkModulus - 2.0*m_shearModulus)/(6.0*m_bulkModulus + 2.0*m_shearModulus);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_deformationGradientXXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXX");
  m_deformationGradientXYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXY");
  m_deformationGradientXZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXZ");
  m_deformationGradientYXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYX");
  m_deformationGradientYYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYY");
  m_deformationGradientYZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYZ");
  m_deformationGradientZXFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZX");
  m_deformationGradientZYFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZY");
  m_deformationGradientZZFieldId   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZZ");
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

  m_fieldIds.push_back(m_deformationGradientXXFieldId);
  m_fieldIds.push_back(m_deformationGradientXYFieldId);
  m_fieldIds.push_back(m_deformationGradientXZFieldId);
  m_fieldIds.push_back(m_deformationGradientYXFieldId);
  m_fieldIds.push_back(m_deformationGradientYYFieldId);
  m_fieldIds.push_back(m_deformationGradientYZFieldId);
  m_fieldIds.push_back(m_deformationGradientZXFieldId);
  m_fieldIds.push_back(m_deformationGradientZYFieldId);
  m_fieldIds.push_back(m_deformationGradientZZFieldId);
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
PeridigmNS::ElasticCorrespondenceMaterial::computeCauchyStress(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* ownedIDs,
                                                               const int* neighborhoodList,
                                                               PeridigmNS::DataManager& dataManager) const
{
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

  CORRESPONDENCE::computeClassicalElasticStress(strainXX, strainXY, strainXZ,
                                                strainYX, strainYY, strainYZ,
                                                strainZX, strainZY, strainZZ,
                                                cauchyStressXX, cauchyStressXY, cauchyStressXZ,
                                                cauchyStressYX, cauchyStressYY, cauchyStressYZ,
                                                cauchyStressZX, cauchyStressZY, cauchyStressZZ,
                                                numOwnedPoints,
                                                m_youngsModulus,
                                                m_poissonsRatio);
}
