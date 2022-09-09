/*! \file Peridigm_ElasticBondBasedMaterial.cpp */

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

#include "Peridigm_ElasticBondBasedMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_bond_based.h"
#include <Teuchos_Assert.hpp>

PeridigmNS::ElasticBondBasedMaterial::ElasticBondBasedMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_bulkModulus(0.0), m_shearModulus(0.0), m_density(0.0), m_horizon(0.0), m_isPlaneStrain(false), m_isPlaneStress(false), m_isPlaneStrainStress(false),
    m_height(0.0), m_volumeFieldId(-1), m_damageFieldId(-1), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1)
{
  Teuchos::ParameterList internalParams = Teuchos::ParameterList(params); // Copy list such that it can be modified internally
  //! \todo Add meaningful asserts on material properties.
  m_density = internalParams.get<double>("Density");
  m_horizon = internalParams.get<double>("Horizon");
  if(internalParams.isParameter("Poisson's Ratio") || internalParams.isParameter("Shear Modulus") || (internalParams.isParameter("Bulk Modulus") && internalParams.isParameter("Young's Modulus"))){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  The Elastic bond based material model supports only one elastic constant, the bulk modulus or Young's modulus.");
  }
  m_isPlaneStrain = internalParams.get<bool>("Plane Strain", false);
  m_isPlaneStress = internalParams.get<bool>("Plane Stress", false);
  if(m_isPlaneStrain && m_isPlaneStress)
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error: 'Plane Strain' and 'Plane Stress' are mutual exclusive!");
  if(m_isPlaneStrain || m_isPlaneStress){
    m_isPlaneStrainStress = true;
    if(!internalParams.isParameter("Height"))
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error: In case of 'Plane Strain' and 'Plane Stress' the height must be given!")
  }
  m_height = internalParams.get<double>("Height", 0.0);

  internalParams.set<double>("Poisson's Ratio", 0.25); // 3d and plane strain case
  if(m_isPlaneStress)
    internalParams.set<double>("Poisson's Ratio", 1./3.); // 2d plane stress
  m_shearModulus = calculateShearModulus(internalParams);
  m_bulkModulus  = calculateBulkModulus(internalParams);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Damage");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Bond_Damage");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
}

PeridigmNS::ElasticBondBasedMaterial::~ElasticBondBasedMaterial()
{
}

void
PeridigmNS::ElasticBondBasedMaterial::initialize(const double dt,
                                                 const int numOwnedPoints,
                                                 const int* ownedIDs,
                                                 const int* neighborhoodList,
                                                 PeridigmNS::DataManager& dataManager)
{
}

void
PeridigmNS::ElasticBondBasedMaterial::computeForce(const double dt,
                                          const int numOwnedPoints,
                                          const int* ownedIDs,
                                          const int* neighborhoodList,
                                          PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *bondDamage, *force;

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);

  MATERIAL_EVALUATION::computeInternalForceElasticBondBased(x,y,cellVolume,bondDamage,force,neighborhoodList,numOwnedPoints,m_bulkModulus,m_horizon,m_isPlaneStrainStress,m_height);
}
