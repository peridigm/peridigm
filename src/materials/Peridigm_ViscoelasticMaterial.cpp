/*! \file Peridigm_ViscoelasticMaterial.cpp */

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

#include "Peridigm_ViscoelasticMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "viscoelastic.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <limits>

PeridigmNS::ViscoelasticMaterial::ViscoelasticMaterial(const Teuchos::ParameterList & params)
 : Material(params),
   m_applyAutomaticDifferentiationJacobian(false),
   volumeFieldId(-1), damageFieldId(-1), weightedVolumeFieldId(-1), dilatationFieldId(-1), modelCoordinatesFieldId(-1),
   coordinatesFieldId(-1), forceDensityFieldId(-1), bondDamageFieldId(-1), deviatoricBackExtensionFieldId(-1)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_shearModulus = params.get<double>("Shear Modulus");
  m_horizon = params.get<double>("Horizon");
  m_density = params.get<double>("Density");
  m_lambda_i = params.get<double>("lambda_i");
  m_tau_b = params.get<double>("tau b");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Automatic Differentiation Jacobian"), "**** Error:  Automatic Differentiation is not supported for the Viscoelastic material.\n");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor is not supported for the Viscoelastic material.\n");

  // set up vector of variable specs
  m_variableSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  m_variableSpecs->push_back(Field_NS::VOLUME);
  m_variableSpecs->push_back(Field_NS::DAMAGE);
  m_variableSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
  m_variableSpecs->push_back(Field_NS::DILATATION);
  m_variableSpecs->push_back(Field_NS::COORD3D);
  m_variableSpecs->push_back(Field_NS::CURCOORD3D);
  m_variableSpecs->push_back(Field_NS::FORCE_DENSITY3D);
  m_variableSpecs->push_back(Field_NS::BOND_DAMAGE);
  m_variableSpecs->push_back(Field_NS::DEVIATORIC_BACK_EXTENSION);


  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();

  volumeFieldId                      = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::CONSTANT, "Volume");
  damageFieldId                      = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  weightedVolumeFieldId              = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::CONSTANT, "Weighted_Volume");
  dilatationFieldId                  = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Dilatation");
  modelCoordinatesFieldId            = fieldManager.getFieldId(PeridigmNS::PeridigmField::NODE,    PeridigmNS::PeridigmField::VECTOR, PeridigmNS::PeridigmField::CONSTANT, "Model_Coordinates");
  coordinatesFieldId                 = fieldManager.getFieldId(PeridigmNS::PeridigmField::NODE,    PeridigmNS::PeridigmField::VECTOR, PeridigmNS::PeridigmField::TWO_STEP, "Coordinates");
  forceDensityFieldId                = fieldManager.getFieldId(PeridigmNS::PeridigmField::NODE,    PeridigmNS::PeridigmField::VECTOR, PeridigmNS::PeridigmField::TWO_STEP, "Force_Density");
  bondDamageFieldId                  = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND,    PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
  deviatoricBackExtensionFieldId     = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND,    PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Deviatoric_Back_Extension");
}

PeridigmNS::ViscoelasticMaterial::~ViscoelasticMaterial()
{
}

void PeridigmNS::ViscoelasticMaterial::initialize(const double dt,
                                                  const int numOwnedPoints,
                                                  const int* ownedIDs,
                                                  const int* neighborhoodList,
                                                  PeridigmNS::DataManager& dataManager) const
{
  double *xOverlap, *cellVolumeOverlap, *weightedVolume;
  // dataManager.getData(Field_NS::COORD3D, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  // dataManager.getData(Field_NS::VOLUME, PeridigmField::STEP_NONE)->ExtractView(&cellVolumeOverlap);
  // dataManager.getData(Field_NS::WEIGHTED_VOLUME, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolumeOverlap);
  dataManager.getData(weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);

  MATERIAL_EVALUATION::computeWeightedVolume(xOverlap,cellVolumeOverlap,weightedVolume,numOwnedPoints,neighborhoodList);
}

void
PeridigmNS::ViscoelasticMaterial::computeForce(const double dt,
                                                          const int numOwnedPoints,
                                                          const int* ownedIDs,
                                                          const int* neighborhoodList,
                                                          PeridigmNS::DataManager& dataManager) const
{
  double *x, *yN, *yNP1, *volume, *dilatationN, *dilatationNp1, *weightedVolume, *bondDamage, *edbN, *edbNP1,  *force;
  // dataManager.getData(Field_NS::COORD3D, PeridigmField::STEP_NONE)->ExtractView(&x);
  // dataManager.getData(Field_NS::VOLUME, PeridigmField::STEP_NONE)->ExtractView(&volume);
  // dataManager.getData(Field_NS::WEIGHTED_VOLUME, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  // dataManager.getData(Field_NS::CURCOORD3D, PeridigmField::STEP_N)->ExtractView(&yN);
  // dataManager.getData(Field_NS::CURCOORD3D, PeridigmField::STEP_NP1)->ExtractView(&yNP1);
  // dataManager.getData(Field_NS::DILATATION, PeridigmField::STEP_N)->ExtractView(&dilatationN);
  // dataManager.getData(Field_NS::DILATATION, PeridigmField::STEP_NP1)->ExtractView(&dilatationNp1);
  // dataManager.getData(Field_NS::DEVIATORIC_BACK_EXTENSION, PeridigmField::STEP_N)->ExtractView(&edbN);
  // dataManager.getData(Field_NS::DEVIATORIC_BACK_EXTENSION, PeridigmField::STEP_NP1)->ExtractView(&edbNP1);
  // dataManager.getData(Field_NS::BOND_DAMAGE, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  // dataManager.getData(Field_NS::FORCE_DENSITY3D, PeridigmField::STEP_NP1)->ExtractView(&force);

  // dataManager.getData(Field_NS::FORCE_DENSITY3D, PeridigmField::STEP_NP1)->PutScalar(0.0);

  dataManager.getData(modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(coordinatesFieldId, PeridigmField::STEP_N)->ExtractView(&yN);
  dataManager.getData(coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&yNP1);
  dataManager.getData(dilatationFieldId, PeridigmField::STEP_N)->ExtractView(&dilatationN);
  dataManager.getData(dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatationNp1);
  dataManager.getData(deviatoricBackExtensionFieldId, PeridigmField::STEP_N)->ExtractView(&edbN);
  dataManager.getData(deviatoricBackExtensionFieldId, PeridigmField::STEP_NP1)->ExtractView(&edbNP1);
  dataManager.getData(bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);

  dataManager.getData(forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  MATERIAL_EVALUATION::computeDilatation(x,yNP1,weightedVolume,volume,bondDamage,dilatationNp1,neighborhoodList,numOwnedPoints);
  MATERIAL_EVALUATION::computeInternalForceViscoelasticStandardLinearSolid(dt,
                                                                           x,
                                                                           yN,
                                                                           yNP1,
                                                                           weightedVolume,
                                                                           volume,
                                                                           dilatationN,
                                                                           dilatationNp1,
                                                                           bondDamage,
                                                                           edbN,
                                                                           edbNP1,
                                                                           force,
                                                                           neighborhoodList,
                                                                           numOwnedPoints,
                                                                           m_bulkModulus,
                                                                           m_shearModulus,
                                                                           m_lambda_i,
                                                                           m_tau_b);
}

