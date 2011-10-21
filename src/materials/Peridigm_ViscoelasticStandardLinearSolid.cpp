/*! \file Peridigm_ViscoelasticStandardLinearSolid.cpp */

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

#include "Peridigm_ViscoelasticStandardLinearSolid.hpp"
#include "Peridigm_DamageModelFactory.hpp"
#include <Teuchos_TestForException.hpp>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include "ordinary_std_linear_visco_solid.h"
#include "ordinary_utilities.h"
#include <limits>

PeridigmNS::ViscoelasticStandardLinearSolid::ViscoelasticStandardLinearSolid(const Teuchos::ParameterList & params)
:
Material(params),
m_damageModel()
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_shearModulus = params.get<double>("Shear Modulus");
  m_horizon = params.get<double>("Material Horizon");
  m_density = params.get<double>("Density");
  m_lambda_i   = params.get<double>("lambda_i");
  m_tau_b = params.get<double>("tau b");

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
  
  // Create the damage model, if any
  if(params.isSublist("Damage Model")){
    DamageModelFactory damageModelFactory;
    m_damageModel = damageModelFactory.create( params.sublist("Damage Model") );
    // Add damage model's variable specs to list of material model's variable specs
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > damageModelSpecs = m_damageModel->VariableSpecs();
    for(unsigned int i=0 ; i<damageModelSpecs->size() ; ++i)
      m_variableSpecs->push_back( (*damageModelSpecs)[i] );
  }
}

PeridigmNS::ViscoelasticStandardLinearSolid::~ViscoelasticStandardLinearSolid()
{
}

void PeridigmNS::ViscoelasticStandardLinearSolid::initialize(const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data
  double *xOverlap, *cellVolumeOverlap, *weightedVolume;
  dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolumeOverlap);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);

  MATERIAL_EVALUATION::computeWeightedVolume(xOverlap,cellVolumeOverlap,weightedVolume,numOwnedPoints,neighborhoodList);
}

void
PeridigmNS::ViscoelasticStandardLinearSolid::updateConstitutiveData(const double dt,
                                                                    const int numOwnedPoints,
                                                                    const int* ownedIDs,
                                                                    const int* neighborhoodList,
                                                                    PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *yNP1, *volume, *dilatation, *damage, *weightedVolume, *bondDamage;
  dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&yNP1);
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);

  // Update the bond damage
  if(!m_damageModel.is_null()){
    m_damageModel->computeDamage(dt,
                                 numOwnedPoints,
                                 ownedIDs,
                                 neighborhoodList,
                                 dataManager);
  }

  //  Update the element damage (percent of bonds broken)
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    int nodeID = ownedIDs[iID];
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
    double totalDamage = 0.0;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      totalDamage += bondDamage[bondIndex++];
    }
    if(numNeighbors > 0)
      totalDamage /= numNeighbors;
    else
      totalDamage = 0.0;
    damage[nodeID] = totalDamage;
  }

  MATERIAL_EVALUATION::computeDilatation(x,yNP1,weightedVolume,volume,bondDamage,dilatation,neighborhoodList,numOwnedPoints);
}

void
PeridigmNS::ViscoelasticStandardLinearSolid::computeForce(const double dt,
                                                          const int numOwnedPoints,
                                                          const int* ownedIDs,
                                                          const int* neighborhoodList,
                                                          PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *yN, *yNP1, *volume, *dilatationN, *dilatationNp1, *weightedVolume, *bondDamage, *edbN, *edbNP1,  *force;
  dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_N)->ExtractView(&yN);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&yNP1);
  dataManager.getData(Field_NS::DILATATION, Field_ENUM::STEP_N)->ExtractView(&dilatationN);
  dataManager.getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1)->ExtractView(&dilatationNp1);
  dataManager.getData(Field_NS::DEVIATORIC_BACK_EXTENSION, Field_ENUM::STEP_N)->ExtractView(&edbN);
  dataManager.getData(Field_NS::DEVIATORIC_BACK_EXTENSION, Field_ENUM::STEP_NP1)->ExtractView(&edbNP1);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->ExtractView(&force);

  // Zero out the force
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->PutScalar(0.0);

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

