/*! \file Peridigm_LinearElasticIsotropicMaterial.cpp */

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

#include "Peridigm_LinearElasticIsotropicMaterial.hpp"
#include "Peridigm_CriticalStretchDamageModel.hpp"
#include <Teuchos_TestForException.hpp>
#include "ordinary_elastic.h"
#include "ordinary_utilities.h"

using namespace std;

PeridigmNS::LinearElasticIsotropicMaterial::LinearElasticIsotropicMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_damageModel()
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_shearModulus = params.get<double>("Shear Modulus");
  m_density = params.get<double>("Density");

  // set up vector of variable specs
  m_variableSpecs = Teuchos::rcp(new vector<Field_NS::FieldSpec>);
  m_variableSpecs->push_back(Field_NS::VOLUME);
  m_variableSpecs->push_back(Field_NS::DAMAGE);
  m_variableSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
  m_variableSpecs->push_back(Field_NS::DILATATION);
  m_variableSpecs->push_back(Field_NS::COORD3D);
  m_variableSpecs->push_back(Field_NS::CURCOORD3D);
  m_variableSpecs->push_back(Field_NS::FORCE_DENSITY3D);
  m_variableSpecs->push_back(Field_NS::BOND_DAMAGE);

  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > damageModelVariableSpecs;
  if(params.isSublist("Damage Model")){
    Teuchos::ParameterList damageParams = params.sublist("Damage Model");
    if(!damageParams.isParameter("Type")){
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         "Damage model \"Type\" not specified in Damage Model parameter list.");
    }
    string& damageModelType = damageParams.get<string>("Type");
    if(damageModelType == "Critical Stretch"){
      m_damageModel = Teuchos::rcp(new PeridigmNS::CriticalStretchDamageModel(damageParams));
      damageModelVariableSpecs = m_damageModel->VariableSpecs();
    }
    else{
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         "Invalid damage model, \"None\" or \"Critical Stretch\" required.");
    }

    // add damage model's variable specs to list of variable specs for this material
    // \todo Avoid duplicate specs here!
    for(unsigned int i=0 ; i<damageModelVariableSpecs->size() ; ++i)
      m_variableSpecs->push_back( (*damageModelVariableSpecs)[i] );
  }

}

PeridigmNS::LinearElasticIsotropicMaterial::~LinearElasticIsotropicMaterial()
{
}

void
PeridigmNS::LinearElasticIsotropicMaterial::initialize(const double dt,
                                                       const int numOwnedPoints,
                                                       const int* ownedIDs,
                                                       const int* neighborhoodList,
                                                       PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *cellVolume, *weightedVolume;
  dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&weightedVolume);

  MATERIAL_EVALUATION::computeWeightedVolume(x,cellVolume,weightedVolume,numOwnedPoints,neighborhoodList);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::updateConstitutiveData(const double dt,
                                                                   const int numOwnedPoints,
                                                                   const int* ownedIDs,
                                                                   const int* neighborhoodList,
                                                                   PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *y, *cellVolume, *weightedVolume, *dilatation, *damage, *bondDamage;
  dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&bondDamage);

  // Update the bond damgae
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

  MATERIAL_EVALUATION::computeDilatation(x,y,weightedVolume,cellVolume,bondDamage,dilatation,neighborhoodList,numOwnedPoints);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::computeForce(const double dt,
                                                         const int numOwnedPoints,
                                                         const int* ownedIDs,
                                                         const int* neighborhoodList,
                                                         PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *weightedVolume, *dilatation, *bondDamage, *force;
  dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&force);

  MATERIAL_EVALUATION::computeInternalForceLinearElastic(x,y,weightedVolume,cellVolume,dilatation,bondDamage,force,neighborhoodList,numOwnedPoints,m_bulkModulus,m_shearModulus);
}
