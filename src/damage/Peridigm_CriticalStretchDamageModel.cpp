/*! \file Peridigm_CriticalStretchDamageModel.cpp */

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

#include "Peridigm_CriticalStretchDamageModel.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

PeridigmNS::CriticalStretchDamageModel::CriticalStretchDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_damageFieldId(-1), m_bondDamageFieldId(-1)
{
  m_criticalStretch = params.get<double>("Critical Stretch");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");

  // set up vector of variable specs
  m_variableSpecs = Teuchos::rcp(new vector<Field_NS::FieldSpec>);
  m_variableSpecs->push_back(Field_NS::DAMAGE);
  m_variableSpecs->push_back(Field_NS::BOND_DAMAGE);
}

PeridigmNS::CriticalStretchDamageModel::~CriticalStretchDamageModel()
{
}

void
PeridigmNS::CriticalStretchDamageModel::initialize(const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   PeridigmNS::DataManager& dataManager) const
{
  double *damage, *bondDamage;
  dataManager.getData(m_damageFieldId, Field_ENUM::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);

  // Initialize damage to zero
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int nodeID = ownedIDs[iID];
    damage[nodeID] = 0.0;
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      bondDamage[bondIndex++] = 0.0;
	}
  }
}

void
PeridigmNS::CriticalStretchDamageModel::computeDamage(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  double *x, *y, *damage, *bondDamage;
  dataManager.getData(m_modelCoordinatesFieldId, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, Field_ENUM::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_damageFieldId, Field_ENUM::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);

  double trialDamage(0.0);
  int neighborhoodListIndex(0), bondIndex(0);
  int nodeId, numNeighbors, neighborID, iID, iNID;
  double nodeInitialX[3], nodeCurrentX[3], initialDistance, currentDistance, relativeExtension, totalDamage;

  // Update the bond damage
  // Break bonds if the extension is greater than the critical extension

  for(iID=0 ; iID<numOwnedPoints ; ++iID){
	nodeId = ownedIDs[iID];
	nodeInitialX[0] = x[nodeId*3];
    nodeInitialX[1] = x[nodeId*3+1];
    nodeInitialX[2] = x[nodeId*3+2];
	nodeCurrentX[0] = y[nodeId*3];
	nodeCurrentX[1] = y[nodeId*3+1];
	nodeCurrentX[2] = y[nodeId*3+2];
	numNeighbors = neighborhoodList[neighborhoodListIndex++];
	for(iNID=0 ; iNID<numNeighbors ; ++iNID){
	  neighborID = neighborhoodList[neighborhoodListIndex++];
      initialDistance = 
        distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
                 x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
      currentDistance = 
        distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
                 y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
      relativeExtension = (currentDistance - initialDistance)/initialDistance;
      trialDamage = 0.0;
      if(relativeExtension > m_criticalStretch)
        trialDamage = 1.0;
      if(trialDamage > bondDamage[bondIndex]){
        bondDamage[bondIndex] = trialDamage;
      }
      bondIndex += 1;
    }
  }

  //  Update the element damage (percent of bonds broken)

  neighborhoodListIndex = 0;
  bondIndex = 0;
  for(iID=0 ; iID<numOwnedPoints ; ++iID){
	nodeId = ownedIDs[iID];
	numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
	totalDamage = 0.0;
	for(iNID=0 ; iNID<numNeighbors ; ++iNID){
	  totalDamage += bondDamage[bondIndex++];
	}
	if(numNeighbors > 0)
	  totalDamage /= numNeighbors;
	else
	  totalDamage = 0.0;
 	damage[nodeId] = totalDamage;
  }
}
