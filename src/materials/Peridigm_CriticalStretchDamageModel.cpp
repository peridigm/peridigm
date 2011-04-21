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
#include <Teuchos_TestForException.hpp>

using namespace std;

PeridigmNS::CriticalStretchDamageModel::CriticalStretchDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params)
{
  m_criticalStretch = params.get<double>("Critical Stretch");

  // set up vector of variable specs
  m_variableSpecs = Teuchos::rcp(new vector<Field_NS::FieldSpec>);
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
}

void
PeridigmNS::CriticalStretchDamageModel::computeDamage(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  int vectorLength = dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->MyLength();
  double *x, *y, *bondDamage;
  dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&bondDamage);

  // Update the bond damage
  // Break bonds if the extension is greater than the critical extension
  double trialDamage = 0.0;
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int nodeID = ownedIDs[iID];
	TEST_FOR_EXCEPT_MSG(nodeID*3+2 >= vectorLength, "Invalid neighbor list / x vector\n");
	double nodeInitialX[3] = { x[nodeID*3],
							   x[nodeID*3+1],
							   x[nodeID*3+2] };
	double nodeCurrentX[3] = { y[nodeID*3],
							   y[nodeID*3+1],
							   y[nodeID*3+2] };
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
	  int neighborID = neighborhoodList[neighborhoodListIndex++];
	  TEST_FOR_EXCEPT_MSG(neighborID < 0, "Invalid neighbor list\n");
	  TEST_FOR_EXCEPT_MSG(neighborID*3+2 >= vectorLength, "Invalid neighbor list / initial x vector\n");
	  double initialDistance = 
		distance(nodeInitialX[0], nodeInitialX[1], nodeInitialX[2],
				 x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);
	  double currentDistance = 
		distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
				 y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
	  double relativeExtension = (currentDistance - initialDistance)/initialDistance;
	  trialDamage = 0.0;
	  if(relativeExtension > m_criticalStretch)
		trialDamage = 1.0;
	  if(trialDamage > bondDamage[bondIndex]){
		bondDamage[bondIndex] = trialDamage;
      }
	  bondIndex += 1;
	}
  }
}
