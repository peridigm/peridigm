/*! \file Peridigm_DamageModel.cpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#include "Peridigm_CriticalStretchDamageModel.hpp"
#include <Teuchos_TestForException.hpp>
#include "PdMaterialUtilities.h"

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
PeridigmNS::CriticalStretchDamageModel::initialize(const Epetra_Vector& u,
                                                   const Epetra_Vector& v,
                                                   const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   double* bondState,
                                                   PeridigmNS::DataManager& dataManager,
                                                   Epetra_Vector& force) const
{
}

void
PeridigmNS::CriticalStretchDamageModel::computeDamage(const Epetra_Vector& u,
                                                      const Epetra_Vector& v,
                                                      const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      double* bondState,
                                                      PeridigmNS::DataManager& dataManager,
                                                      Epetra_Vector& force) const
{
  int vectorLength = dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->MyLength();
  double *x, *y;
  dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);

  // Update the bondState
  // Break bonds if the extension is greater than the critical extension
  double trialDamage = 0.0;
  int neighborhoodListIndex = 0;
  int bondStateIndex = 0;
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
	  if(trialDamage > bondState[bondStateIndex]){
		bondState[bondStateIndex] = trialDamage;
      }
	  bondStateIndex += 1;
	}
  }
}
