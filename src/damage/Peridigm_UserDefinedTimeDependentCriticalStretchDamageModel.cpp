/*! \file Peridigm_UserDefinedTimeDependentCriticalStretchDamageModel.cpp */

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

#include "Peridigm_UserDefinedTimeDependentCriticalStretchDamageModel.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel::UserDefinedTimeDependentCriticalStretchDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params), m_applyThermalStrains(false), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_damageFieldId(-1), m_bondDamageFieldId(-1), m_deltaTemperatureFieldId(-1)
{

  functiondmg = params.get<string>("Time Dependent Critical Stretch");
  
  // set up RTCompiler
  rtcFunction = Teuchos::rcp<PG_RuntimeCompiler::Function>(new PG_RuntimeCompiler::Function(2, "rtcUserDefinedTimeDependentShortRangeForceContactModel"));
  rtcFunction->addVar("double", "t");
  rtcFunction->addVar("double", "value");
    

  if(params.isParameter("Thermal Expansion Coefficient")){
    m_alpha = params.get<double>("Thermal Expansion Coefficient");
    m_applyThermalStrains = true;
  }

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  //m_stepFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Step");
  m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");
  if(m_applyThermalStrains)
    m_deltaTemperatureFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
    
  //PeridigmNS::PeridigmField::Step step = PeridigmField::STEP_NONE;  
  //step = PeridigmField::STEP_NP1;

  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  //m_fieldIds.push_back(m_stepFieldId);
  
  if(m_applyThermalStrains)
    m_fieldIds.push_back(m_deltaTemperatureFieldId);
}

PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel::~UserDefinedTimeDependentCriticalStretchDamageModel()
{
}

void PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel::evaluateParserDmg(double & currentValue, double & previousValue, const double & timeCurrent, const double & timePrevious){
  
  string rtcFunctionString = functiondmg;
  if(rtcFunctionString.find("value") == string::npos)
    rtcFunctionString = "value = " + rtcFunctionString;
  bool success = rtcFunction->addBody(rtcFunctionString);
  if(!success){
    string msg = "\n**** Error:  rtcFunction->addBody(functiondmg) returned nonzero error code in UserDefinedTimeDependentCriticalStretchDamageModel::evaluateParserDmg().\n";
    msg += "**** " + rtcFunction->getErrors() + "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
  }

  // set the return value to 0.0
  if(success)
    success = rtcFunction->varValueFill(1, 0.0);
  // evaluate at previous time
  if(success)
    success = rtcFunction->varValueFill(0, timePrevious);
  if(success)
    success = rtcFunction->execute();
  if(success)
    previousValue = rtcFunction->getValueOfVar("value");
  // evaluate at current time
  if(success)
    success = rtcFunction->varValueFill(0, timeCurrent);
  if(success)
    success = rtcFunction->execute();
  if(success)
    currentValue = rtcFunction->getValueOfVar("value");
  if(!success){
    string msg = "\n**** Error in UserDefinedTimeDependentCriticalStretchDamageModel::evaluateParserDmg().\n";
    msg += "**** " + rtcFunction->getErrors() + "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!success, msg);
  }

  m_criticalStretch = currentValue;
  
}

void
PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel::initialize(const double dt,
                                                   const int numOwnedPoints,
                                                   const int* ownedIDs,
                                                   const int* neighborhoodList,
                                                   PeridigmNS::DataManager& dataManager) const
{

  double *damage, *bondDamage, *step;
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  //dataManager.getData(m_stepFieldId, PeridigmField::STEP_NP1)->ExtractView(&step);
  
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
PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel::computeDamage(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* neighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const 
{
  double *x, *y, *damage, *bondDamageN, *bondDamageNP1, *deltaTemperature, *step;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamageN);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);
  //dataManager.getData(m_stepFieldId, PeridigmField::STEP_NP1)->ExtractView(&step);
  deltaTemperature = NULL;
  if(m_applyThermalStrains)
    dataManager.getData(m_deltaTemperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&deltaTemperature);
    

  double trialDamage(0.0);
  int neighborhoodListIndex(0), bondIndex(0);
  int nodeId, numNeighbors, neighborID, iID, iNID;
  double nodeInitialX[3], nodeCurrentX[3], initialDistance, currentDistance, relativeExtension, totalDamage;

  // Set the bond damage to the previous value
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1) =
    dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N);

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
      if(m_applyThermalStrains)
        currentDistance -= m_alpha*deltaTemperature[nodeId]*initialDistance;
      relativeExtension = (currentDistance - initialDistance)/initialDistance;
      trialDamage = 0.0;
      if(relativeExtension > m_criticalStretch)
        trialDamage = 1.0;
      if(trialDamage > bondDamageNP1[bondIndex]){
        bondDamageNP1[bondIndex] = trialDamage;
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
      totalDamage += bondDamageNP1[bondIndex++];
    }
    if(numNeighbors > 0)
      totalDamage /= numNeighbors;
    else
      totalDamage = 0.0;
    damage[nodeId] = totalDamage;
  }
}
