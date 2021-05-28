/*! \file Peridigm_Compute_JIntegral.cpp */

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
// Michael L. Parks      parks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include <vector>

#include "Peridigm_Compute_JIntegral.hpp"
#include "Peridigm_Field.hpp"
#include <iostream>

//! Standard constructor.
PeridigmNS::Compute_JIntegral::Compute_JIntegral(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                 Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                 Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_), m_blockId(-1), m_volumeFieldId(-1), m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1),
    m_velocityFieldId(-1), m_velocityGradientFieldId(-1), m_storedElasticEnergyDensityFieldId(-1), m_storedElasticEnergyFieldId(-1), m_bondDamageFieldId(-1)
{
  std::cout << "WARNING:  The J-Integral compute class is currently not implemented (work in progress).\n" << std::endl;

  m_blockId = params->get<int>("Block ID");
  m_dt = params->get<double>("dt");

  FieldManager& fieldManager = FieldManager::self();
  m_volumeFieldId = fieldManager.getFieldId("Volume");
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  m_velocityFieldId = fieldManager.getFieldId("Velocity");
  m_velocityGradientFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::TWO_STEP, "Velocity_Gradient");
  m_storedElasticEnergyDensityFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Stored_Elastic_Energy_Density");
  m_storedElasticEnergyFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Stored_Elastic_Energy");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_velocityFieldId);
  m_fieldIds.push_back(m_velocityGradientFieldId);
  m_fieldIds.push_back(m_storedElasticEnergyDensityFieldId);
  m_fieldIds.push_back(m_storedElasticEnergyFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
}

//! Destructor.
PeridigmNS::Compute_JIntegral::~Compute_JIntegral(){}

//! Compute the J integral
int PeridigmNS::Compute_JIntegral::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const
{
  // Isolate a single block over which the J integral will be computed
  std::vector<PeridigmNS::Block>::iterator block = blocks->end();
  for(std::vector<PeridigmNS::Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    if (blockIt->getID() == m_blockId) {
      block = blockIt;
    }
  }
  TEUCHOS_TEST_FOR_TERMINATION(block == blocks->end(), "**** Error:  Compute_JIntegral, specified block_id not found.\n");

  // Extract data for the block
  Teuchos::RCP<const PeridigmNS::Material> materialModel = block->getMaterialModel();
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = block->getNeighborhoodData();
  Teuchos::RCP<PeridigmNS::DataManager> dataManager = block->getDataManager();
  int numOwnedPoints = neighborhoodData->NumOwnedPoints();
  int* const ownedIDs = neighborhoodData->OwnedIDs();
  int* const neighborhoodList = neighborhoodData->NeighborhoodList();

  // Extract pointers to the underlying data
  double *volume, *modelCoord, *coord, *velocityStepN, *velocityStepNP1, *velocityGradient, *storedElasticEnergyDensity, *storedElasticEnergy, *bondDamage;
  dataManager->getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
  dataManager->getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&coord);
  dataManager->getData(m_velocityFieldId, PeridigmField::STEP_N)->ExtractView(&velocityStepN);
  dataManager->getData(m_velocityFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocityStepNP1);
  dataManager->getData(m_velocityGradientFieldId, PeridigmField::STEP_NP1)->ExtractView(&velocityGradient);
  dataManager->getData(m_storedElasticEnergyDensityFieldId, PeridigmField::STEP_NONE)->ExtractView(&storedElasticEnergyDensity);
  dataManager->getData(m_storedElasticEnergyFieldId, PeridigmField::STEP_NONE)->ExtractView(&storedElasticEnergy);
  dataManager->getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);

  // Compute the stored elastic energy (strain energy) in the block
  // This calls through to the material model's computeStoredElasticEnergyDensity() function, which must be implemented
  // for the specific material model that has been assigned to the block
  block->getMaterialModel()->computeStoredElasticEnergyDensity(m_dt, numOwnedPoints, ownedIDs, neighborhoodList, *dataManager);

  // Convert stored elastic energy density to stored elastic energy
  for(int i=0 ; i<numOwnedPoints ; ++i) {
    storedElasticEnergy[i] = storedElasticEnergyDensity[i] * volume[i];
  }

  // Compute the velocity gradient
  int neighborhoodListIndex = 0;
  int bondListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    double nodeModelCoord[3];
    nodeModelCoord[0] = modelCoord[iID*3];
    nodeModelCoord[1] = modelCoord[iID*3+1];
    nodeModelCoord[2] = modelCoord[iID*3+2];
    double nodeCoord[3];
    nodeCoord[0] = coord[iID*3];
    nodeCoord[1] = coord[iID*3+1];
    nodeCoord[2] = coord[iID*3+2];
    double nodeVelocityStepN[3];
    nodeVelocityStepN[0] = velocityStepN[iID*3];
    nodeVelocityStepN[1] = velocityStepN[iID*3+1];
    nodeVelocityStepN[2] = velocityStepN[iID*3+2];
    double nodeVelocityStepNP1[3];
    nodeVelocityStepNP1[0] = velocityStepNP1[iID*3];
    nodeVelocityStepNP1[1] = velocityStepNP1[iID*3+1];
    nodeVelocityStepNP1[2] = velocityStepNP1[iID*3+2];
    double* nodeVelocityGradient = &velocityGradient[iID*9];
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      double neighborModelCoord[3];
      neighborModelCoord[0] = modelCoord[neighborID*3];
      neighborModelCoord[1] = modelCoord[neighborID*3+1];
      neighborModelCoord[2] = modelCoord[neighborID*3+2];
      double neighborCoord[3];
      neighborCoord[0] = coord[neighborID*3];
      neighborCoord[1] = coord[neighborID*3+1];
      neighborCoord[2] = coord[neighborID*3+2];
      double neighborVelocityStepN[3];
      neighborVelocityStepN[0] = velocityStepN[neighborID*3];
      neighborVelocityStepN[1] = velocityStepN[neighborID*3+1];
      neighborVelocityStepN[2] = velocityStepN[neighborID*3+2];
      double neighborVelocityStepNP1[3];
      neighborVelocityStepNP1[0] = velocityStepNP1[neighborID*3];
      neighborVelocityStepNP1[1] = velocityStepNP1[neighborID*3+1];
      neighborVelocityStepNP1[2] = velocityStepNP1[neighborID*3+2];

      // COMPUTE AND STORE THE VELOCITY GRADIENT
      // delta t is stored in m_dt, which unfortunately must be passed in through the input deck for now

      nodeVelocityGradient[0] = 0.0; // xx component
      nodeVelocityGradient[1] = 0.0; // xy component
      nodeVelocityGradient[2] = 0.0; // xz component
      nodeVelocityGradient[3] = 0.0; // yx component
      nodeVelocityGradient[4] = 0.0; // yy component
      nodeVelocityGradient[5] = 0.0; // yx component
      nodeVelocityGradient[6] = 0.0; // zx component
      nodeVelocityGradient[7] = 0.0; // zy component
      nodeVelocityGradient[8] = 0.0; // zz component
    }
  }

  return(0);
}
