/*! \file Peridigm_DiffusionMaterial.cpp */

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

#include "Peridigm_DiffusionMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_Constants.hpp"
#ifdef PERIDIGM_IMPROVED_QUADRATURE
  #include <gsl/gsl_linalg.h>
  #include <gsl/gsl_cblas.h>
  #include "util/nathelpers.h"
  #include "scalarDiffusionQuad/nonlocQuad.h"
#endif

PeridigmNS::DiffusionMaterial::DiffusionMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_horizon(0.0),
    m_useImprovedQuadrature(false),
    m_coefficient(0.0),
    m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1),
    m_temperatureFieldId(-1),
    m_fluxDivergenceFieldId(-1),
    m_quadratureWeightsFieldId(-1)
{
  m_horizon = params.get<double>("Horizon");
  m_coefficient = params.get<double>("Coefficient");
  if (params.isParameter("Use Improved Quadrature")) {
    m_useImprovedQuadrature = params.get<bool>("Use Improved Quadrature");
  }

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_temperatureFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Temperature");
  m_fluxDivergenceFieldId          = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Flux_Divergence");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_temperatureFieldId);
  m_fieldIds.push_back(m_fluxDivergenceFieldId);

  if (m_useImprovedQuadrature) {
    m_quadratureWeightsFieldId     = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Quadrature_Weights");
    m_fieldIds.push_back(m_quadratureWeightsFieldId);
  }

#ifndef PERIDIGM_IMPROVED_QUADRATURE
  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_useImprovedQuadrature, "**** Error:  Improved quadrature not available.  Recompile Peridigm with USE_IMPROVED_QUADRATURE:BOOL=ON\n");
#endif
}

PeridigmNS::DiffusionMaterial::~DiffusionMaterial()
{
}

void
PeridigmNS::DiffusionMaterial::initialize(const double dt,
                                          const int numOwnedPoints,
                                          const int* ownedIDs,
                                          const int* neighborhoodList,
                                          PeridigmNS::DataManager& dataManager)
{
#ifdef PERIDIGM_IMPROVED_QUADRATURE
  if (m_useImprovedQuadrature) {
    double *x, *volume, *quadratureWeights;
    dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    dataManager.getData(m_quadratureWeightsFieldId, PeridigmField::STEP_NONE)->ExtractView(&quadratureWeights);

    int nodeId, numNeighbors, neighborId, neighborhoodListIndex, bondListIndex;
    double neighborVolume;

    triple<double> X;
    double delta = m_horizon;
    int Porder = 2;

    neighborhoodListIndex = 0;
    bondListIndex = 0;
    for(int iID=0 ; iID<numOwnedPoints ; ++iID){
      nodeId = ownedIDs[iID];
      numNeighbors = neighborhoodList[neighborhoodListIndex++];
      X = vec3(x[nodeId*3], x[nodeId*3+1], x[nodeId*3+2]);
      std::vector< triple<double> > bondList;
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
        neighborId = neighborhoodList[neighborhoodListIndex++];
        vec3 neighbor(x[neighborId*3], x[neighborId*3+1], x[neighborId*3+2]);
        bondList.push_back(neighbor);
      }
      bondList.push_back(X);
      // Instantiate class for generating quadrature weights
      nonlocalDiffusionQuad quad(bondList, X, delta, Porder);
      std::vector<double> weights = quad.getWeights();
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
        quadratureWeights[bondListIndex++] = weights.at(iNID);
      }
    }
  }
#endif
}

void
PeridigmNS::DiffusionMaterial::computeFluxDivergence(const double dt,
                                                     const int numOwnedPoints,
                                                     const int* ownedIDs,
                                                     const int* neighborhoodList,
                                                     PeridigmNS::DataManager& dataManager) const
{
  Teuchos::RCP< std::map< std::string, std::vector<int> > > node_sets = bc_manager_->getNodeSets();
  //std::cout << "DEBUGGING number of node sets " << node_sets->size() << std::endl;

  // Zero out the flux divergence
  dataManager.getData(m_fluxDivergenceFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *volume, *modelCoord, *temperature, *fluxDivergence, *quadratureWeights;

  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
  dataManager.getData(m_temperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
  dataManager.getData(m_fluxDivergenceFieldId, PeridigmField::STEP_NP1)->ExtractView(&fluxDivergence);
  if (m_useImprovedQuadrature) {
    dataManager.getData(m_quadratureWeightsFieldId, PeridigmField::STEP_NONE)->ExtractView(&quadratureWeights);
  }

  int neighborhoodListIndex(0);
  int numNeighbors, neighborID, iID, iNID, bondListIndex;
  double nodeInitialPosition[3], initialDistance, quadWeight;
  double kernel, nodeTemperature, temperatureDifference, nodeFluxDivergence, neighborFluxDivergence;

  const double pi = value_of_pi();

  bondListIndex = 0;
  for(iID=0 ; iID<numOwnedPoints ; ++iID){
    nodeTemperature = temperature[iID];
    nodeInitialPosition[0] = modelCoord[iID*3];
    nodeInitialPosition[1] = modelCoord[iID*3+1];
    nodeInitialPosition[2] = modelCoord[iID*3+2];
    numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(iNID=0 ; iNID<numNeighbors ; ++iNID){
      neighborID = neighborhoodList[neighborhoodListIndex++];
      quadWeight = volume[neighborID];
      if (m_useImprovedQuadrature) {
        quadWeight = quadratureWeights[bondListIndex++];
      }
      initialDistance = distance(nodeInitialPosition[0], nodeInitialPosition[1], nodeInitialPosition[2],
                                 modelCoord[neighborID*3], modelCoord[neighborID*3+1], modelCoord[neighborID*3+2]);
      kernel = 6.0/(pi*m_horizon*m_horizon*m_horizon*m_horizon*initialDistance);
      temperatureDifference = temperature[neighborID] - nodeTemperature;
      nodeFluxDivergence = m_coefficient*kernel*temperatureDifference*quadWeight;
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(nodeFluxDivergence), "**** NaN detected in DiffusionMaterial::computeFluxDivergence().\n");
      fluxDivergence[iID] += nodeFluxDivergence;
    }
  }
}
