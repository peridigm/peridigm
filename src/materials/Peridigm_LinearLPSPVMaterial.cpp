/*! \file Peridigm_LinearLPSPVMaterial.cpp */

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

#include "Peridigm_LinearLPSPVMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "elastic_pv.h"     // for weighted volume
#include "linear_lps_pv.h"  // for internal force
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_Comm.h>
#ifdef PERIDIGM_IMPROVED_QUADRATURE
  #include <gsl/gsl_linalg.h>
  #include <gsl/gsl_cblas.h>
  #include "util/nathelpers.h"
  #include "stateBasedQuad/nonlocQuadStateBasedFast.h"
#endif

PeridigmNS::LinearLPSPVMaterial::LinearLPSPVMaterial(const Teuchos::ParameterList& params)
  : Material(params), m_pid(-1), m_verbose(false),
    m_bulkModulus(0.0), m_shearModulus(0.0), m_density(0.0), m_horizon(0.0), m_useImprovedQuadrature(false),
    m_omega(PeridigmNS::InfluenceFunction::self().getInfluenceFunction()),
    m_useAnalyticWeightedVolume(false), m_analyticWeightedVolume(0.0), m_usePartialVolume(false),
    m_volumeFieldId(-1), m_damageFieldId(-1), m_weightedVolumeFieldId(-1), m_dilatationFieldId(-1), m_modelCoordinatesFieldId(-1),
    m_coordinatesFieldId(-1), m_forceDensityFieldId(-1), m_bondDamageFieldId(-1),
    m_selfVolumeFieldId(-1), m_neighborVolumeFieldId(-1), m_quadratureWeightsFieldId(-1)
{
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  if(params.isParameter("Verbose")) {
    m_verbose = params.get<bool>("Verbose");
  }
  if (params.isParameter("Use Improved Quadrature")) {
    m_useImprovedQuadrature = params.get<bool>("Use Improved Quadrature");
  }
  if(params.isParameter("Use Partial Volume")) {
    m_usePartialVolume = params.get<bool>("Use Partial Volume");
  }
  if(params.isParameter("Use Analytic Weighted Volume")) {
    m_useAnalyticWeightedVolume = params.get<bool>("Use Analytic Weighted Volume");
    if(m_useAnalyticWeightedVolume) {
      m_analyticWeightedVolume = params.get<double>("Analytic Weighted Volume");
    }
  }

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_damageFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_weightedVolumeFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
  m_dilatationFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_influenceFunctionFieldId       = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::CONSTANT, "Influence_Function");
  m_bondDamageFieldId              = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);
  m_fieldIds.push_back(m_dilatationFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_influenceFunctionFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);

  if (m_useImprovedQuadrature) {
    m_quadratureWeightsFieldId     = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Quadrature_Weights");
    m_fieldIds.push_back(m_quadratureWeightsFieldId);
  }

#ifndef PERIDIGM_PV
  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_usePartialVolume, "**** Error:  Partial volumes not available.  Recompile Peridigm with USE_PV:BOOL=ON\n");
#endif

#ifndef PERIDIGM_IMPROVED_QUADRATURE
  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_useImprovedQuadrature, "**** Error:  Improved quadrature not available.  Recompile Peridigm with USE_IMPROVED_QUADRATURE:BOOL=ON\n");
#endif
}

PeridigmNS::LinearLPSPVMaterial::~LinearLPSPVMaterial()
{
}

void
PeridigmNS::LinearLPSPVMaterial::initialize(const double dt,
					    const int numOwnedPoints,
					    const int* ownedIDs,
					    const int* neighborhoodList,
					    PeridigmNS::DataManager& dataManager)
{
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_pid = dataManager.getOwnedScalarPointMap()->Comm().MyPID();

  // Extract pointers to the underlying data
  double *xOverlap, *influenceFunctionValues;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(m_influenceFunctionFieldId, PeridigmField::STEP_NONE)->ExtractView(&influenceFunctionValues);
  MATERIAL_EVALUATION::computeAndStoreInfluenceFunctionValues(xOverlap,
                                                              influenceFunctionValues,
                                                              numOwnedPoints,
                                                              neighborhoodList,
                                                              m_horizon);

  if(m_usePartialVolume){
    m_selfVolumeFieldId = fieldManager.getFieldId("Self_Volume");
    m_neighborVolumeFieldId = fieldManager.getFieldId("Neighbor_Volume");
  }

  // Extract pointers to the underlying data
  double *cellVolumeOverlap, *weightedVolume;
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolumeOverlap);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);

  // Extract pointers to partial volume data
  double *selfVolume(0), *neighborVolume(0);
  if(m_usePartialVolume){
    if(m_verbose && m_pid == 0)
      std::cout << "\nLinearLPS Material Model applying partial volume.\n" << std::endl;
    dataManager.getData(m_selfVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&selfVolume);
    dataManager.getData(m_neighborVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&neighborVolume);
  }
  else{
    if(m_verbose && m_pid == 0)
      std::cout << "\nLinearLPS Material Model not applying partial volume.\n" << std::endl;
  }

  if(m_useAnalyticWeightedVolume){
    if(m_verbose && m_pid == 0)
      std::cout << "\nLinearLPS Material Model applying analytic weighted volume.\n" << std::endl;
    dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->PutScalar(m_analyticWeightedVolume);
  }
  else{
    if(m_verbose && m_pid == 0)
      std::cout << "\nLinearLPS Material Model applying numerical weighted volume.\n" << std::endl;
    MATERIAL_EVALUATION::computeWeightedVolumePV(xOverlap,
						 cellVolumeOverlap,
						 selfVolume,
						 neighborVolume,
						 influenceFunctionValues,
						 weightedVolume,
						 numOwnedPoints,
						 neighborhoodList,
						 m_horizon);
  }

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
      nonlocalQuadStateBased quad(bondList, X, delta, Porder);
      std::vector<double> weights = quad.getWeightsState();
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
        quadratureWeights[bondListIndex++] = weights.at(iNID);
      }
    }
  }
#endif
}

void
PeridigmNS::LinearLPSPVMaterial::precompute(const double dt,
                                            const int numOwnedPoints,
                                            const int* ownedIDs,
                                            const int* neighborhoodList,
                                            PeridigmNS::DataManager& dataManager) const
{
  // Zero out the dilatation
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  double *x, *y, *volume, *weightedVolume, *dilatation, *influenceFunctionValues, *damage, *neighborVolume(0), *quadratureWeights(0);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(m_influenceFunctionFieldId, PeridigmField::STEP_NONE)->ExtractView(&influenceFunctionValues);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  if(m_usePartialVolume){
    dataManager.getData(m_neighborVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&neighborVolume);
  }
  if (m_useImprovedQuadrature) {
    dataManager.getData(m_quadratureWeightsFieldId, PeridigmField::STEP_NONE)->ExtractView(&quadratureWeights);
  }

  const double *xSelf, *ySelf, *xNeighbor, *yNeighbor;
  double u[3], uNeighbor[3], dotProduct, zeta[3], volNeighbor, omega;
  int numNeighbors, neighborId, neighborhoodListIndex(0), bondIndex(0);

  for(int i_pt=0; i_pt<numOwnedPoints; i_pt++){
    dilatation[i_pt] = 0.0;
    numNeighbors = neighborhoodList[neighborhoodListIndex++];
    xSelf = &x[3*i_pt];
    ySelf = &y[3*i_pt];
    for(int i_neighbor=0; i_neighbor<numNeighbors; i_neighbor++){
      neighborId = neighborhoodList[neighborhoodListIndex++];
      xNeighbor = &x[3*neighborId];
      yNeighbor = &y[3*neighborId];
      if(m_usePartialVolume)
        volNeighbor = neighborVolume[bondIndex];
      else
        volNeighbor = volume[neighborId];
      for(int i=0 ; i<3 ; ++i){
        zeta[i] = xNeighbor[i] - xSelf[i];
        u[i] = ySelf[i] - xSelf[i];
        uNeighbor[i] = yNeighbor[i] - xNeighbor[i];
      }
      omega = influenceFunctionValues[bondIndex];
      dotProduct = zeta[0]*(uNeighbor[0]-u[0]) + zeta[1]*(uNeighbor[1]-u[1]) + zeta[2]*(uNeighbor[2]-u[2]);
      dilatation[i_pt] += omega*(1.0 - damage[bondIndex])*dotProduct*volNeighbor;
      bondIndex += 1;
    }
    if(numNeighbors > 0){
      dilatation[i_pt] *= 3.0/weightedVolume[i_pt];
    }
  }
}

void
PeridigmNS::LinearLPSPVMaterial::computeForce(const double dt,
                                              const int numOwnedPoints,
                                              const int* ownedIDs,
                                              const int* neighborhoodList,
                                              PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *volume, *weightedVolume, *dilatation, *influenceFunctionValues, *damage, *force, *selfVolume(0), *neighborVolume(0), *quadratureWeights(0);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(m_influenceFunctionFieldId, PeridigmField::STEP_NONE)->ExtractView(&influenceFunctionValues);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
  if(m_usePartialVolume){
    dataManager.getData(m_selfVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&selfVolume);
    dataManager.getData(m_neighborVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&neighborVolume);
  }
  if (m_useImprovedQuadrature) {
    dataManager.getData(m_quadratureWeightsFieldId, PeridigmField::STEP_NONE)->ExtractView(&quadratureWeights);
  }

  const double *xSelf, *ySelf, *xNeighbor, *yNeighbor;
  double u[3], uNeighbor[3], dotProduct, zeta[3], volSelf, volNeighbor, matVec[3], dyadicProduct[3][3];
  int numNeighbors, neighborId, neighborhoodListIndex(0), bondIndex(0);

  for(int i_pt=0; i_pt<numOwnedPoints; i_pt++){
    numNeighbors = neighborhoodList[neighborhoodListIndex++];
    xSelf = &x[3*i_pt];
    ySelf = &y[3*i_pt];
    for(int i_neighbor=0; i_neighbor<numNeighbors; i_neighbor++){
      neighborId = neighborhoodList[neighborhoodListIndex++];
      xNeighbor = &x[3*neighborId];
      yNeighbor = &y[3*neighborId];
      if(m_usePartialVolume){
        volSelf = selfVolume[bondIndex];
        volNeighbor = neighborVolume[bondIndex];
      }
      else{
        volSelf = volume[i_pt];
        volNeighbor = volume[neighborId];
      }
      for(int i=0 ; i<3 ; ++i){
        zeta[i] = xNeighbor[i] - xSelf[i];
        u[i] = ySelf[i] - xSelf[i];
        uNeighbor[i] = yNeighbor[i] - xNeighbor[i];
      }
      double normZetaSquared = zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2];
      double omega = influenceFunctionValues[bondIndex];
      double temp1 = (9.0*m_bulkModulus - 15.0*m_shearModulus)*omega*(dilatation[i_pt])/(3.0*(weightedVolume[i_pt]));
      double temp2 = 15.0*m_shearModulus*omega/((weightedVolume[i_pt])*normZetaSquared);
      for(int i=0 ; i<3 ; i++){
        for(int j=0 ; j<3 ; j++){
          dyadicProduct[i][j] = zeta[i]*zeta[j];
        }
      }
      for(int i=0 ; i<3 ; ++i){
        matVec[i] = dyadicProduct[i][0]*(uNeighbor[0]-u[0]) + dyadicProduct[i][1]*(uNeighbor[1]-u[1]) + dyadicProduct[i][2]*(uNeighbor[2]-u[2]);
      }
      double fx = (1.0 - damage[bondIndex])*(temp1*zeta[0] + temp2*matVec[0]);
      double fy = (1.0 - damage[bondIndex])*(temp1*zeta[1] + temp2*matVec[1]);
      double fz = (1.0 - damage[bondIndex])*(temp1*zeta[2] + temp2*matVec[2]);
      force[3*i_pt]   += fx*volNeighbor;
      force[3*i_pt+1] += fy*volNeighbor;
      force[3*i_pt+2] += fz*volNeighbor;

#if 0
      force[3*neighborId]   -= fx*volSelf;
      force[3*neighborId+1] -= fy*volSelf;
      force[3*neighborId+2] -= fz*volSelf;
#else
      temp1 = (9.0*m_bulkModulus - 15.0*m_shearModulus)*omega*(dilatation[neighborId])/(3.0*(weightedVolume[neighborId]));
      temp2 = 15.0*m_shearModulus*omega/((weightedVolume[neighborId])*normZetaSquared);
      fx = -temp1*zeta[0] - temp2*matVec[0];
      fy = -temp1*zeta[1] - temp2*matVec[1];
      fz = -temp1*zeta[2] - temp2*matVec[2];
      force[3*i_pt]   -= fx*volNeighbor;
      force[3*i_pt+1] -= fy*volNeighbor;
      force[3*i_pt+2] -= fz*volNeighbor;
#endif

      bondIndex += 1;
    }
  }
}
