/* \file Peridigm_SpeciesConcentrationMaterial.cpp */

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

#include "Peridigm_SpeciesConcentrationMaterial.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_Constants.hpp"

PeridigmNS::SpeciesConcentrationMaterial::SpeciesConcentrationMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_horizon(0.0),
    m_coefficient(0.0),
    m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1),
    m_concentrationFieldId(-1),
    m_fluxDivergenceFieldId(-1)
{
  m_horizon = params.get<double>("Horizon");
  m_coefficient = params.get<double>("Coefficient");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_temperatureFieldId             = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_concentrationFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Concentration");
  m_fluxDivergenceFieldId          = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR,      PeridigmField::TWO_STEP, "Flux_Divergence");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_temperatureFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_concentrationFieldId);
  m_fieldIds.push_back(m_fluxDivergenceFieldId);
}

PeridigmNS::SpeciesConcentrationMaterial::~SpeciesConcentrationMaterial()
{
}

void
PeridigmNS::SpeciesConcentrationMaterial::initialize(const double dt,
                                          const int numOwnedPoints,
                                          const int* ownedIDs,
                                          const int* neighborhoodList,
                                          PeridigmNS::DataManager& dataManager)
{}

void
PeridigmNS::SpeciesConcentrationMaterial::computeFluxDivergence(const double dt,
                                                     const int numOwnedPoints,
                                                     const int* ownedIDs,
                                                     const int* neighborhoodList,
                                                     PeridigmNS::DataManager& dataManager) const
{
  // Zero out the flux divergence
  dataManager.getData(m_fluxDivergenceFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *volume, *temperature, *modelCoord, *concentration, *fluxDivergence;

  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(m_temperatureFieldId, PeridigmField::STEP_NP1)->ExtractView(&temperature);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&modelCoord);
  dataManager.getData(m_concentrationFieldId, PeridigmField::STEP_NP1)->ExtractView(&concentration);
  dataManager.getData(m_fluxDivergenceFieldId, PeridigmField::STEP_NP1)->ExtractView(&fluxDivergence);

  const double pi = value_of_pi();

  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    double nodeConcentration = concentration[iID];
    double nodeModelCoordX = modelCoord[iID*3];
    double nodeModelCoordY = modelCoord[iID*3+1];
    double nodeModelCoordZ = modelCoord[iID*3+2];
    double nodeTemperature = temperature[iID];
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      double neighborVolume = volume[neighborID];
      double initialDistance = distance(nodeModelCoordX, nodeModelCoordY, nodeModelCoordZ, modelCoord[neighborID*3], modelCoord[neighborID*3+1], modelCoord[neighborID*3+2]);
      double concentrationDifference = concentration[neighborID] - nodeConcentration;
      double mu = 1.0; // PLACEHOLDER
      double coefficient = m_coefficient * (1.0 + 0.0001*nodeTemperature); // PLACEHOLDER
      double contribution_to_flux_divergence = mu * coefficient * (concentrationDifference / initialDistance) * neighborVolume;
      TEUCHOS_TEST_FOR_TERMINATION(!std::isfinite(contribution_to_flux_divergence), "**** NaN detected in SpeciesConcentrationMaterial::computeFluxDivergence().\n");
      fluxDivergence[iID] += contribution_to_flux_divergence;
    }
  }
}
