/*! \file Peridigm_VectorPoissonMaterial.cpp */

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

#include "Peridigm_VectorPoissonMaterial.hpp"
#include "Peridigm_Field.hpp"
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

PeridigmNS::VectorPoissonMaterial::VectorPoissonMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_horizon(0.0),
    m_coefficient(0.0),
    m_volumeFieldId(-1),
    m_modelCoordinatesFieldId(-1),
    m_coordinatesFieldId(-1),
    m_forceDensityFieldId(-1)
{
  m_horizon = params.get<double>("Horizon");
  m_coefficient = params.get<double>("Coefficient");

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR,      PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId             = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR,      PeridigmField::TWO_STEP, "Force_Density");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
}

PeridigmNS::VectorPoissonMaterial::~VectorPoissonMaterial()
{
}

void
PeridigmNS::VectorPoissonMaterial::initialize(const double dt,
					      const int numOwnedPoints,
					      const int* ownedIDs,
					      const int* neighborhoodList,
					      PeridigmNS::DataManager& dataManager)
{
}

void
PeridigmNS::VectorPoissonMaterial::computeForce(const double dt,
						const int numOwnedPoints,
						const int* ownedIDs,
						const int* neighborhoodList,
						PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *force;

  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);

  int neighborhoodListIndex(0);
  int numNeighbors, neighborID, iID, iNID;
  double nodeInitialPosition[3], nodeDisplacement[3], initialDistance;

  for(iID=0 ; iID<numOwnedPoints ; ++iID){

    nodeInitialPosition[0] = x[iID*3];
    nodeInitialPosition[1] = x[iID*3+1];
    nodeInitialPosition[2] = x[iID*3+2];
    nodeDisplacement[0] = y[iID*3]   - nodeInitialPosition[0];
    nodeDisplacement[1] = y[iID*3+1] - nodeInitialPosition[1];
    nodeDisplacement[2] = y[iID*3+2] - nodeInitialPosition[2];

    numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(iNID=0 ; iNID<numNeighbors ; ++iNID){
      neighborID = neighborhoodList[neighborhoodListIndex++];
      initialDistance = distance(nodeInitialPosition[0], nodeInitialPosition[1], nodeInitialPosition[2],
				 x[neighborID*3], x[neighborID*3+1], x[neighborID*3+2]);

      // TODO:  Write me.


    }
  }

}
