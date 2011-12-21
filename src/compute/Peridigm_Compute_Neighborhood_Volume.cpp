/*! \file Peridigm_Compute_Neighborhood_Volume.cpp */

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

#include <vector>

#include "Peridigm_Compute_Neighborhood_Volume.hpp"
#include "../core/Peridigm.hpp"

PeridigmNS::Compute_Neighborhood_Volume::Compute_Neighborhood_Volume(PeridigmNS::Peridigm *peridigm_ ){peridigm = peridigm_;}

PeridigmNS::Compute_Neighborhood_Volume::~Compute_Neighborhood_Volume(){}

std::vector<Field_NS::FieldSpec> PeridigmNS::Compute_Neighborhood_Volume::getFieldSpecs() const {
  std::vector<Field_NS::FieldSpec> myFieldSpecs;
  myFieldSpecs.push_back(Field_NS::NEIGHBORHOOD_VOLUME);
  return myFieldSpecs;
}

int PeridigmNS::Compute_Neighborhood_Volume::compute(Teuchos::RCP<PeridigmNS::DataManager> dataManager) const {

  // \todo Move to an initialization function, since this only needs to be computed once.

//   double *volume, *partialVolume, *neighborhoodVolume;
//   dataManager->getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&volume);
//   dataManager->getData(Field_NS::PARTIAL_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&partialVolume);
//   dataManager->getData(Field_NS::NEIGHBORHOOD_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&neighborhoodVolume);

//   Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = dataManager->getNeighborhoodData();
//   const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
//   const int* ownedIDs = neighborhoodData->OwnedIDs();
//   const int* neighborhoodList = neighborhoodData->NeighborhoodList();

//   int neighborhoodListIndex = 0;
//   int bondIndex = 0;
//   for(int iID=0 ; iID<numOwnedPoints ; ++iID){
// 	int nodeID = ownedIDs[iID];
// 	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
//     neighborhoodVolume[iID] = 0.0;
// 	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
// 	  int neighborID = neighborhoodList[neighborhoodListIndex++];
//       double neighborVolume = volume[neighborID];
//       double partialVolume = partialVolume[bondIndex++];
//       neighborhoodVolume[iID] += neighborVolume*partialVolume;
// 	}
//   }
  return 0;
}
