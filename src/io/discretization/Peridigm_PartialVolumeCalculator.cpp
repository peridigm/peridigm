/*! \file Peridigm_PartialVolumeCalculator.cpp */

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

#include "Peridigm_PartialVolumeCalculator.hpp"

void PeridigmNS::computePartialVolume(Teuchos::RCP<PeridigmNS::Block> block,
                                      Teuchos::RCP<PeridigmNS::AbstractDiscretization> discretization)
{
  double *x, *partialVolume;
  block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
  block->getData(Field_NS::PARTIAL_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&partialVolume);

  int vectorLength = block->getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->MyLength();

  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = block->getNeighborhoodData();
  const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
  const int* ownedIDs = neighborhoodData->OwnedIDs();
  const int* neighborhoodList = neighborhoodData->NeighborhoodList();
  
  double horizon = discretization->getHorizon();

  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int nodeID = ownedIDs[iID];
	TEST_FOR_EXCEPT_MSG(nodeID*3+2 >= vectorLength, "Invalid neighbor list / x vector\n");
	double nodeInitialX[3] = { x[nodeID*3],
							   x[nodeID*3+1],
							   x[nodeID*3+2] };
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
	  int neighborLocalID = neighborhoodList[neighborhoodListIndex++];
	  TEST_FOR_EXCEPT_MSG(neighborLocalID < 0, "Invalid neighbor list\n");
	  TEST_FOR_EXCEPT_MSG(neighborLocalID*3+2 >= vectorLength, "Invalid neighbor list / initial x vector\n");
      int neighborGlobalID = block->getOverlapScalarPointMap()->GID(neighborLocalID);

      Teuchos::RCP< std::vector<double> > exodusNodePositions =
        discretization->getExodusMeshNodePositions(neighborGlobalID);

      if(!exodusNodePositions.is_null())
        partialVolume[bondIndex] = computePartialVolume(nodeInitialX, &(*exodusNodePositions)[0], horizon);

	  bondIndex += 1;
	}
  }
}

double PeridigmNS::computePartialVolume(const double* const pt,
                                        const double* const neighborElementNodes,
                                        const double horizon)
{
  return 1.0;
}
