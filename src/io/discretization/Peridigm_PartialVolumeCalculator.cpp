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
using namespace std;

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
  int numSubdivisions = 16;
  double upperBound = 1.0;
  double lowerBound = -1.0;
  double stepSize = (upperBound - lowerBound)/numSubdivisions;

  double naturalX, naturalY, naturalZ, x, y, z, distance;
  vector<double> shape(8), dShape_dNaturalX(8), dShape_dNaturalY(8), dShape_dNaturalZ(8);
  double jacobian[3][3], jacobianDeterminate;

  double referenceVolume = 8.0/(numSubdivisions*numSubdivisions*numSubdivisions);
  double subVolume;
  double volume = 0.0;
  double volumeInNeighborhood = 0.0;

  // First, check to see if the neighbor's nodes are either all
  // within the neighborhood, or all outside the neighborhood.
  bool allIn = true;
  bool allOut = true;
  for(int i=0 ; i<8 ; ++i){
    x = neighborElementNodes[3*i];
    y = neighborElementNodes[3*i+1];
    z = neighborElementNodes[3*i+2];
    distance = (pt[0] - x)*(pt[0] - x) + (pt[1] - y)*(pt[1] - y) + (pt[2] - z)*(pt[2] - z);
    if(distance > horizon*horizon)
      allIn = false;
    else
      allOut = false;
  }
  if(allIn)
    return 1.0;
  else if(allOut)
    return 0.0;

  // Loop over a grid on the reference element
  for(int k=0 ; k<numSubdivisions ; ++k){
    for(int j=0 ; j<numSubdivisions ; ++j){
      for(int i=0 ; i<numSubdivisions ; ++i){

        naturalX = lowerBound + stepSize/2.0 + i*stepSize;
        naturalY = lowerBound + stepSize/2.0 + j*stepSize;
        naturalZ = lowerBound + stepSize/2.0 + k*stepSize;

        shape[0] = 0.125*(1.0 - naturalX)*(1.0 - naturalY)*(1.0 - naturalZ);
        shape[1] = 0.125*(1.0 + naturalX)*(1.0 - naturalY)*(1.0 - naturalZ);
        shape[2] = 0.125*(1.0 + naturalX)*(1.0 + naturalY)*(1.0 - naturalZ);
        shape[3] = 0.125*(1.0 - naturalX)*(1.0 + naturalY)*(1.0 - naturalZ);
        shape[4] = 0.125*(1.0 - naturalX)*(1.0 - naturalY)*(1.0 + naturalZ);
        shape[5] = 0.125*(1.0 + naturalX)*(1.0 - naturalY)*(1.0 + naturalZ);
        shape[6] = 0.125*(1.0 + naturalX)*(1.0 + naturalY)*(1.0 + naturalZ);
        shape[7] = 0.125*(1.0 - naturalX)*(1.0 + naturalY)*(1.0 + naturalZ);

        x = 0.0;
        y = 0.0;
        z = 0.0;
        for(int m=0 ; m<8 ; ++m){
          x += neighborElementNodes[3*m  ]*shape[m];
          y += neighborElementNodes[3*m+1]*shape[m];
          z += neighborElementNodes[3*m+2]*shape[m];
        }

        dShape_dNaturalX[0] = -0.125*(1.0 - naturalY)*(1.0 - naturalZ);
        dShape_dNaturalX[1] =  0.125*(1.0 - naturalY)*(1.0 - naturalZ);
        dShape_dNaturalX[2] =  0.125*(1.0 + naturalY)*(1.0 - naturalZ);
        dShape_dNaturalX[3] = -0.125*(1.0 + naturalY)*(1.0 - naturalZ);
        dShape_dNaturalX[4] = -0.125*(1.0 - naturalY)*(1.0 + naturalZ);
        dShape_dNaturalX[5] =  0.125*(1.0 - naturalY)*(1.0 + naturalZ);
        dShape_dNaturalX[6] =  0.125*(1.0 + naturalY)*(1.0 + naturalZ);
        dShape_dNaturalX[7] = -0.125*(1.0 + naturalY)*(1.0 + naturalZ);

        dShape_dNaturalY[0] = -0.125*(1.0 - naturalX)*(1.0 - naturalZ);
        dShape_dNaturalY[1] = -0.125*(1.0 + naturalX)*(1.0 - naturalZ);
        dShape_dNaturalY[2] =  0.125*(1.0 + naturalX)*(1.0 - naturalZ);
        dShape_dNaturalY[3] =  0.125*(1.0 - naturalX)*(1.0 - naturalZ);
        dShape_dNaturalY[4] = -0.125*(1.0 - naturalX)*(1.0 + naturalZ);
        dShape_dNaturalY[5] = -0.125*(1.0 + naturalX)*(1.0 + naturalZ);
        dShape_dNaturalY[6] =  0.125*(1.0 + naturalX)*(1.0 + naturalZ);
        dShape_dNaturalY[7] =  0.125*(1.0 - naturalX)*(1.0 + naturalZ);

        dShape_dNaturalZ[0] = -0.125*(1.0 - naturalX)*(1.0 - naturalY);
        dShape_dNaturalZ[1] = -0.125*(1.0 + naturalX)*(1.0 - naturalY);
        dShape_dNaturalZ[2] = -0.125*(1.0 + naturalX)*(1.0 + naturalY);
        dShape_dNaturalZ[3] = -0.125*(1.0 - naturalX)*(1.0 + naturalY);
        dShape_dNaturalZ[4] =  0.125*(1.0 - naturalX)*(1.0 - naturalY);
        dShape_dNaturalZ[5] =  0.125*(1.0 + naturalX)*(1.0 - naturalY);
        dShape_dNaturalZ[6] =  0.125*(1.0 + naturalX)*(1.0 + naturalY);
        dShape_dNaturalZ[7] =  0.125*(1.0 - naturalX)*(1.0 + naturalY);

        for(int m=0 ; m<3 ; ++m){
          for(int n=0 ; n<3 ; ++n){
            jacobian[m][n] = 0.0;
          }
        }
        for(int m=0 ; m<8 ; ++m){
          jacobian[0][0] += neighborElementNodes[3*m  ]*dShape_dNaturalX[m];
          jacobian[0][1] += neighborElementNodes[3*m+1]*dShape_dNaturalX[m];
          jacobian[0][2] += neighborElementNodes[3*m+2]*dShape_dNaturalX[m];
          jacobian[1][0] += neighborElementNodes[3*m  ]*dShape_dNaturalY[m];
          jacobian[1][1] += neighborElementNodes[3*m+1]*dShape_dNaturalY[m];
          jacobian[1][2] += neighborElementNodes[3*m+2]*dShape_dNaturalY[m];
          jacobian[2][0] += neighborElementNodes[3*m  ]*dShape_dNaturalZ[m];
          jacobian[2][1] += neighborElementNodes[3*m+1]*dShape_dNaturalZ[m];
          jacobian[2][2] += neighborElementNodes[3*m+2]*dShape_dNaturalZ[m];
        }

        jacobianDeterminate =
          jacobian[0][0]*(jacobian[1][1]*jacobian[2][2] - jacobian[1][2]*jacobian[2][1]) -
          jacobian[0][1]*(jacobian[1][0]*jacobian[2][2] - jacobian[1][2]*jacobian[2][0]) +
          jacobian[0][2]*(jacobian[1][0]*jacobian[2][1] - jacobian[1][1]*jacobian[2][0]);

        subVolume = referenceVolume*jacobianDeterminate;

        volume += subVolume;

        distance = (pt[0] - x)*(pt[0] - x) + (pt[1] - y)*(pt[1] - y) + (pt[2] - z)*(pt[2] - z);
        if(distance <= horizon*horizon)
          volumeInNeighborhood += subVolume;
      }
    }
  }
  double volumeFraction = volumeInNeighborhood/volume;

  return volumeFraction;
}
