/*! \file Peridigm_ShortRangeForceContactModel.cpp */

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

#include "Peridigm_ShortRangeForceContactModel.hpp"
#include <Teuchos_TestForException.hpp>

PeridigmNS::ShortRangeForceContactModel::ShortRangeForceContactModel(const Teuchos::ParameterList& params)
  : ContactModel(params),
    m_contactRadius(0.0),
    m_springConstant(0.0),
    m_horizon(0.0)
{
  //! \todo Add meaningful asserts on material properties.
  if(!params.isParameter("Contact Radius"))
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Contact Radius\" not specified.");
  m_contactRadius = params.get<double>("Contact Radius");
  if(!params.isParameter("Spring Constant"))
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Spring Constant\" not specified.");
  m_springConstant = params.get<double>("Spring Constant");
  if(!params.isParameter("Horizon"))
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Horizon\" not specified.");
  m_horizon = params.get<double>("Horizon");
}

PeridigmNS::ShortRangeForceContactModel::~ShortRangeForceContactModel()
{
}

void
PeridigmNS::ShortRangeForceContactModel::computeForce(const double dt,
                                                      const int numOwnedPoints,
                                                      const int* ownedIDs,
                                                      const int* contactNeighborhoodList,
                                                      PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(Field_NS::CONTACT_FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->PutScalar(0.0);

  double *cellVolume, *y, *contactForce;
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::CONTACT_FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->ExtractView(&contactForce);

  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int numNeighbors = contactNeighborhoodList[neighborhoodListIndex++];
    if(numNeighbors > 0){
      int nodeID = ownedIDs[iID];
      double nodeCurrentX[3] = { y[nodeID*3],
                                 y[nodeID*3+1],
                                 y[nodeID*3+2] };
      double nodeVolume = cellVolume[nodeID];
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
        int neighborID = contactNeighborhoodList[neighborhoodListIndex++];
        TEST_FOR_EXCEPT_MSG(neighborID < 0, "Invalid neighbor list\n");
        double currentDistance =
          distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
                   y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
        if(currentDistance < m_contactRadius){
          double c = 9.0*m_springConstant/(3.1415*m_horizon*m_horizon*m_horizon*m_horizon);
          double temp = c*(m_contactRadius - currentDistance)/m_horizon;
          double neighborVolume = cellVolume[neighborID];

          // compute contributions to force density
          contactForce[nodeID*3]       -= temp*neighborVolume*(y[neighborID*3]   - nodeCurrentX[0])/currentDistance;
          contactForce[nodeID*3+1]     -= temp*neighborVolume*(y[neighborID*3+1] - nodeCurrentX[1])/currentDistance;
          contactForce[nodeID*3+2]     -= temp*neighborVolume*(y[neighborID*3+2] - nodeCurrentX[2])/currentDistance;
          contactForce[neighborID*3]   += temp*nodeVolume*(y[neighborID*3]   - nodeCurrentX[0])/currentDistance;
          contactForce[neighborID*3+1] += temp*nodeVolume*(y[neighborID*3+1] - nodeCurrentX[1])/currentDistance;
          contactForce[neighborID*3+2] += temp*nodeVolume*(y[neighborID*3+2] - nodeCurrentX[2])/currentDistance;
        }
      }
	}
  }
}
