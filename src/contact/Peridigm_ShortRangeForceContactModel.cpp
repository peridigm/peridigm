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
#include <Teuchos_Assert.hpp>

PeridigmNS::ShortRangeForceContactModel::ShortRangeForceContactModel(const Teuchos::ParameterList& params)
  : ContactModel(params),
    m_contactRadius(0.0),
    m_springConstant(0.0),
    m_frictionCoefficient(0.0),
    m_horizon(0.0)
{
  //! \todo Add meaningful asserts on material properties.
  if(!params.isParameter("Contact Radius"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Contact Radius\" not specified.");
  m_contactRadius = params.get<double>("Contact Radius");
  if(!params.isParameter("Spring Constant"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Spring Constant\" not specified.");
  m_springConstant = params.get<double>("Spring Constant");
  if(!params.isParameter("Friction Coefficient"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Friction Coefficient\" not specified.");
  m_frictionCoefficient = params.get<double>("Friction Coefficient");
  if(!params.isParameter("Horizon"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Short range force contact parameter \"Horizon\" not specified.");
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

  double *cellVolume, *y, *contactForce, *velocity;
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::VELOC3D, Field_ENUM::STEP_NP1)->ExtractView(&velocity);
  dataManager.getData(Field_NS::CONTACT_FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->ExtractView(&contactForce);

  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    int numNeighbors = contactNeighborhoodList[neighborhoodListIndex++];
    if(numNeighbors > 0){
      int nodeID = ownedIDs[iID];
      double nodeCurrentX[3] = { y[nodeID*3],
                                 y[nodeID*3+1],
                                 y[nodeID*3+2] };
      double nodeCurrentV[3] = { velocity[nodeID*3],
                                 velocity[nodeID*3+1],
                                 velocity[nodeID*3+2] };
      double nodeVolume = cellVolume[nodeID];
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
        int neighborID = contactNeighborhoodList[neighborhoodListIndex++];
        TEUCHOS_TEST_FOR_EXCEPT_MSG(neighborID < 0, "Invalid neighbor list\n");
        double currentDistance =
          distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
                   y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
        if(currentDistance < m_contactRadius){
          double c = 9.0*m_springConstant/(3.1415*m_horizon*m_horizon*m_horizon*m_horizon);	// half value (of 18) due to force being applied to both nodes
          double temp = c*(m_contactRadius - currentDistance)/m_horizon;
          double neighborVolume = cellVolume[neighborID];
          
          if (m_frictionCoefficient != 0.0){
            // calculate the perpendicular velocity of the current node wrt the vector between the nodes 
            double l[3] = { y[neighborID*3] - nodeCurrentX[0], 
                            y[neighborID*3+1] - nodeCurrentX[1], 
                            y[neighborID*3+2] - nodeCurrentX[2] };
            double normal[3] = { l[0]/currentDistance,
                                 l[1]/currentDistance,
                                 l[2]/currentDistance };
            double currentDotNormal = nodeCurrentV[0]*normal[0] + 
                                      nodeCurrentV[1]*normal[1] + 
                                      nodeCurrentV[2]*normal[2];
            double currentDotNeighbor = velocity[neighborID*3]*normal[0] + 
                                        velocity[neighborID*3+1]*normal[1] + 
                                        velocity[neighborID*3+2]*normal[2]; 
            double nodeCurrentVperp[3] = { nodeCurrentV[0] - currentDotNormal*normal[0],
                                           nodeCurrentV[1] - currentDotNormal*normal[1],
                                           nodeCurrentV[2] - currentDotNormal*normal[2] };           
            double nodeNeighborVperp[3] = { velocity[neighborID*3] - currentDotNormal*normal[0],
                                            velocity[neighborID*3+1] - currentDotNormal*normal[1],
                                            velocity[neighborID*3+2] - currentDotNormal*normal[2] };

            // calculate frame of reference for the perpendicular velocities
            double Vcm[3] = { 0.5*(nodeCurrentVperp[0] + nodeNeighborVperp[0]),
                              0.5*(nodeCurrentVperp[1] + nodeNeighborVperp[1]),
                              0.5*(nodeCurrentVperp[2] + nodeNeighborVperp[2]) };
          
            // calculate the relative velocity of the current node wrt the neighboring node and vice versa
            double nodeCurrentVrel[3] = { nodeCurrentVperp[0] - Vcm[0],
                                          nodeCurrentVperp[1] - Vcm[1],
                                          nodeCurrentVperp[2] - Vcm[2] };
            double nodeNeighborVrel[3] = { nodeNeighborVperp[0] - Vcm[0],
                                           nodeNeighborVperp[1] - Vcm[1],
                                           nodeNeighborVperp[2] - Vcm[2] };
            double normCurrentVrel = sqrt(nodeCurrentVrel[0]*nodeCurrentVrel[0] + 
                                          nodeCurrentVrel[1]*nodeCurrentVrel[1] + 
                                          nodeCurrentVrel[2]*nodeCurrentVrel[2]); 
            double normNeighborVrel = sqrt(nodeNeighborVrel[0]*nodeNeighborVrel[0] + 
                                           nodeNeighborVrel[1]*nodeNeighborVrel[1] + 
                                           nodeNeighborVrel[2]*nodeNeighborVrel[2]);         
            
            // calculate the normal forces
            double currentNormalForce[3] = { -(temp*neighborVolume*(y[neighborID*3]   - nodeCurrentX[0])/currentDistance),
                                             -(temp*neighborVolume*(y[neighborID*3+1] - nodeCurrentX[1])/currentDistance),
                                             -(temp*neighborVolume*(y[neighborID*3+2] - nodeCurrentX[2])/currentDistance) };
            double neighborNormalForce[3] = { (temp*nodeVolume*(y[neighborID*3]   - nodeCurrentX[0])/currentDistance),
                                              (temp*nodeVolume*(y[neighborID*3+1] - nodeCurrentX[1])/currentDistance),
                                              (temp*nodeVolume*(y[neighborID*3+2] - nodeCurrentX[2])/currentDistance) };
            double normCurrentNormalForce = sqrt(currentNormalForce[0]*currentNormalForce[0] + 
                                                 currentNormalForce[1]*currentNormalForce[1] + 
                                                 currentNormalForce[2]*currentNormalForce[2]);
            double normNeighborNormalForce = sqrt(neighborNormalForce[0]*neighborNormalForce[0] + 
                                                  neighborNormalForce[1]*neighborNormalForce[1] + 
                                                  neighborNormalForce[2]*neighborNormalForce[2]);
            
            // calculate the friction forces
            double currentFrictionForce[3] = {0,0,0};
            double neighborFrictionForce[3] = {0,0,0};
            if (normCurrentVrel != 0.0) {
              currentFrictionForce[0] = -m_frictionCoefficient*normCurrentNormalForce*nodeCurrentVrel[0]/normCurrentVrel;
              currentFrictionForce[1] = -m_frictionCoefficient*normCurrentNormalForce*nodeCurrentVrel[1]/normCurrentVrel;
              currentFrictionForce[2] = -m_frictionCoefficient*normCurrentNormalForce*nodeCurrentVrel[2]/normCurrentVrel;
            }
            if (normNeighborVrel != 0.0) {
              neighborFrictionForce[0] = -m_frictionCoefficient*normNeighborNormalForce*nodeNeighborVrel[0]/normNeighborVrel;
              neighborFrictionForce[1] = -m_frictionCoefficient*normNeighborNormalForce*nodeNeighborVrel[1]/normNeighborVrel;
              neighborFrictionForce[2] = -m_frictionCoefficient*normNeighborNormalForce*nodeNeighborVrel[2]/normNeighborVrel;
            }

            // compute total contributions to force density
            contactForce[nodeID*3]       += currentNormalForce[0] + currentFrictionForce[0];
            contactForce[nodeID*3+1]     += currentNormalForce[1] + currentFrictionForce[1];
            contactForce[nodeID*3+2]     += currentNormalForce[2] + currentFrictionForce[2];
            contactForce[neighborID*3]   += neighborNormalForce[0] + neighborFrictionForce[0];
            contactForce[neighborID*3+1] += neighborNormalForce[1] + neighborFrictionForce[1];
            contactForce[neighborID*3+2] += neighborNormalForce[2] + neighborFrictionForce[2];
          }
          else {          
          // compute contributions to force density (Normal Force Only)
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
}
