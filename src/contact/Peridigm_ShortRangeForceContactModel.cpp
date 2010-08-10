/*! \file Peridigm_ShortRangeForceContactModel.cpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#include "Peridigm_ShortRangeForceContactModel.hpp"
#include <Teuchos_TestForException.hpp>

Peridigm::ShortRangeForceContactModel::ShortRangeForceContactModel(const Teuchos::ParameterList& params)
  : ContactModel(params),
    m_decompStates(),
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

Peridigm::ShortRangeForceContactModel::~ShortRangeForceContactModel()
{
}

void
Peridigm::ShortRangeForceContactModel::computeForce(const Epetra_Vector& x,
                                                    const Epetra_Vector& u,
                                                    const Epetra_Vector& v,
                                                    const double dt,
                                                    const Epetra_Vector& cellVolume,
                                                    const int numOwnedPoints,
                                                    const int* ownedIDs,
                                                    const int* contactNeighborhoodList,
                                                    Epetra_MultiVector& scalarConstitutiveData,
                                                    Epetra_MultiVector& vectorConstitutiveData,
                                                    Epetra_Vector& force) const
{
  std::pair<int,double*> vectorView = m_decompStates.extractStrideView(vectorConstitutiveData);
  double *y = m_decompStates.extractCurrentPositionView(vectorView);

  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int numNeighbors = contactNeighborhoodList[neighborhoodListIndex++];
    if(numNeighbors > 0){
      int nodeID = ownedIDs[iID];
      TEST_FOR_EXCEPT_MSG(nodeID*3+2 >= force.MyLength(), "Invalid neighbor list / contact force vector\n");
      double nodeCurrentX[3] = { y[nodeID*3],
                                 y[nodeID*3+1],
                                 y[nodeID*3+2] };
      double nodeVolume = cellVolume[nodeID];
      for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
        int neighborID = contactNeighborhoodList[neighborhoodListIndex++];
        TEST_FOR_EXCEPT_MSG(neighborID < 0, "Invalid neighbor list\n");
        TEST_FOR_EXCEPT_MSG(neighborID*3+2 >= force.MyLength(), "Invalid neighbor list / contact force vector\n");
        double currentDistance =
          distance(nodeCurrentX[0], nodeCurrentX[1], nodeCurrentX[2],
                   y[neighborID*3], y[neighborID*3+1], y[neighborID*3+2]);
        if(currentDistance < m_contactRadius){
          double c = 9.0*m_springConstant/(3.1415*m_horizon*m_horizon*m_horizon*m_horizon);
          double temp = c*(m_contactRadius - currentDistance)/m_horizon;
          double neighborVolume = cellVolume[neighborID];

          // compute contributions to force density
            
          force[nodeID*3]       -= temp*neighborVolume*(y[neighborID*3]   - nodeCurrentX[0])/currentDistance;
          force[nodeID*3+1]     -= temp*neighborVolume*(y[neighborID*3+1] - nodeCurrentX[1])/currentDistance;
          force[nodeID*3+2]     -= temp*neighborVolume*(y[neighborID*3+2] - nodeCurrentX[2])/currentDistance;
          force[neighborID*3]   += temp*nodeVolume*(y[neighborID*3]   - nodeCurrentX[0])/currentDistance;
          force[neighborID*3+1] += temp*nodeVolume*(y[neighborID*3+1] - nodeCurrentX[1])/currentDistance;
          force[neighborID*3+2] += temp*nodeVolume*(y[neighborID*3+2] - nodeCurrentX[2])/currentDistance;
        }
      }
	}
  }
}
