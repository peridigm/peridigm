/*! \file Peridigm_Material.cpp */

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

#include "Peridigm_Material.hpp"
#include <Teuchos_TestForException.hpp>

using namespace std;

void PeridigmNS::Material::computeJacobian(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager,
                                           PeridigmNS::DenseMatrix& jacobian) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *y, *force;
  dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&force);

  // The Jacobian is of the form:
  //
  // dF_0x/dx_0  dF_0x/dy_0  dF_0x/dz_0  dF_0x/dx_1  dF_0x/dy_1  dF_0x/dz_1  ...  dF_0x/dx_n  dF_0x/dy_n  dF_0x/dz_n  
  // dF_0y/dx_0  dF_0y/dy_0  dF_0y/dz_0  dF_0y/dx_1  dF_0y/dy_1  dF_0y/dz_1  ...  dF_0y/dx_n  dF_0y/dy_n  dF_0y/dz_n  
  // dF_0z/dx_0  dF_0z/dy_0  dF_0z/dz_0  dF_0z/dx_1  dF_0z/dy_1  dF_0z/dz_1  ...  dF_0z/dx_n  dF_0z/dy_n  dF_0z/dz_n  
  // dF_1x/dx_0  dF_1x/dy_0  dF_1x/dz_0  dF_1x/dx_1  dF_1x/dy_1  dF_1x/dz_1  ...  dF_1x/dx_n  dF_1x/dy_n  dF_1x/dz_n  
  // dF_1y/dx_0  dF_1y/dy_0  dF_1y/dz_0  dF_1y/dx_1  dF_1y/dy_1  dF_1y/dz_1  ...  dF_1y/dx_n  dF_1y/dy_n  dF_1y/dz_n  
  // dF_1z/dx_0  dF_1z/dy_0  dF_1z/dz_0  dF_1z/dx_1  dF_1z/dy_1  dF_1z/dz_1  ...  dF_1z/dx_n  dF_1z/dy_n  dF_1z/dz_n  
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  //     .           .           .           .           .           .                .           .           .
  // dF_nx/dx_0  dF_nx/dy_0  dF_nx/dz_0  dF_nx/dx_1  dF_nx/dy_1  dF_nx/dz_1  ...  dF_nx/dx_n  dF_nx/dy_n  dF_nx/dz_n  
  // dF_ny/dx_0  dF_ny/dy_0  dF_ny/dz_0  dF_ny/dx_1  dF_ny/dy_1  dF_ny/dz_1  ...  dF_ny/dx_n  dF_ny/dy_n  dF_ny/dz_n  
  // dF_nz/dx_0  dF_nz/dy_0  dF_nz/dz_0  dF_nz/dx_1  dF_nz/dy_1  dF_nz/dz_1  ...  dF_nz/dx_n  dF_nz/dy_n  dF_nz/dz_n  

  // Each entry is computed by finite difference:
  //
  // dF_0x/dx_0 = ( F_0x(perturbed x_0) - F_0x(unperturbed) ) / epsilon

  double epsilon = 1.0e-8; // \todo Instead, use 1.0e-6 * smallest_radius_in_model

  // Create a temporary vector for storing the unperturbed force
  Teuchos::RCP<Epetra_Vector> forceVector = dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1);
  Teuchos::RCP<Epetra_Vector> unperturbedForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
  double* unperturbedForce;
  unperturbedForceVector->ExtractView(&unperturbedForce);

  // loop over all points
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // create a neighborhood consisting of a single point and its neighbors
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    int tempNumOwnedPoints = 1;
    int* tempOwnedIDs = &iID;
    vector<int> tempNeighborhoodList(numNeighbors+1); 
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex + iNID];
      tempNeighborhoodList[iNID+1] = neighborID;
    }

    // compute and store the unperturbed force
    updateConstitutiveData(dt, tempNumOwnedPoints, tempOwnedIDs, &tempNeighborhoodList[0], dataManager);
    computeForce(dt, tempNumOwnedPoints, tempOwnedIDs, &tempNeighborhoodList[0], dataManager);
    for(int i=0 ; i<forceVector->MyLength() ; ++i)
      unperturbedForce[i] = force[i];

    // perturb one dof in the neighborhood at a time and compute the force
    // the point itself plus each of its neighbors must be perturbed
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID];
      else
        perturbID = iID;

      for(int dof=0 ; dof<3 ; ++dof){

        // perturb a dof and compute the forces
        // \todo This will work but may be more expensive than it needs to be, may be able to compute unperturbed force only once.
        y[3*perturbID+dof] += 0.0; // += epsilon;
        updateConstitutiveData(dt, tempNumOwnedPoints, tempOwnedIDs, &tempNeighborhoodList[0], dataManager);
        computeForce(dt, tempNumOwnedPoints, tempOwnedIDs, &tempNeighborhoodList[0], dataManager);

        // fill in the corresponding Jacobian entries
        for(int i=0 ; i<numNeighbors+1 ; ++i){

          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i];
          else
            forceID = iID;

          for(int j=0 ; j<3 ; ++j){
            double value = ( force[3*forceID+j] - unperturbedForce[3*forceID+j] ) / epsilon;
            int row = i*forceID + j;
            int col = 3*perturbID + dof;
            jacobian.addValue(row, col, value);
          }
        }

        // unperturb the dof
        y[3*perturbID+dof] -= 0.0; // -= epsilon;
      }
    }

    neighborhoodListIndex += numNeighbors;
  }
}
