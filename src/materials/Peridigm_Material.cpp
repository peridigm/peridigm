/*! \file Peridigm_Material.cpp */

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

#include "Peridigm_Material.hpp"
#include <Teuchos_TestForException.hpp>

using namespace std;

void PeridigmNS::Material::computeJacobian(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager,
                                           PeridigmNS::SerialMatrix& jacobian) const
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

  double epsilon = 1.0e-6; // \todo Instead, use 1.0e-6 * smallest_radius_in_model

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
      int neighborID = neighborhoodList[neighborhoodListIndex++];
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
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = iID;

      for(int dof=0 ; dof<3 ; ++dof){

        // perturb a dof and compute the forces
        double oldY = y[3*perturbID+dof];
        y[3*perturbID+dof] += epsilon;
        updateConstitutiveData(dt, tempNumOwnedPoints, tempOwnedIDs, &tempNeighborhoodList[0], dataManager);
        computeForce(dt, tempNumOwnedPoints, tempOwnedIDs, &tempNeighborhoodList[0], dataManager);

        // fill in the corresponding Jacobian entries
        for(int i=0 ; i<numNeighbors+1 ; ++i){

          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = iID;

          for(int d=0 ; d<3 ; ++d){
            double value = ( force[3*forceID+d] - unperturbedForce[3*forceID+d] ) / epsilon;
            int row = 3*forceID + d;
            int col = 3*perturbID + dof;
            jacobian.addValue(row, col, value);
          }
        }

        // unperturb the dof
        y[3*perturbID+dof] = oldY;
      }
    }
  }
}
