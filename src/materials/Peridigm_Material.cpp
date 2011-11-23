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
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialComm.h>

using namespace std;

void PeridigmNS::Material::computeJacobian(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager,
                                           PeridigmNS::SerialMatrix& jacobian) const
{
  // Compute a finite-difference Jacobian using either FORWARD_DIFFERENCE or CENTRAL_DIFFERENCE
  computeFiniteDifferenceJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian, CENTRAL_DIFFERENCE);
}

void PeridigmNS::Material::computeFiniteDifferenceJacobian(const double dt,
                                                           const int numOwnedPoints,
                                                           const int* ownedIDs,
                                                           const int* neighborhoodList,
                                                           PeridigmNS::DataManager& dataManager,
                                                           PeridigmNS::SerialMatrix& jacobian,
                                                           FiniteDifferenceScheme finiteDifferenceScheme) const
{
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
  // Forward difference:
  // dF_0x/dx_0 = ( F_0x(perturbed x_0) - F_0x(unperturbed) ) / epsilon
  //
  // Central difference:
  // dF_0x/dx_0 = ( F_0x(positive perturbed x_0) - F_0x(negative perturbed x_0) ) / ( 2.0*epsilon )

  double epsilon = 1.0e-6; // \todo Instead, use 1.0e-6 * smallest_radius_in_model

  // loop over all points
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // create a temporary neighborhood consisting of a single point and its neighbors
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    vector<int> tempMyGlobalIDs(numNeighbors+1);
    tempMyGlobalIDs[0] = dataManager.getOverlapScalarPointMap()->GID(iID);
    vector<int> tempNeighborhoodList(numNeighbors+1); 
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // create a temporary DataManager containing data for this point and its neighborhood
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // the temporary data manager will have the same field specs and data as the real data manager
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs = dataManager.getFieldSpecs();
    tempDataManager.allocateData(fieldSpecs);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // set up numOwnedPoints and ownedIDs
    // there is only one owned ID, and it has local ID zero in the tempDataManager
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data in the constitutiveData array
    double *y, *force;
    tempDataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->ExtractView(&force);

    // Create a temporary vector for storing force
    Teuchos::RCP<Epetra_Vector> forceVector = tempDataManager.getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1);
    Teuchos::RCP<Epetra_Vector> tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
    double* tempForce;
    tempForceVector->ExtractView(&tempForce);

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < 3*(numNeighbors+1))
      scratchMatrix.Resize(3*(numNeighbors+1));

    // Create a list of indices for the rows/columns in the scratch matrix
    // These indices correspond to the DataManager's three-dimensional overlap map
    vector<int> indices(3*(numNeighbors+1));
    for(int i=0 ; i<numNeighbors+1 ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      int localID = dataManager.getOverlapScalarPointMap()->LID(globalID);
      for(int j=0 ; j<3 ; ++j)
        indices[3*i+j] = 3*localID+j;
    }

    if(finiteDifferenceScheme == FORWARD_DIFFERENCE){
      // compute and store the unperturbed force
      updateConstitutiveData(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      for(int i=0 ; i<forceVector->MyLength() ; ++i)
        tempForce[i] = force[i];
    }

    // perturb one dof in the neighborhood at a time and compute the force
    // the point itself plus each of its neighbors must be perturbed
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = 0;

      for(int dof=0 ; dof<3 ; ++dof){

        // perturb a dof and compute the forces
        double oldY = y[3*perturbID+dof];

        if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
          // compute and store the negatively perturbed force
          y[3*perturbID+dof] -= epsilon;
          updateConstitutiveData(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          y[3*perturbID+dof] = oldY;
          for(int i=0 ; i<forceVector->MyLength() ; ++i)
            tempForce[i] = force[i];
        }

        // compute the positively perturbed force
        y[3*perturbID+dof] += epsilon;
        updateConstitutiveData(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        y[3*perturbID+dof] = oldY;

        for(int i=0 ; i<numNeighbors+1 ; ++i){

          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = 0;

          for(int d=0 ; d<3 ; ++d){
            double value = ( force[3*forceID+d] - tempForce[3*forceID+d] ) / epsilon;
            if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
              value *= 0.5;
            scratchMatrix(3*forceID+d, 3*perturbID+dof) = value;
          }
        }
      }
    }

    // sum the values into the global tangent matrix (this is expensive)
    jacobian.addValues((int)indices.size(), &indices[0], scratchMatrix.Data());
  }
}
