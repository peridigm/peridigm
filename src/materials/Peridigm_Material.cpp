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
#include <Epetra_SerialComm.h>

using namespace std;

void PeridigmNS::Material::computeJacobian(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager,
                                           PeridigmNS::SerialMatrix& jacobian) const
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
  // dF_0x/dx_0 = ( F_0x(perturbed x_0) - F_0x(unperturbed) ) / epsilon

  double epsilon = 1.0e-6; // \todo Instead, use 1.0e-6 * smallest_radius_in_model

  // loop over all points
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // create a temporary neighborhood consisting of a single point and its neighbors
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    vector<int> tempMyGlobalIDs(numNeighbors+1);
    tempMyGlobalIDs[0] = dataManager.getOverlapIDScalarMap()->GID(iID);
    vector<int> tempNeighborhoodList(numNeighbors+1); 
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapIDScalarMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // create a temporary DataManager containing data for this point and its neighborhood
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            tempThreeDimensionalMap,
                            tempBondMap);

    // the temporary data manager will have the same field specs as the real data manager
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs = dataManager.getFieldSpecs();
    tempDataManager.allocateData(fieldSpecs);

    // copy data into the tempDataManager
    // \todo Can this be done with an import (parallel comm -> serial comm)?
    // \todo Write a DataManager function to perform this copy operation.
    for(unsigned int iFieldSpec=0 ; iFieldSpec<fieldSpecs->size() ; ++iFieldSpec){

      Field_NS::FieldSpec spec = (*fieldSpecs)[iFieldSpec];
      Field_NS::FieldSpec::FieldStateArchitecture arch = spec.getStateArchitecture();
      Field_NS::FieldSpec::FieldLength length = spec.getLength();

      // for scalar and vector data, copy data for the given point and all its neighbors
      if(length == Field_NS::FieldSpec::SCALAR || length == Field_NS::FieldSpec::VECTOR3D){

        if(arch == Field_NS::FieldSpec::STATELESS){
          Epetra_Vector& source = *(dataManager.getData(spec, Field_NS::FieldSpec::STEP_NONE));
          Epetra_Vector& target = *(tempDataManager.getData(spec, Field_NS::FieldSpec::STEP_NONE));
          for(unsigned int i=0 ; i<tempMyGlobalIDs.size() ; ++i){
            int globalID = tempMyGlobalIDs[i];
            int sourceLocalID = source.Map().LID(globalID);
            int targetLocalID = target.Map().LID(globalID);
            if(spec.getLength() == Field_NS::FieldSpec::SCALAR){
              target[targetLocalID] = source[sourceLocalID];
            }
            else if(spec.getLength() == Field_NS::FieldSpec::VECTOR3D){
              target[3*targetLocalID] = source[3*sourceLocalID];
              target[3*targetLocalID+1] = source[3*sourceLocalID+1];
              target[3*targetLocalID+2] = source[3*sourceLocalID+2];
            }
          }
        }
        else if(arch == Field_NS::FieldSpec::STATEFUL){
          {
            Epetra_Vector& source = *(dataManager.getData(spec, Field_NS::FieldSpec::STEP_N));
            Epetra_Vector& target = *(tempDataManager.getData(spec, Field_NS::FieldSpec::STEP_N));
            for(unsigned int i=0 ; i<tempMyGlobalIDs.size() ; ++i){
              int globalID = tempMyGlobalIDs[i];
              int sourceLocalID = source.Map().LID(globalID);
              int targetLocalID = target.Map().LID(globalID);
              if(spec.getLength() == Field_NS::FieldSpec::SCALAR){
                target[targetLocalID] = source[sourceLocalID];
              }
              else if(spec.getLength() == Field_NS::FieldSpec::VECTOR3D){
                target[3*targetLocalID] = source[3*sourceLocalID];
                target[3*targetLocalID+1] = source[3*sourceLocalID+1];
                target[3*targetLocalID+2] = source[3*sourceLocalID+2];
              }
            }
          }
          {
            Epetra_Vector& source = *(dataManager.getData(spec, Field_NS::FieldSpec::STEP_NP1));
            Epetra_Vector& target = *(tempDataManager.getData(spec, Field_NS::FieldSpec::STEP_NP1));
            for(unsigned int i=0 ; i<tempMyGlobalIDs.size() ; ++i){
              int globalID = tempMyGlobalIDs[i];
              int sourceLocalID = source.Map().LID(globalID);
              int targetLocalID = target.Map().LID(globalID);
              if(spec.getLength() == Field_NS::FieldSpec::SCALAR){
                target[targetLocalID] = source[sourceLocalID];
              }
              else if(spec.getLength() == Field_NS::FieldSpec::VECTOR3D){
                target[3*targetLocalID] = source[3*sourceLocalID];
                target[3*targetLocalID+1] = source[3*sourceLocalID+1];
                target[3*targetLocalID+2] = source[3*sourceLocalID+2];
              }
            }
          }
        }
      }
      else if(spec.getLength() == Field_NS::FieldSpec::BOND){

        if(arch == Field_NS::FieldSpec::STATELESS){
          Epetra_Vector& source = *(dataManager.getData(spec, Field_NS::FieldSpec::STEP_NONE));
          Epetra_Vector& target = *(tempDataManager.getData(spec, Field_NS::FieldSpec::STEP_NONE));
          // there is bond data only for the owned ID
          int globalID = tempMyGlobalIDs[0];
          int sourceLocalID = source.Map().LID(globalID);
          int targetLocalID = target.Map().LID(globalID);
          int sourceFirstPointInElement = source.Map().FirstPointInElement(sourceLocalID);
          int targetFirstPointInElement = target.Map().FirstPointInElement(targetLocalID);
          int elementSize = source.Map().ElementSize(sourceLocalID);
          for(int i=0 ; i<elementSize ; ++i)
            target[targetFirstPointInElement+i] = source[sourceFirstPointInElement+i];
        }
        else if(arch == Field_NS::FieldSpec::STATEFUL){
          Epetra_Vector& source = *(dataManager.getData(spec, Field_NS::FieldSpec::STEP_N));
          Epetra_Vector& target = *(tempDataManager.getData(spec, Field_NS::FieldSpec::STEP_N));
          // there is bond data only for the owned ID
          int globalID = tempMyGlobalIDs[0];
          int sourceLocalID = source.Map().LID(globalID);
          int targetLocalID = target.Map().LID(globalID);
          int sourceFirstPointInElement = source.Map().FirstPointInElement(sourceLocalID);
          int targetFirstPointInElement = target.Map().FirstPointInElement(targetLocalID);
          int elementSize = source.Map().ElementSize(sourceLocalID);
          for(int i=0 ; i<elementSize ; ++i)
            target[targetFirstPointInElement+i] = source[sourceFirstPointInElement+i];
          source = *(dataManager.getData(spec, Field_NS::FieldSpec::STEP_NP1));
          target = *(tempDataManager.getData(spec, Field_NS::FieldSpec::STEP_NP1));
          globalID = tempMyGlobalIDs[0];
          sourceLocalID = source.Map().LID(globalID);
          targetLocalID = target.Map().LID(globalID);
          sourceFirstPointInElement = source.Map().FirstPointInElement(sourceLocalID);
          targetFirstPointInElement = target.Map().FirstPointInElement(targetLocalID);
          elementSize = source.Map().ElementSize(sourceLocalID);
          for(int i=0 ; i<elementSize ; ++i)
            target[targetFirstPointInElement+i] = source[sourceFirstPointInElement+i];
        }
      }
    }

    // set up numOwnedPoints and ownedIDs
    // there is only one owned ID, and it has local ID zero in the tempDataManager
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data in the constitutiveData array
    double *y, *force;
    tempDataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&force);

    // Create a temporary vector for storing the unperturbed force
    Teuchos::RCP<Epetra_Vector> forceVector = tempDataManager.getData(Field_NS::FORCE_DENSITY3D, Field_NS::FieldSpec::STEP_NP1);
    Teuchos::RCP<Epetra_Vector> unperturbedForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
    double* unperturbedForce;
    unperturbedForceVector->ExtractView(&unperturbedForce);

    // compute and store the unperturbed force
    updateConstitutiveData(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
    computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
    for(int i=0 ; i<forceVector->MyLength() ; ++i)
      unperturbedForce[i] = force[i];

    // perturb one dof in the neighborhood at a time and compute the force
    // the point itself plus each of its neighbors must be perturbed
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = 0;//iID;

      for(int dof=0 ; dof<3 ; ++dof){

        // perturb a dof and compute the forces
        double oldY = y[3*perturbID+dof];
        y[3*perturbID+dof] += epsilon;
        updateConstitutiveData(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);

        // fill in the corresponding Jacobian entries
        for(int i=0 ; i<numNeighbors+1 ; ++i){

          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = 0;//iID;

          for(int d=0 ; d<3 ; ++d){
            double value = ( force[3*forceID+d] - unperturbedForce[3*forceID+d] ) / epsilon;
            int globalForceID = tempOneDimensionalMap->GID(forceID);
            int globalPerturbID = tempOneDimensionalMap->GID(perturbID);
            int localForceID = dataManager.getOwnedIDScalarMap()->LID(globalForceID);
            int localPerturbID = dataManager.getOwnedIDScalarMap()->LID(globalPerturbID);
            int row = 3*localForceID + d;
            int col = 3*localPerturbID + dof;
            jacobian.addValue(row, col, value);
          }
        }

        // unperturb the dof
        y[3*perturbID+dof] = oldY;
      }
    }
  }
}
