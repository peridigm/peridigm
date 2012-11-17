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
#include "Peridigm_Field.hpp"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialComm.h>
#include <boost/math/special_functions/fpclassify.hpp>

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


  TEUCHOS_TEST_FOR_EXCEPT_MSG(m_finiteDifferenceProbeLength == DBL_MAX, "**** Finite-difference Jacobian requires that the \"Finite Difference Probe Length\" parameter be set.\n");
  double epsilon = m_finiteDifferenceProbeLength;

  // Get field ids for all relevant data
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId = fieldManager.getFieldId("Volume");
  int coordinatesFId = fieldManager.getFieldId("Coordinates");
  int tangentReferenceCoordinatesFId = fieldManager.getFieldId("Tangent_Reference_Coordinates");
  int forceDensityFId = fieldManager.getFieldId("Force_Density");

  // Loop over all points.
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    // The temporary neighborhood is sorted by global ID to somewhat increase the chance
    // that the downstream Epetra_CrsMatrix::SumIntoMyValues() calls will be as efficient
    // as possible.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    vector<int> tempMyGlobalIDs(numNeighbors+1);
    // Create a placeholder that will appear at the beginning of the sorted list.
    tempMyGlobalIDs[0] = -1;
    vector<int> tempNeighborhoodList(numNeighbors+1); 
    tempNeighborhoodList[0] = numNeighbors;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      tempMyGlobalIDs[iNID+1] = dataManager.getOverlapScalarPointMap()->GID(neighborID);
      tempNeighborhoodList[iNID+1] = iNID+1;
    }
    sort(tempMyGlobalIDs.begin(), tempMyGlobalIDs.end());
    // Put the node at the center of the neighborhood at the beginning of the list.
    tempMyGlobalIDs[0] = dataManager.getOwnedScalarPointMap()->GID(iID);

    Epetra_SerialComm serialComm;
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numNeighbors+1, numNeighbors+1, &tempMyGlobalIDs[0], 3, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempBondMap = Teuchos::rcp(new Epetra_BlockMap(1, 1, &tempMyGlobalIDs[0], numNeighbors, 0, serialComm));

    // Create a temporary DataManager containing data for this point and its neighborhood.
    PeridigmNS::DataManager tempDataManager;
    tempDataManager.setMaps(Teuchos::RCP<const Epetra_BlockMap>(),
                            tempOneDimensionalMap,
                            Teuchos::RCP<const Epetra_BlockMap>(),
                            tempThreeDimensionalMap,
                            tempBondMap);

    // The temporary data manager will have the same field specs and data as the real data manager.
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs = dataManager.getFieldSpecs();
    tempDataManager.allocateData(fieldSpecs);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data in the constitutiveData array.
    double *volume, *yReference, *y, *force;
    tempDataManager.getData(volumeFId, Field_ENUM::STEP_NONE)->ExtractView(&volume);
    tempDataManager.getData(tangentReferenceCoordinatesFId, Field_ENUM::STEP_NONE)->ExtractView(&yReference);
    tempDataManager.getData(coordinatesFId, Field_ENUM::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(forceDensityFId, Field_ENUM::STEP_NP1)->ExtractView(&force);

    // Create a temporary vector for storing force
    Teuchos::RCP<Epetra_Vector> forceVector = tempDataManager.getData(forceDensityFId, Field_ENUM::STEP_NP1);
    Teuchos::RCP<Epetra_Vector> tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
    double* tempForce;
    tempForceVector->ExtractView(&tempForce);

    Teuchos::RCP<Epetra_Vector> coordStepNP1;
    if(finiteDifferenceScheme == CONSISTENT_FORWARD_DIFFERENCE){
      // Store the current coordinates in a scratch vector
      coordStepNP1 = Teuchos::rcp(new Epetra_Vector(*tempDataManager.getData(coordinatesFId, Field_ENUM::STEP_NP1)));
      // Set the coordinates to the tangent reference coordinates
      for(int i=0 ; i<coordStepNP1->MyLength() ; ++i)
        y[i] = yReference[i];
    }
    
    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < 3*(numNeighbors+1))
      scratchMatrix.Resize(3*(numNeighbors+1));

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(3*(numNeighbors+1));
    for(int i=0 ; i<numNeighbors+1 ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<3 ; ++j)
        globalIndices[3*i+j] = 3*globalID+j;
    }

    if(finiteDifferenceScheme == FORWARD_DIFFERENCE || finiteDifferenceScheme == CONSISTENT_FORWARD_DIFFERENCE){
      // Compute and store the unperturbed force.
      computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
      for(int i=0 ; i<forceVector->MyLength() ; ++i)
        tempForce[i] = force[i];
    }

    // Perturb one dof in the neighborhood at a time and compute the force.
    // The point itself plus each of its neighbors must be perturbed.
    for(int iNID=0 ; iNID<numNeighbors+1 ; ++iNID){

      int perturbID;
      if(iNID < numNeighbors)
        perturbID = tempNeighborhoodList[iNID+1];
      else
        perturbID = 0;

      for(int dof=0 ; dof<3 ; ++dof){

        // Perturb a dof and compute the forces.
        double oldY = y[3*perturbID+dof];

        if(finiteDifferenceScheme == CENTRAL_DIFFERENCE){
          // Compute and store the negatively perturbed force.
          y[3*perturbID+dof] -= epsilon;
          computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
          y[3*perturbID+dof] = oldY;
          for(int i=0 ; i<forceVector->MyLength() ; ++i)
            tempForce[i] = force[i];
        }

        // Compute the positively perturbed force.
        double perturbation = epsilon;
        if(finiteDifferenceScheme == CONSISTENT_FORWARD_DIFFERENCE){
          perturbation = (*coordStepNP1)[3*perturbID+dof] - y[3*perturbID+dof];
          if(perturbation < 0.0)
            perturbation -= epsilon;
          else
            perturbation += epsilon;
        }

        y[3*perturbID+dof] += perturbation;
        computeForce(dt, tempNumOwnedPoints, &tempOwnedIDs[0], &tempNeighborhoodList[0], tempDataManager);
        y[3*perturbID+dof] = oldY;

        for(int i=0 ; i<numNeighbors+1 ; ++i){

          int forceID;
          if(i < numNeighbors)
            forceID = tempNeighborhoodList[i+1];
          else
            forceID = 0;

          for(int d=0 ; d<3 ; ++d){
            double value = ( force[3*forceID+d] - tempForce[3*forceID+d] ) / perturbation;
            if(finiteDifferenceScheme == CENTRAL_DIFFERENCE)
              value *= 0.5;
            scratchMatrix(3*forceID+d, 3*perturbID+dof) = value;
          }
        }
      }
    }

    // Convert force density to force
    // \todo Create utility function for this in ScratchMatrix
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        scratchMatrix(row, col) *= volume[row/3];
      }
    }

    // Check for NaNs
    for(unsigned int row=0 ; row<globalIndices.size() ; ++row){
      for(unsigned int col=0 ; col<globalIndices.size() ; ++col){
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(scratchMatrix(row, col)), "**** NaN detected in finite-difference Jacobian.\n");
      }
    }

    // Sum the values into the global tangent matrix (this is expensive).
    jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());
  }
}

void PeridigmNS::Material::computeApproximateDeformationGradient(const double dt,
                                                                 const int numOwnedPoints,
                                                                 const int* ownedIDs,
                                                                 const int* neighborhoodList,
                                                                 PeridigmNS::DataManager& dataManager) const
{

  // Get field ids for all relevant data
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId = fieldManager.getFieldId("Volume");
  int modelCoordinatesFId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFId = fieldManager.getFieldId("Coordinates");
  int weightedVolumeFId = fieldManager.getFieldId("WeightedVolume");
  int bondDamageFId = fieldManager.getFieldId("Bond_Damage");

  // Extract pointers to the underlying data
  double *volume, *x, *y, *weightedVolume, *bondDamage;
  dataManager.getData(volumeFId, Field_ENUM::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(modelCoordinatesFId, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(coordinatesFId, Field_ENUM::STEP_NP1)->ExtractView(&y);

  dataManager.getData(weightedVolumeFId, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(bondDamageFId, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);

  double X[3], Y[3];
  double shapeTensor[3][3];
  double shapeTensorInverse[3][3];
  double deformationGradientFirstTerm[3][3];
  double deformationGradient[3][3];

  int neighborhoodListIndex = 0;
  for(int ID=0 ; ID<numOwnedPoints ; ++ID){

    // Zero out data
    for(int i=0 ; i<3 ; ++i){
      for(int j=0 ; j<3 ; ++j){
        shapeTensor[i][j] = 0.0;
        deformationGradientFirstTerm[i][j] = 0.0;
      }
    }
    
    // \todo Include bond damage in this calculation.

    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      for(int i=0 ; i<3 ; ++i){
        X[i] = x[3*neighborID+i] - x[3*ID+i];
        Y[i] = y[3*neighborID+i] - y[3*ID+i];
      }
      double temp = volume[neighborID];
      for(int i=0 ; i<3 ; ++i){
        for(int j=0 ; j<3 ; ++j){
          shapeTensor[i][j] += temp*X[i]*X[j];
          deformationGradientFirstTerm[i][j] += temp*Y[i]*X[j];
        }
      }
      
      // \todo Check condition number of shape tensor.

      // Invert the shape tensor
      double minor[9] ;
      minor[0] = shapeTensor[1][1]*shapeTensor[2][2] - shapeTensor[1][2]*shapeTensor[2][1] ;
      minor[1] = shapeTensor[1][0]*shapeTensor[2][2] - shapeTensor[1][2]*shapeTensor[2][0] ;
      minor[2] = shapeTensor[1][0]*shapeTensor[2][1] - shapeTensor[1][1]*shapeTensor[2][0] ;
      minor[3] = shapeTensor[0][1]*shapeTensor[2][2] - shapeTensor[0][2]*shapeTensor[2][1] ;
      minor[4] = shapeTensor[0][0]*shapeTensor[2][2] - shapeTensor[2][0]*shapeTensor[0][2] ;
      minor[5] = shapeTensor[0][0]*shapeTensor[2][1] - shapeTensor[0][1]*shapeTensor[2][0] ;
      minor[6] = shapeTensor[0][1]*shapeTensor[1][2] - shapeTensor[0][2]*shapeTensor[1][1] ;
      minor[7] = shapeTensor[0][0]*shapeTensor[1][2] - shapeTensor[0][2]*shapeTensor[1][0] ;
      minor[8] = shapeTensor[0][0]*shapeTensor[1][1] - shapeTensor[0][1]*shapeTensor[1][0] ;
      double det = shapeTensor[0][0]*minor[0] - shapeTensor[0][1]*minor[1] + shapeTensor[0][2]*minor[2] ;
      shapeTensorInverse[0][0] = minor[0]/det;
      shapeTensorInverse[0][1] = -1.0*minor[3]/det;
      shapeTensorInverse[0][2] = minor[6]/det;
      shapeTensorInverse[1][0] = -1.0*minor[1]/det;
      shapeTensorInverse[1][1] = minor[4]/det;
      shapeTensorInverse[1][2] = -1.0*minor[7]/det;
      shapeTensorInverse[2][0] = minor[2]/det;
      shapeTensorInverse[2][1] = -1.0*minor[5]/det;
      shapeTensorInverse[2][2] = minor[8]/det;

      // Compute deformation gradient
      for(int i=0 ; i<3 ; ++i){
        for(int j=0 ; j<3 ; ++j){
          deformationGradient[i][j] = 
            deformationGradientFirstTerm[i][0]*shapeTensorInverse[0][i] +
            deformationGradientFirstTerm[i][1]*shapeTensorInverse[1][i] +
            deformationGradientFirstTerm[i][2]*shapeTensorInverse[2][i];
        }
      }

    }

  }
}
