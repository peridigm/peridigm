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

    // The temporary data manager will have the same fields and data as the real data manager.
    vector<int> fieldIds = dataManager.getFieldIds();
    tempDataManager.allocateData(fieldIds);
    tempDataManager.copyLocallyOwnedDataFromDataManager(dataManager);

    // Set up numOwnedPoints and ownedIDs.
    // There is only one owned ID, and it has local ID zero in the tempDataManager.
    int tempNumOwnedPoints = 1;
    vector<int> tempOwnedIDs(1);
    tempOwnedIDs[0] = 0;

    // Extract pointers to the underlying data in the constitutiveData array.
    double *volume, *y, *force;
    tempDataManager.getData(volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    tempDataManager.getData(coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1)->ExtractView(&force);

    // Create a temporary vector for storing force
    Teuchos::RCP<Epetra_Vector> forceVector = tempDataManager.getData(forceDensityFId, PeridigmField::STEP_NP1);
    Teuchos::RCP<Epetra_Vector> tempForceVector = Teuchos::rcp(new Epetra_Vector(*forceVector));
    double* tempForce;
    tempForceVector->ExtractView(&tempForce);

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

    if(finiteDifferenceScheme == FORWARD_DIFFERENCE){
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


	// Compute the purturbed force
        y[3*perturbID+dof] += epsilon;
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

double PeridigmNS::Material::calculateBulkModulus(const Teuchos::ParameterList & params) const
{
  bool bulkModulusDefined(false), shearModulusDefined(false), youngsModulusDefined(false), poissonsRatioDefined(false);
  double bulkModulus(0.0), shearModulus(0.0), youngsModulus(0.0), poissonsRatio(0.0);
  double computedValue;

  if( params.isParameter("Bulk Modulus") ){
    bulkModulusDefined = true;
    bulkModulus = params.get<double>("Bulk Modulus");
  }
  if( params.isParameter("Shear Modulus") ){
    shearModulus = params.get<double>("Shear Modulus");
    shearModulusDefined = true;
  }
  if( params.isParameter("Young's Modulus") ){
    youngsModulus = params.get<double>("Young's Modulus");
    youngsModulusDefined = true;
  }
  if( params.isParameter("Poisson's Ratio") ){
    poissonsRatio = params.get<double>("Poisson's Ratio");
    poissonsRatioDefined = true;
  }

  int numDefinedConstants = static_cast<int>(bulkModulusDefined) + 
    static_cast<int>(shearModulusDefined) + 
    static_cast<int>(youngsModulusDefined) + 
    static_cast<int>(poissonsRatioDefined);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(numDefinedConstants != 2, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  if(bulkModulusDefined)
    computedValue = bulkModulus;
  else if(youngsModulusDefined && shearModulusDefined)
    computedValue = (youngsModulus * shearModulus) / (3.0*(3.0*shearModulus - youngsModulus));
  else if(youngsModulusDefined && poissonsRatioDefined)
    computedValue = youngsModulus / (3.0*(1.0 - 2.0*poissonsRatio));
  else if(shearModulusDefined && poissonsRatioDefined)
    computedValue = (2.0*shearModulus*(1.0 + poissonsRatio)) / (3.0*(1.0 - 2.0*poissonsRatio));
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  return computedValue;
}


double PeridigmNS::Material::calculateShearModulus(const Teuchos::ParameterList & params) const
{
  bool bulkModulusDefined(false), shearModulusDefined(false), youngsModulusDefined(false), poissonsRatioDefined(false);
  double bulkModulus(0.0), shearModulus(0.0), youngsModulus(0.0), poissonsRatio(0.0);
  double computedValue;

  if( params.isParameter("Bulk Modulus") ){
    bulkModulusDefined = true;
    bulkModulus = params.get<double>("Bulk Modulus");
  }
  if( params.isParameter("Shear Modulus") ){
    shearModulus = params.get<double>("Shear Modulus");
    shearModulusDefined = true;
  }
  if( params.isParameter("Young's Modulus") ){
    youngsModulus = params.get<double>("Young's Modulus");
    youngsModulusDefined = true;
  }
  if( params.isParameter("Poisson's Ratio") ){
    poissonsRatio = params.get<double>("Poisson's Ratio");
    poissonsRatioDefined = true;
  }

  int numDefinedConstants = static_cast<int>(bulkModulusDefined) + 
    static_cast<int>(shearModulusDefined) + 
    static_cast<int>(youngsModulusDefined) + 
    static_cast<int>(poissonsRatioDefined);

  TEUCHOS_TEST_FOR_EXCEPT_MSG(numDefinedConstants != 2, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  if(shearModulusDefined)
    computedValue = shearModulus;
  else if(bulkModulusDefined && youngsModulusDefined)
    computedValue = (3.0*bulkModulus*youngsModulus) / (9.0*bulkModulus - youngsModulus);
  else if(bulkModulusDefined & poissonsRatioDefined)
    computedValue = (3.0*bulkModulus*(1.0 - 2.0*poissonsRatio)) / (2.0*(1.0 + poissonsRatio));
  else if(youngsModulusDefined && poissonsRatioDefined)
    computedValue = youngsModulus / (2.0*(1.0 + poissonsRatio));
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  Exactly two elastic constants must be provided.  Allowable constants are \"Bulk Modulus\", \"Shear Modulus\", \"Young's Modulus\", \"Poisson's Ratio\".\n");

  return computedValue;
}

void PeridigmNS::Material::computeApproximateDeformationGradient(const int numOwnedPoints,
                                                                 const int* ownedIDs,
                                                                 const int* neighborhoodList,
                                                                 PeridigmNS::DataManager& dataManager) const
{

  // Get field ids for all relevant data
  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  int volumeFId = fieldManager.getFieldId("Volume");
  int modelCoordinatesFId = fieldManager.getFieldId("Model_Coordinates");
  int coordinatesFId = fieldManager.getFieldId("Coordinates");
  int weightedVolumeFId = fieldManager.getFieldId("Weighted_Volume");
  int bondDamageFId = fieldManager.getFieldId("Bond_Damage");
  int deformationGradientXXFId = fieldManager.getFieldId("Deformation_GradientXX");
  int deformationGradientXYFId = fieldManager.getFieldId("Deformation_GradientXY");
  int deformationGradientXZFId = fieldManager.getFieldId("Deformation_GradientXZ");
  int deformationGradientYXFId = fieldManager.getFieldId("Deformation_GradientYX");
  int deformationGradientYYFId = fieldManager.getFieldId("Deformation_GradientYY");
  int deformationGradientYZFId = fieldManager.getFieldId("Deformation_GradientYZ");
  int deformationGradientZXFId = fieldManager.getFieldId("Deformation_GradientZX");
  int deformationGradientZYFId = fieldManager.getFieldId("Deformation_GradientZY");
  int deformationGradientZZFId = fieldManager.getFieldId("Deformation_GradientZZ");

  // Extract pointers to the underlying data
  double *volume, *x, *y, *weightedVolume, *bondDamage;
  double *deformationGradientXX, *deformationGradientXY, *deformationGradientXZ;
  double *deformationGradientYX, *deformationGradientYY, *deformationGradientYZ;
  double *deformationGradientZX, *deformationGradientZY, *deformationGradientZZ;
  dataManager.getData(volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(modelCoordinatesFId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(weightedVolumeFId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(bondDamageFId, PeridigmField::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(deformationGradientXXFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientXX);
  dataManager.getData(deformationGradientXYFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientXY);
  dataManager.getData(deformationGradientXZFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientXZ);
  dataManager.getData(deformationGradientYXFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientYX);
  dataManager.getData(deformationGradientYYFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientYY);
  dataManager.getData(deformationGradientYZFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientYZ);
  dataManager.getData(deformationGradientZXFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZX);
  dataManager.getData(deformationGradientZYFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZY);
  dataManager.getData(deformationGradientZZFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradientZZ);

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
    // \todo Include influence function in this calculation.

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
          deformationGradientFirstTerm[i][0]*shapeTensorInverse[0][j] +
          deformationGradientFirstTerm[i][1]*shapeTensorInverse[1][j] +
          deformationGradientFirstTerm[i][2]*shapeTensorInverse[2][j];
      }
    }

    // Copy result into data manager
    deformationGradientXX[ID] = deformationGradient[0][0];
    deformationGradientXY[ID] = deformationGradient[0][1];
    deformationGradientXZ[ID] = deformationGradient[0][2];
    deformationGradientYX[ID] = deformationGradient[1][0];
    deformationGradientYY[ID] = deformationGradient[1][1];
    deformationGradientYZ[ID] = deformationGradient[1][2];
    deformationGradientZX[ID] = deformationGradient[2][0];
    deformationGradientZY[ID] = deformationGradient[2][1];
    deformationGradientZZ[ID] = deformationGradient[2][2];
  }
}
