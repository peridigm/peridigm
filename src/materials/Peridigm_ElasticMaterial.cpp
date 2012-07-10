/*! \file Peridigm_LinearElasticIsotropicMaterial.cpp */

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

#include "Peridigm_ElasticMaterial.hpp"
#include "Peridigm_Timer.hpp"
#include "elastic.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>

using namespace std;

PeridigmNS::LinearElasticIsotropicMaterial::LinearElasticIsotropicMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_applyAutomaticDifferentiationJacobian(true)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_shearModulus = params.get<double>("Shear Modulus");
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");
  if(params.isParameter("Apply Automatic Differentiation Jacobian"))
    m_applyAutomaticDifferentiationJacobian = params.get<bool>("Apply Automatic Differentiation Jacobian");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(params.isParameter("Apply Shear Correction Factor"), "**** Error:  Shear Correction Factor not supported for Viscoelastic material.\n");

  // set up vector of variable specs
  m_variableSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
  m_variableSpecs->push_back(Field_NS::VOLUME);
  m_variableSpecs->push_back(Field_NS::DAMAGE);
  m_variableSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
  m_variableSpecs->push_back(Field_NS::DILATATION);
  m_variableSpecs->push_back(Field_NS::COORD3D);
  m_variableSpecs->push_back(Field_NS::CURCOORD3D);
  m_variableSpecs->push_back(Field_NS::FORCE_DENSITY3D);
  m_variableSpecs->push_back(Field_NS::BOND_DAMAGE);
}

PeridigmNS::LinearElasticIsotropicMaterial::~LinearElasticIsotropicMaterial()
{
}

void
PeridigmNS::LinearElasticIsotropicMaterial::initialize(const double dt,
                                                       const int numOwnedPoints,
                                                       const int* ownedIDs,
                                                       const int* neighborhoodList,
                                                       PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *cellVolume, *weightedVolume;
  dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);

  MATERIAL_EVALUATION::computeWeightedVolume(x,cellVolume,weightedVolume,numOwnedPoints,neighborhoodList);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::updateConstitutiveData(const double dt,
                                                                   const int numOwnedPoints,
                                                                   const int* ownedIDs,
                                                                   const int* neighborhoodList,
                                                                   PeridigmNS::DataManager& dataManager) const
{
  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *y, *cellVolume, *weightedVolume, *dilatation, *damage, *bondDamage;
  dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);

  MATERIAL_EVALUATION::computeDilatationAD(x,y,weightedVolume,cellVolume,bondDamage,dilatation,neighborhoodList,numOwnedPoints);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::computeForce(const double dt,
                                                         const int numOwnedPoints,
                                                         const int* ownedIDs,
                                                         const int* neighborhoodList,
                                                         PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *weightedVolume, *dilatation, *damage, *bondDamage, *force;
  dataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(Field_NS::DILATATION, Field_ENUM::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);
  dataManager.getData(Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1)->ExtractView(&force);

  MATERIAL_EVALUATION::computeInternalForceLinearElasticAD(x,y,weightedVolume,cellVolume,dilatation,bondDamage,force,neighborhoodList,numOwnedPoints,m_bulkModulus,m_shearModulus);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::computeJacobian(const double dt,
                                                            const int numOwnedPoints,
                                                            const int* ownedIDs,
                                                            const int* neighborhoodList,
                                                            PeridigmNS::DataManager& dataManager,
                                                            PeridigmNS::SerialMatrix& jacobian) const
{
  if(m_applyAutomaticDifferentiationJacobian){
    // Compute the Jacobian via automatic differentiation
    computeAutomaticDifferentiationJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian);  
  }
  else{
    // Call the base class function, which computes the Jacobian by finite difference
    PeridigmNS::Material::computeJacobian(dt, numOwnedPoints, ownedIDs, neighborhoodList, dataManager, jacobian);
  }
}


void
PeridigmNS::LinearElasticIsotropicMaterial::computeAutomaticDifferentiationJacobian(const double dt,
                                                                                    const int numOwnedPoints,
                                                                                    const int* ownedIDs,
                                                                                    const int* neighborhoodList,
                                                                                    PeridigmNS::DataManager& dataManager,
                                                                                    PeridigmNS::SerialMatrix& jacobian) const
{
  // Compute contributions to the tangent matrix on an element-by-element basis

  // To reduce memory re-allocation, use static variable to store Fad types for
  // current coordinates (independent variables).
  static vector<Sacado::Fad::DFad<double> > y_AD;

  // Loop over all points.
  int neighborhoodListIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){

    PeridigmNS::Timer::self().startTimer("AD Jacobian General Set Up");

    // Create a temporary neighborhood consisting of a single point and its neighbors.
    // The temporary neighborhood is sorted by global ID to somewhat increase the chance
    // that the downstream Epetra_CrsMatrix::SumIntoMyValues() calls will be as efficient
    // as possible.
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    int numEntries = numNeighbors+1;
    int numDof = 3*numEntries;
    vector<int> tempMyGlobalIDs(numEntries);
    // Create a placeholder that will appear at the beginning of the sorted list.
    tempMyGlobalIDs[0] = -1;
    vector<int> tempNeighborhoodList(numEntries); 
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
    Teuchos::RCP<Epetra_BlockMap> tempOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 1, 0, serialComm));
    Teuchos::RCP<Epetra_BlockMap> tempThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(numEntries, numEntries, &tempMyGlobalIDs[0], 3, 0, serialComm));
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
    vector<int> tempOwnedIDs(tempNumOwnedPoints);
    tempOwnedIDs[0] = 0;

    // Use the scratchMatrix as sub-matrix for storing tangent values prior to loading them into the global tangent matrix.
    // Resize scratchMatrix if necessary
    if(scratchMatrix.Dimension() < numDof)
      scratchMatrix.Resize(numDof);

    // Create a list of global indices for the rows/columns in the scratch matrix.
    vector<int> globalIndices(numDof);
    for(int i=0 ; i<numEntries ; ++i){
      int globalID = tempOneDimensionalMap->GID(i);
      for(int j=0 ; j<3 ; ++j)
        globalIndices[3*i+j] = 3*globalID+j;
    }

    // Extract pointers to the underlying data in the constitutiveData array.
    double *x, *y, *cellVolume, *weightedVolume, *damage, *bondDamage;
    tempDataManager.getData(Field_NS::COORD3D, Field_ENUM::STEP_NONE)->ExtractView(&x);
    tempDataManager.getData(Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1)->ExtractView(&y);
    tempDataManager.getData(Field_NS::VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&cellVolume);
    tempDataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_ENUM::STEP_NONE)->ExtractView(&weightedVolume);
    tempDataManager.getData(Field_NS::DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&damage);
    tempDataManager.getData(Field_NS::BOND_DAMAGE, Field_ENUM::STEP_NP1)->ExtractView(&bondDamage);

    PeridigmNS::Timer::self().stopTimer("AD Jacobian General Set Up");
    PeridigmNS::Timer::self().startTimer("AD Jacobian FAD Set Up");

    // Create arrays of Fad objects for the current coordinates, dilatation, and force density
    // Modify the existing vector of Fad objects for the current coordinates
    if((int)y_AD.size() < numDof)
      y_AD.resize(numDof);
    for(int i=0 ; i<numDof ; ++i){
      y_AD[i].diff(i, numDof);
      y_AD[i].val() = y[i];
    }
    // Create vectors of empty AD types for the dependent variables
    vector<Sacado::Fad::DFad<double> > dilatation_AD(numEntries);
    vector<Sacado::Fad::DFad<double> > force_AD(numDof);

    PeridigmNS::Timer::self().stopTimer("AD Jacobian FAD Set Up");
    PeridigmNS::Timer::self().startTimer("AD Jacobian Constitutive Model");

    // Evaluate the constitutive model using the AD types
    MATERIAL_EVALUATION::computeDilatationAD(x,&y_AD[0],weightedVolume,cellVolume,bondDamage,&dilatation_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints);
    MATERIAL_EVALUATION::computeInternalForceLinearElasticAD(x,&y_AD[0],weightedVolume,cellVolume,&dilatation_AD[0],bondDamage,&force_AD[0],&tempNeighborhoodList[0],tempNumOwnedPoints,m_bulkModulus,m_shearModulus);

    PeridigmNS::Timer::self().stopTimer("AD Jacobian Constitutive Model");
    PeridigmNS::Timer::self().startTimer("AD Jacobian Global Fill");

    // Load derivative values into scratch matrix
    // Multiply by volume along the way to convert force density to force
    for(int row=0 ; row<numDof ; ++row){
      for(int col=0 ; col<numDof ; ++col){
        scratchMatrix(row, col) = force_AD[row].dx(col) * cellVolume[row/3];
      }
    }

    // Sum the values into the global tangent matrix (this is expensive).
    jacobian.addValues((int)globalIndices.size(), &globalIndices[0], scratchMatrix.Data());

    PeridigmNS::Timer::self().stopTimer("AD Jacobian Global Fill");
  }
}
