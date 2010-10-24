/*! \file Peridigm_LinearElasticIsotropicMaterial.cpp */

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

#include "Peridigm_LinearElasticIsotropicMaterial.hpp"
#include "Peridigm_CriticalStretchDamageModel.hpp"
#include <Teuchos_TestForException.hpp>
#include "PdMaterialUtilities.h"

using namespace std;

PeridigmNS::LinearElasticIsotropicMaterial::LinearElasticIsotropicMaterial(const Teuchos::ParameterList& params)
  : Material(params),
    m_decompStates(),
    m_damageModel()
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = params.get<double>("Bulk Modulus");
  m_shearModulus = params.get<double>("Shear Modulus");
  m_density = params.get<double>("Density");

  if(params.isSublist("Damage Model")){
    Teuchos::ParameterList damageParams = params.sublist("Damage Model");
    if(!damageParams.isParameter("Type")){
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         "Damage model \"Type\" not specified in Damage Model parameter list.");
    }
    string& damageModelType = damageParams.get<string>("Type");
    if(damageModelType == "Critical Stretch"){
      m_damageModel = Teuchos::rcp(new PeridigmNS::CriticalStretchDamageModel(damageParams));
    }
    else{
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         "Invalid damage model, \"None\" or \"Critical Stretch\" required.");
    }

    // query damage model for number of constitutive data and the names of these data

    // -> critical stretch says it needs scalar y

    // check to see if new storage is required

    // -> mat model sees that it already has y

    // add new storage if required

    // -> result is no new storage
  }

  // set up vector of variable specs
  m_variableSpecs = Teuchos::rcp(new vector<Field_NS::FieldSpec>);
  m_variableSpecs->push_back(Field_NS::VOLUME);
  m_variableSpecs->push_back(Field_NS::DAMAGE);
  m_variableSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
  m_variableSpecs->push_back(Field_NS::DILATATION);
  m_variableSpecs->push_back(Field_NS::COORD3D);
  m_variableSpecs->push_back(Field_NS::DISPL3D);
  m_variableSpecs->push_back(Field_NS::VELOC3D);
  m_variableSpecs->push_back(Field_NS::ACCEL3D);
  m_variableSpecs->push_back(Field_NS::FORCE3D);
}

PeridigmNS::LinearElasticIsotropicMaterial::~LinearElasticIsotropicMaterial()
{
}

void
PeridigmNS::LinearElasticIsotropicMaterial::initialize(const Epetra_Vector& x,
                                                     const Epetra_Vector& u,
                                                     const Epetra_Vector& v,
                                                     const double dt,
                                                     const Epetra_Vector& cellVolume,
                                                     const int numOwnedPoints,
                                                     const int* ownedIDs,
                                                     const int* neighborhoodList,
                                                     double* bondState,
                                                     Epetra_MultiVector& scalarConstitutiveData,
                                                     Epetra_MultiVector& vectorConstitutiveData,
                                                     Epetra_MultiVector& bondConstitutiveData,
                                                     Epetra_Vector& force) const
{
  // Sanity checks on vector sizes
  TEST_FOR_EXCEPT_MSG(x.MyLength() != u.MyLength(), 
					  "x and u vector lengths do not match\n");
  TEST_FOR_EXCEPT_MSG(x.MyLength() != v.MyLength(), 
					  "x and v vector lengths do not match\n");
  TEST_FOR_EXCEPT_MSG(x.MyLength() != vectorConstitutiveData.MyLength(), 
					  "x and vector constitutive data vector lengths do not match\n");
  TEST_FOR_EXCEPT_MSG(x.MyLength() != force.MyLength(), 
					  "x and force vector lengths do not match\n");
  TEST_FOR_EXCEPT_MSG(cellVolume.MyLength() != scalarConstitutiveData.MyLength(), 
					  "cellVolume and scalar constitutive data vector lengths do not match\n");

  //! \todo Create structure for storing influence function values.
//   double omega = 1.0;

  // Initialize data fields
  scalarConstitutiveData.PutScalar(0.0);
  vectorConstitutiveData.PutScalar(0.0);
  bondConstitutiveData.PutScalar(0.0);
  force.PutScalar(0.0);
  int neighborhoodListIndex = 0;
  int bondStateIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
	  bondState[bondStateIndex++] = 0.0;
      neighborhoodListIndex++;
    }
  }

  TEST_FOR_EXCEPT_MSG(bondConstitutiveData.MyLength() != bondStateIndex,
					  "bondConstitutiveData vector and bondState array lengths do not match\n");

  // Extract pointers to the underlying data in the constitutiveData array
  std::pair<int,double*> scalarView = m_decompStates.extractStrideView(scalarConstitutiveData);
  double* weightedVolume = m_decompStates.extractWeightedVolumeView(scalarView);

  PdMaterialUtilities::computeWeightedVolume(x.Values(),cellVolume.Values(),weightedVolume,numOwnedPoints,neighborhoodList);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::updateConstitutiveData(const Epetra_Vector& x,
																 const Epetra_Vector& u,
																 const Epetra_Vector& v,
																 const double dt,
																 const Epetra_Vector& cellVolume,
																 const int numOwnedPoints,
																 const int* ownedIDs,
																 const int* neighborhoodList,
																 double* bondState,
																 Epetra_MultiVector& scalarConstitutiveData,
																 Epetra_MultiVector& vectorConstitutiveData,
																 Epetra_MultiVector& bondConstitutiveData,
																 Epetra_Vector& force) const
{
  //! \todo Create structure for storing influence function values.
//   double omega = 1.0;
  
  // Extract pointers to the underlying data in the constitutiveData array
  std::pair<int,double*> scalarView = m_decompStates.extractStrideView(scalarConstitutiveData);
  double* weightedVolume = m_decompStates.extractWeightedVolumeView(scalarView);
  double* dilatation = m_decompStates.extractDilatationView(scalarView);
  double* damage = m_decompStates.extractDamageView(scalarView);

  std::pair<int,double*> vectorView = m_decompStates.extractStrideView(vectorConstitutiveData);
  double *y = m_decompStates.extractCurrentPositionView(vectorView);

  // Update the geometry
  PdMaterialUtilities::updateGeometry(x.Values(),u.Values(),v.Values(),y,x.MyLength(),dt);

  // Update the bondState
  if(!m_damageModel.is_null()){
    m_damageModel->computeDamage(x,
								 u,
								 v,
                                 dt,
                                 cellVolume,
                                 numOwnedPoints,
                                 ownedIDs,
                                 neighborhoodList,
                                 bondState,
                                 scalarConstitutiveData,
                                 vectorConstitutiveData,
                                 bondConstitutiveData,
                                 force);
  }

  //  Update the element damage (percent of bonds broken)
  int neighborhoodListIndex = 0;
  int bondStateIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
	int nodeID = ownedIDs[iID];
	int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
	double totalDamage = 0.0;
	for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
	  totalDamage += bondState[bondStateIndex++];
	}
	if(numNeighbors > 0)
	  totalDamage /= numNeighbors;
	else
	  totalDamage = 0.0;
 	damage[nodeID] = totalDamage;
  }

PdMaterialUtilities::computeDilatation(x.Values(),y,weightedVolume,cellVolume.Values(),bondState,dilatation,neighborhoodList,numOwnedPoints);
}

void
PeridigmNS::LinearElasticIsotropicMaterial::computeForce(const Epetra_Vector& x,
													   const Epetra_Vector& u,
													   const Epetra_Vector& v,
													   const double dt,
													   const Epetra_Vector& cellVolume,
													   const int numOwnedPoints,
													   const int* ownedIDs,
													   const int* neighborhoodList,
													   double* bondState,
													   Epetra_MultiVector& scalarConstitutiveData,
													   Epetra_MultiVector& vectorConstitutiveData,
													   Epetra_MultiVector& bondConstitutiveData,
													   Epetra_Vector& force) const
{
//   double omega = 1.0;

  // Extract pointers to the underlying data in the constitutiveData array
  std::pair<int,double*> scalarView = m_decompStates.extractStrideView(scalarConstitutiveData);
  double* weightedVolume = m_decompStates.extractWeightedVolumeView(scalarView);
  double* dilatation = m_decompStates.extractDilatationView(scalarView);
  std::pair<int,double*> vectorView = m_decompStates.extractStrideView(vectorConstitutiveData);
  double *y = m_decompStates.extractCurrentPositionView(vectorView);

  // Compute the force on each particle that results from interactions
  // with locally-owned nodes
  force.PutScalar(0.0);

  PdMaterialUtilities::computeInternalForceLinearElastic(x.Values(),y,weightedVolume,cellVolume.Values(),dilatation,bondState,force.Values(),neighborhoodList,numOwnedPoints,m_bulkModulus,m_shearModulus);

}
