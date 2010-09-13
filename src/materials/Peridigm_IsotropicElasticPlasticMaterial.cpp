/*
 * Peridigm_IsotropicElasticPlasticMaterial.cxx
 *
 */
#include "Peridigm_IsotropicElasticPlasticMaterial.hpp"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "PdMaterialUtilities.h"


PeridigmNS::IsotropicElasticPlasticMaterial::IsotropicElasticPlasticMaterial(const Teuchos::ParameterList & params)
:
Material(params),
decompStates()
{

	decompStates.addScalarStateBondVariable("scalarPlasticExtensionState");

	//! \todo Add meaningful asserts on material properties.
	m_bulkModulus = params.get<double>("Bulk Modulus");
	m_shearModulus = params.get<double>("Shear Modulus");
	m_horizon = params.get<double>("Material Horizon");
	m_density = params.get<double>("Density");
	m_yieldStress = params.get<double>("Yield Stress");
}


PeridigmNS::IsotropicElasticPlasticMaterial::~IsotropicElasticPlasticMaterial()
{
}


void PeridigmNS::IsotropicElasticPlasticMaterial::initialize(const Epetra_Vector& x,
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
//	  double omega = 1.0;

	  // Initialize data fields
	  scalarConstitutiveData.PutScalar(0.0);
	  vectorConstitutiveData.PutScalar(0.0);
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

	  // Extract pointers to the underlying data in the constitutiveData array
	  std::pair<int,double*> scalarView = decompStates.extractStrideView(scalarConstitutiveData);
	  double* weightedVolume = decompStates.extractWeightedVolumeView(scalarView);

	  PdMaterialUtilities::computeWeightedVolume(x.Values(),cellVolume.Values(),weightedVolume,numOwnedPoints,neighborhoodList);

}

void
PeridigmNS::IsotropicElasticPlasticMaterial::updateConstitutiveData(const Epetra_Vector& x,
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

	// Extract pointers to the underlying data in the constitutiveData array
	std::pair<int,double*> scalarView = decompStates.extractStrideView(scalarConstitutiveData);
	double* weightedVolume = decompStates.extractWeightedVolumeView(scalarView);
	double* dilatation = decompStates.extractDilatationView(scalarView);

	std::pair<int,double*> vectorView = decompStates.extractStrideView(vectorConstitutiveData);
	double *y = decompStates.extractCurrentPositionView(vectorView);

	// Update the geometry
	PdMaterialUtilities::updateGeometry(x.Values(),u.Values(),v.Values(),y,x.MyLength(),dt);


	PdMaterialUtilities::computeDilatation(x.Values(),y,weightedVolume,cellVolume.Values(),bondState,dilatation,neighborhoodList,numOwnedPoints);
}

void
PeridigmNS::IsotropicElasticPlasticMaterial::computeForce(const Epetra_Vector& x,
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

	  // Extract pointers to the underlying data in the constitutiveData array
	  std::pair<int,double*> scalarView = decompStates.extractStrideView(scalarConstitutiveData);
	  double* weightedVolume = decompStates.extractWeightedVolumeView(scalarView);
	  double* dilatation = decompStates.extractDilatationView(scalarView);
	//	double* damage = decompStates.extractDamageView(scalarView);
	  std::pair<int,double*> vectorView = decompStates.extractStrideView(vectorConstitutiveData);
	  double *y = decompStates.extractCurrentPositionView(vectorView);
	  std::pair<int,double*> scalarBondView = decompStates.extractStrideView(bondConstitutiveData);
	  double* edpN = decompStates.extractScalarBondVariable(scalarBondView,"scalarPlasticExtensionState");

	  // Compute the force on each particle that results from interactions
	  // with locally-owned nodes
	  force.PutScalar(0.0);

//	  PdMaterialUtilities::computeInternalForceIsotropicElasticPlastic
//	  (
//			  x.Values(),
//			  y,
//			  weightedVolume,
//			  cellVolume.Values(),
//			  dilatation,
//			  bondState,
//			  edpN,
//			  force.Values(),
//			  neighborhoodList,
//			  numOwnedPoints,
//			  m_bulkModulus,
//			  m_shearModulus,
//			  m_horizon,
//			  m_yieldStress
//	  );

}

