/*
 * Peridigm_IsotropicElasticPlasticMaterial.cxx
 *
 */
#include "Peridigm_IsotropicElasticPlasticMaterial.hpp"
#include "Peridigm_CriticalStretchDamageModel.hpp"
#include <Teuchos_TestForException.hpp>
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "PdMaterialUtilities.h"
#include <limits>


PeridigmNS::IsotropicElasticPlasticMaterial::IsotropicElasticPlasticMaterial(const Teuchos::ParameterList & params)
:
Material(params),
m_damageModel()
{
	//! \todo Add meaningful asserts on material properties.
	m_bulkModulus = params.get<double>("Bulk Modulus");
	m_shearModulus = params.get<double>("Shear Modulus");
	m_horizon = params.get<double>("Material Horizon");
	m_density = params.get<double>("Density");
	m_yieldStress = params.get<double>("Yield Stress");

	/*
	 * Set the yield stress to a very large value: this in effect makes the model run elastic -- useful for testing
	 */
	if(params.isType<string>("Test"))
		m_yieldStress = std::numeric_limits<double>::max( );

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
	}

    // set up vector of variable specs
    m_variableSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);
    m_variableSpecs->push_back(Field_NS::VOLUME);
    m_variableSpecs->push_back(Field_NS::DAMAGE);
    m_variableSpecs->push_back(Field_NS::WEIGHTED_VOLUME);
    m_variableSpecs->push_back(Field_NS::DILATATION);
    m_variableSpecs->push_back(Field_NS::COORD3D);
    m_variableSpecs->push_back(Field_NS::CURCOORD3D);
    m_variableSpecs->push_back(Field_NS::FORCE3D);
    m_variableSpecs->push_back(Field_NS::DEVIATORIC_PLASTIC_EXTENSION);
    m_variableSpecs->push_back(Field_NS::LAMBDA);

}


PeridigmNS::IsotropicElasticPlasticMaterial::~IsotropicElasticPlasticMaterial()
{
}


void PeridigmNS::IsotropicElasticPlasticMaterial::initialize(const Epetra_Vector& u,
                                                             const Epetra_Vector& v,
                                                             const double dt,
                                                             const int numOwnedPoints,
                                                             const int* ownedIDs,
                                                             const int* neighborhoodList,
                                                             double* bondState,
                                                             PeridigmNS::DataManager& dataManager,
                                                             Epetra_Vector& force) const
{
	  // Initialize data fields
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

	  // Extract pointers to the underlying data
      double *x, *cellVolume, *weightedVolume;
      dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
      dataManager.getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&cellVolume);
      dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&weightedVolume);

	  PdMaterialUtilities::computeWeightedVolume(x,cellVolume,weightedVolume,numOwnedPoints,neighborhoodList);
}

void
PeridigmNS::IsotropicElasticPlasticMaterial::updateConstitutiveData(const Epetra_Vector& u,
                                                                    const Epetra_Vector& v,
                                                                    const double dt,
                                                                    const int numOwnedPoints,
                                                                    const int* ownedIDs,
                                                                    const int* neighborhoodList,
                                                                    double* bondState,
                                                                    PeridigmNS::DataManager& dataManager,
                                                                    Epetra_Vector& force) const
{
  int vectorLength = dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->MyLength();

  // Extract pointers to the underlying data in the constitutiveData array
  double *x, *y, *volume, *dilatation, *damage, *weightedVolume;
  dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
  dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);
  dataManager.getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&volume);
  dataManager.getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(Field_NS::DAMAGE, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&damage);
  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&weightedVolume);

	// Update the geometry
	PdMaterialUtilities::updateGeometry(x,u.Values(),v.Values(),y,vectorLength,dt);

	// Update the bondState
	if(!m_damageModel.is_null()){
		m_damageModel->computeDamage(u,
                                     v,
                                     dt,
                                     numOwnedPoints,
                                     ownedIDs,
                                     neighborhoodList,
                                     bondState,
                                     dataManager,
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


	PdMaterialUtilities::computeDilatation(x,y,weightedVolume,volume,bondState,dilatation,neighborhoodList,numOwnedPoints);
}

void
PeridigmNS::IsotropicElasticPlasticMaterial::computeForce(const Epetra_Vector& u,
                                                          const Epetra_Vector& v,
                                                          const double dt,
                                                          const int numOwnedPoints,
                                                          const int* ownedIDs,
                                                          const int* neighborhoodList,
                                                          double* bondState,
                                                          PeridigmNS::DataManager& dataManager,
                                                          Epetra_Vector& force) const
{

	  // Extract pointers to the underlying data in the constitutiveData array
      double *x, *y;
      dataManager.getData(Field_NS::COORD3D, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&x);
      dataManager.getData(Field_NS::CURCOORD3D, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&y);

      double* dilatation;
      dataManager.getData(Field_NS::DILATATION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&dilatation);

	  double* edpN;
	  double* edpNP1;
	  dataManager.getData(Field_NS::DEVIATORIC_PLASTIC_EXTENSION, Field_NS::FieldSpec::STEP_N)->ExtractView(&edpN);
	  dataManager.getData(Field_NS::DEVIATORIC_PLASTIC_EXTENSION, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&edpNP1);

	  double* weightedVolume;
	  dataManager.getData(Field_NS::WEIGHTED_VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&weightedVolume);

      double *volume;
      dataManager.getData(Field_NS::VOLUME, Field_NS::FieldSpec::STEP_NONE)->ExtractView(&volume);

	  double* lambdaN;
	  double* lambdaNP1;
	  dataManager.getData(Field_NS::LAMBDA, Field_NS::FieldSpec::STEP_N)->ExtractView(&lambdaN);
	  dataManager.getData(Field_NS::LAMBDA, Field_NS::FieldSpec::STEP_NP1)->ExtractView(&lambdaNP1);


	  // Compute the force on each particle that results from interactions
	  // with locally-owned nodes
	  force.PutScalar(0.0);

	  PdMaterialUtilities::computeInternalForceIsotropicElasticPlastic
        (x,
         y,
         weightedVolume,
         volume,
         dilatation,
         bondState,
         edpN,
         edpNP1,
         lambdaN,
         lambdaNP1,
         force.Values(),
         neighborhoodList,
         numOwnedPoints,
         m_bulkModulus,
         m_shearModulus,
         m_horizon,
         m_yieldStress);

}

