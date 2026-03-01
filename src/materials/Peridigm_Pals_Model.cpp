/*! \file Peridigm_Pals_Model.cpp */

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

#include "Peridigm_Pals_Model.hpp"
#include "Peridigm_Field.hpp"
#include "pals.h"
#include "material_utilities.h"
#include <Teuchos_Assert.hpp>
#include <Epetra_SerialComm.h>
#include <Sacado.hpp>
#include <string>
#include <cstdio>

using std::vector;

using MATERIAL_EVALUATION::PALS::NUM_LAGRANGE_MULTIPLIERS;

PeridigmNS::Pals_Model::Pals_Model(const Teuchos::ParameterList& params)
  : Material(params),
    m_bulkModulus(0.0), m_shearModulus(0.0), m_density(0.0), m_horizon(0.0),
    m_OMEGA_0(&PeridigmInfluenceFunction::one), m_SIGMA_0(&PeridigmInfluenceFunction::one),
    m_volumeFieldId(-1), m_weightedVolumeFieldId(-1), m_normalizedWeightedVolumeFieldId(-1),
    m_dilatationFieldId(-1), m_palsPressureFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_coordinatesFieldId(-1), m_forceDensityFieldId(-1),
    m_bondDamageFieldId(-1),
    num_lagrange_multipliers(NUM_LAGRANGE_MULTIPLIERS),
    m_dilatationNormalizationFieldId(-1), m_deviatoricNormalizationFieldId(-1),
    m_dilatationLagrangeMultiplersFieldIds(NUM_LAGRANGE_MULTIPLIERS),
    m_deviatoricLagrangeMultiplersFieldIds(NUM_LAGRANGE_MULTIPLIERS)
{
  //! \todo Add meaningful asserts on material properties.
  m_bulkModulus = calculateBulkModulus(params);
  m_shearModulus = calculateShearModulus(params);
  m_density = params.get<double>("Density");
  m_horizon = params.get<double>("Horizon");


  if(params.isParameter("Dilatation Influence Function")){
    string type = params.get<string>("Dilatation Influence Function");
    m_OMEGA_0=InfluenceFunction::getPredefinedInfluenceFunction(type);
  }

  if(params.isParameter("Deviatoric Influence Function")){
    string type = params.get<string>("Deviatoric Influence Function");
    m_SIGMA_0=InfluenceFunction::getPredefinedInfluenceFunction(type);
  }

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_volumeFieldId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_weightedVolumeFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Weighted_Volume");
  m_normalizedWeightedVolumeFieldId        = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Normalized_Weighted_Volume");
  m_dilatationFieldId            = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Dilatation");
  m_palsPressureFieldId          = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Pals_Pressure");
  m_modelCoordinatesFieldId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFieldId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_forceDensityFieldId          = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  m_bondDamageFieldId            = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Bond_Damage");
  m_dilatationNormalizationFieldId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Dilatation_Normalization");
  m_deviatoricNormalizationFieldId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deviatoric_Normalization");

  m_fieldIds.push_back(m_volumeFieldId);
  m_fieldIds.push_back(m_weightedVolumeFieldId);
  m_fieldIds.push_back(m_normalizedWeightedVolumeFieldId);
  m_fieldIds.push_back(m_dilatationFieldId);
  m_fieldIds.push_back(m_palsPressureFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_coordinatesFieldId);
  m_fieldIds.push_back(m_forceDensityFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
  m_fieldIds.push_back(m_dilatationNormalizationFieldId);
  m_fieldIds.push_back(m_deviatoricNormalizationFieldId);

  // Initialize field ids for array of Lagrange mulitpliers
  for(int i=0;i<num_lagrange_multipliers;i++){
    std::ostringstream dev,dil;
    dil << "Lagrange_Multiplier_Dilatation_"<< i+1;
    dev << "Lagrange_Multiplier_Deviatoric_"<< i+1;
    int m_dilId= fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, dil.str());
    int m_devId= fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, dev.str());
    m_fieldIds.push_back(m_dilId);
    m_fieldIds.push_back(m_devId);
    m_dilatationLagrangeMultiplersFieldIds[i]=m_dilId;
    m_deviatoricLagrangeMultiplersFieldIds[i]=m_devId;
  }

}

PeridigmNS::Pals_Model::~Pals_Model()
{
}

void
PeridigmNS::Pals_Model::
initialize(const double dt,
           const int numOwnedPoints,
           const int* ownedIDs,
           const int* neighborhoodList,
           PeridigmNS::DataManager& dataManager)
{

  // Extract pointers to the underlying data
  double *xOverlap,  *cellVolumeOverlap, *weightedVolume;
  double *omega_constants, *sigma_constants;
  double *normalized_weighted_volume;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolumeOverlap);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_normalizedWeightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&normalized_weighted_volume);
  dataManager.getData(m_dilatationNormalizationFieldId, PeridigmField::STEP_NONE)->ExtractView(&omega_constants);
  dataManager.getData(m_deviatoricNormalizationFieldId, PeridigmField::STEP_NONE)->ExtractView(&sigma_constants);

  // namespace PALS
  using namespace MATERIAL_EVALUATION::PALS;
  {
    vector<double *> omega_multipliers(num_lagrange_multipliers), sigma_multipliers(num_lagrange_multipliers);
    for(int i=0;i<num_lagrange_multipliers;i++){
      double *dil, *dev;
      dataManager.getData(m_dilatationLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dil);
      dataManager.getData(m_deviatoricLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dev);
      omega_multipliers[i]=dil;
      sigma_multipliers[i]=dev;
    }

    compute_lagrange_multipliers
    (
      xOverlap,
      cellVolumeOverlap,
      numOwnedPoints,
      neighborhoodList,
      m_horizon,
      omega_multipliers,
      omega_constants,
      sigma_multipliers,
      sigma_constants,
      m_OMEGA_0,
      m_SIGMA_0
    );
  }


  {
    vector<const double *> sigma_multipliers(num_lagrange_multipliers);
    for(int i=0;i<num_lagrange_multipliers;i++){
      double *dev;
      dataManager.getData(m_deviatoricLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dev);
      sigma_multipliers[i]=dev;
    }
    /*
     * DEBUGGING
     * Compute normalized weighted volume: all values should be 3.0
    computeNormalizedWeightedVolume
    (
      xOverlap,
      cellVolumeOverlap,
      omega_constants,
      bondDamage,
      normalized_weighted_volume,
      numOwnedPoints,
      neighborhoodList,
      m_horizon,
      m_SIGMA_0
    );
    */
    /*
     * Weighted volume with influence function $\sigma$ rather than $\omega$
     */
    computeWeightedVolume
    (
      xOverlap,
      cellVolumeOverlap,
      sigma_multipliers,
      sigma_constants,
      weightedVolume,
      numOwnedPoints,
      neighborhoodList,
      m_horizon,
      m_SIGMA_0
    );
  }
}


void
PeridigmNS::Pals_Model::computeForce(const double dt,
                                          const int numOwnedPoints,
                                          const int* ownedIDs,
                                          const int* neighborhoodList,
                                          PeridigmNS::DataManager& dataManager) const
{
  // Zero out the forces
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->PutScalar(0.0);

  // Extract pointers to the underlying data
  double *x, *y, *cellVolume, *pals_pressure, *force, *weightedVolume;
  double *omega_constants, *sigma_constants;
  double *dilatation;
  vector<const double *> omega_multipliers(num_lagrange_multipliers), sigma_multipliers(num_lagrange_multipliers);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&y);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&cellVolume);
  dataManager.getData(m_weightedVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&weightedVolume);
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(m_palsPressureFieldId, PeridigmField::STEP_NP1)->ExtractView(&pals_pressure);
  dataManager.getData(m_forceDensityFieldId, PeridigmField::STEP_NP1)->ExtractView(&force);
  dataManager.getData(m_dilatationNormalizationFieldId, PeridigmField::STEP_NONE)->ExtractView(&omega_constants);
  dataManager.getData(m_deviatoricNormalizationFieldId, PeridigmField::STEP_NONE)->ExtractView(&sigma_constants);

  for(int i=0;i<num_lagrange_multipliers;i++){
    double *dil, *dev;
    dataManager.getData(m_dilatationLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dil);
    dataManager.getData(m_deviatoricLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dev);
    omega_multipliers[i]=dil;
    sigma_multipliers[i]=dev;
  }

  // namespace PALS
  using namespace MATERIAL_EVALUATION::PALS;
  computeDilatationAndPalsPressure
  (
    x,
    y,
    cellVolume,
    omega_multipliers,
    omega_constants,
    sigma_multipliers,
    sigma_constants,
    weightedVolume,
    dilatation,
    pals_pressure,
    neighborhoodList,
    numOwnedPoints,
    m_bulkModulus,
    m_shearModulus,
    m_horizon,
    m_OMEGA_0,
    m_SIGMA_0
  );

  computeInternalForcePals
  (
    x,
    y,
    cellVolume,
    omega_multipliers,
    omega_constants,
    sigma_multipliers,
    sigma_constants,
    dilatation,
    pals_pressure,
    force,
    neighborhoodList,
    numOwnedPoints,
    m_bulkModulus,
    m_shearModulus,
    m_horizon,
    m_OMEGA_0,
    m_SIGMA_0
  );


}

void
PeridigmNS::Pals_Model::computeStoredElasticEnergyDensity(const double dt,
                                                               const int numOwnedPoints,
                                                               const int* ownedIDs,
                                                               const int* neighborhoodList,
                                                               PeridigmNS::DataManager& dataManager) const
{
  // namespace PALS
  using namespace MATERIAL_EVALUATION::PALS;

  double K=m_bulkModulus;
  double MU=m_shearModulus;
  double horizon=m_horizon;

  // This function is intended to be called from a compute class.
  // The compute class should have already created the Stored_Elastic_Energy_Density field id.
  int storedElasticEnergyDensityFieldId = PeridigmNS::FieldManager::self().getFieldId("Stored_Elastic_Energy_Density");

  // Extract pointers to the underlying data
  double *xOverlap, *yOverlap, *volumeOverlap;
  double *omega_constants, *sigma_constants;
  double *dilatation, *energy_density;
  vector<const double *> omega_multipliers(num_lagrange_multipliers), sigma_multipliers(num_lagrange_multipliers);
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&xOverlap);
  dataManager.getData(m_coordinatesFieldId, PeridigmField::STEP_NP1)->ExtractView(&yOverlap);
  dataManager.getData(m_volumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&volumeOverlap);
  dataManager.getData(m_dilatationFieldId, PeridigmField::STEP_NP1)->ExtractView(&dilatation);
  dataManager.getData(m_dilatationNormalizationFieldId, PeridigmField::STEP_NONE)->ExtractView(&omega_constants);
  dataManager.getData(m_deviatoricNormalizationFieldId, PeridigmField::STEP_NONE)->ExtractView(&sigma_constants);
  dataManager.getData(storedElasticEnergyDensityFieldId , PeridigmField::STEP_NONE)->ExtractView(&energy_density);
  for(int i=0;i<num_lagrange_multipliers;i++){
    double *dil, *dev;
    dataManager.getData(m_dilatationLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dil);
    dataManager.getData(m_deviatoricLagrangeMultiplersFieldIds[i],PeridigmField::STEP_NONE)->ExtractView(&dev);
    omega_multipliers[i]=dil;
    sigma_multipliers[i]=dev;
  }

  double bond[3];
  const double *xOwned = xOverlap;
  const double *yOwned = yOverlap;
  const double *oc=omega_constants;
  const double *sc=sigma_constants;
  const double *theta=dilatation;
  const int *neighPtr = neighborhoodList;
  double cell_volume;
  for(int q=0; q<numOwnedPoints; q++, xOwned+=3, yOwned+=3, oc++, sc++, theta++){

    // Collect computed Lagrange multipliers for this point 'q'
    double sigma_X[NUM_LAGRANGE_MULTIPLIERS];
    for(int i=0;i<NUM_LAGRANGE_MULTIPLIERS;i++){
      sigma_X[i]=sigma_multipliers[i][q];
    }

    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const double *Y = yOwned;
    double sum=0.0;
    for(int n=0;n<numNeigh;n++,neighPtr++){
      int localId = *neighPtr;
      cell_volume = volumeOverlap[localId];
      const double *XP = &xOverlap[3*localId];
      const double *YP = &yOverlap[3*localId];
      bond[0]=XP[0]-X[0];
      bond[1]=XP[1]-X[1];
      bond[2]=XP[2]-X[2];
      double a = bond[0];
      double b = bond[1];
      double c = bond[2];
      double xi = sqrt(a*a+b*b+c*c);
      a = YP[0]-Y[0];
      b = YP[1]-Y[1];
      c = YP[2]-Y[2];
      double dY = sqrt(a*a+b*b+c*c);
      double e = dY-xi;
      double eps=e-xi*(*theta)/3.0;
      pals_influence<deviatoric_influence> SIGMA(m_SIGMA_0,*sc,sigma_X);
      double sigma = SIGMA(bond,horizon);
      sum += sigma * eps * eps * cell_volume;
    }
    double d=*theta;
    /*
     * Deviatoric term does not include factor of 1/2
     */
    energy_density[q] = 0.5 * K * d * d + MU * sum;
  }
}
