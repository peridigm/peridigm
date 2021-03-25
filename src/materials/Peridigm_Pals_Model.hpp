/*
 * Peridigm_PALS_Model.hpp
 *
 *  Created on: Aug 21, 2013
 *      Author: jamitch
 */

#ifndef PERIDIGM_PALS_MODEL_HPP_
#define PERIDIGM_PALS_MODEL_HPP_

#include "Peridigm_Material.hpp"
#include "Peridigm_InfluenceFunction.hpp"

namespace PeridigmNS {

class Pals_Model : public Material {
public:

  typedef PeridigmNS::InfluenceFunction::functionPointer FunctionPointer;

  //! Constructor.
  Pals_Model(const Teuchos::ParameterList & params);

  //! Destructor.
  virtual ~Pals_Model();

  //! Return name of material type
  virtual std::string Name() const { return("Pals"); }

  //! Returns the density of the material.
  virtual double Density() const { return m_density; }

  //! Returns the bulk modulus of the material.
  virtual double BulkModulus() const { return m_bulkModulus; }

  //! Returns the shear modulus of the material.
  virtual double ShearModulus() const { return m_shearModulus; }

  //! Returns the horizon.
  virtual double Horizon() const { return m_horizon; }

  //! Returns a vector of field IDs corresponding to the variables associated with the material.
  virtual std::vector<int> FieldIds() const { return m_fieldIds; }

  //! Initialized data containers and computes weighted volume.
  virtual void
  initialize(const double dt,
             const int numOwnedPoints,
             const int* ownedIDs,
             const int* neighborhoodList,
             PeridigmNS::DataManager& dataManager);

  //! Evaluate the internal force.
  virtual void
  computeForce(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) const;

    //! Compute stored elastic density energy.
    virtual void
    computeStoredElasticEnergyDensity(const double dt,
                                      const int numOwnedPoints,
                                      const int* ownedIDs,
                                      const int* neighborhoodList,
                                      PeridigmNS::DataManager& dataManager) const;

 protected:


   // material parameters
   double m_bulkModulus;
   double m_shearModulus;
   double m_density;
   double m_horizon;

   // Influence functions
   FunctionPointer m_OMEGA_0;
   FunctionPointer m_SIGMA_0;

   // field spec ids for all relevant data
   std::vector<int> m_fieldIds;
   int m_volumeFieldId;
   int m_weightedVolumeFieldId;
   int m_normalizedWeightedVolumeFieldId;
   int m_dilatationFieldId;
   int m_palsPressureFieldId;
   int m_modelCoordinatesFieldId;
   int m_coordinatesFieldId;
   int m_forceDensityFieldId;
   int m_bondDamageFieldId;

   const int num_lagrange_multipliers;
   int m_dilatationNormalizationFieldId;
   int m_deviatoricNormalizationFieldId;
   std::vector<int> m_dilatationLagrangeMultiplersFieldIds;
   std::vector<int> m_deviatoricLagrangeMultiplersFieldIds;

 };


}


#endif /* PERIDIGM_PALS_MODEL_HPP_ */
