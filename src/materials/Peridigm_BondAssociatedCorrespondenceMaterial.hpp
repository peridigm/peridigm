//! \file Peridigm_BondAssociatedCorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_BONDASSOCIATEDCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_BONDASSOCIATEDCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_Material.hpp"
#include "Peridigm_InfluenceFunction.hpp"

namespace PeridigmNS {

  class BondAssociatedCorrespondenceMaterial : public Material{
  public:

    //! Constructor.
    BondAssociatedCorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~BondAssociatedCorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Bond Associated Correspondence Base Class"); }

    //! Returns the density of the material.
    virtual double Density() const { return m_density; }

    //! Returns the bulk modulus of the material.
    virtual double BulkModulus() const { return m_bulkModulus; }

    //! Returns the shear modulus of the material.
    virtual double ShearModulus() const { return m_shearModulus; }

    //! Returns a vector of field IDs corresponding to the variables associated with the material.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Returns a vector of field IDs that need to be synchronized across block boundaries and MPI boundaries after initialize().
    virtual std::vector<int> FieldIdsForSynchronizationAfterInitialize() const {
      std::vector<int> fieldIds;
      fieldIds.push_back(m_weightedVolumeFieldId);
      return fieldIds;
    }

    //! Returns a vector of field IDs that need to be synchronized across block boundaries and MPI boundaries after precompute().
    virtual std::vector<int> FieldIdsForSynchronizationAfterPrecompute() const {
      std::vector<int> fieldIds;
      fieldIds.push_back(m_deformationGradientXFieldId);
      fieldIds.push_back(m_deformationGradientYFieldId);
      fieldIds.push_back(m_deformationGradientZFieldId);
      fieldIds.push_back(m_deformationGradientDotXFieldId);
      fieldIds.push_back(m_deformationGradientDotYFieldId);
      fieldIds.push_back(m_deformationGradientDotZFieldId);
      return fieldIds;
    }

    //! Initialize the material model.
    virtual void initialize(const double dt,
                            const int numOwnedPoints,
                            const int* ownedIDs,
                            const int* neighborhoodList,
                            PeridigmNS::DataManager& dataManager);

    //! Evaluate the Cauchy stress (pure virtual function, must be implemented by derived correspondence material models).
    virtual void computeCauchyStress(const double dt,
                                     const int numOwnedPoints,
                                     const int* neighborhoodList,
                                     PeridigmNS::DataManager& dataManager) const = 0;

    //! Evaluate the internal force.
    virtual void computeForce(const double dt,
                              const int numOwnedPoints,
                              const int* ownedIDs,
                              const int* neighborhoodList,
                              PeridigmNS::DataManager& dataManager) const;

    //! Evaluate the node-level fields for symmetric bond-associated models
    virtual void 
    precompute(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) const;

  protected:

    // material parameters
    double m_bulkModulus;
    double m_shearModulus;
    double m_density;
    int m_accuracyOrder;
    PeridigmNS::InfluenceFunction::functionPointer m_OMEGA;

    // field spec ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_horizonFieldId;
    int m_volumeFieldId;
    int m_modelCoordinatesFieldId;
    int m_displacementsFieldId;
    int m_coordinatesFieldId;
    int m_velocitiesFieldId;
    int m_forceDensityFieldId;
    int m_damageFieldId;
    int m_bondDamageFieldId;
    int m_influenceStateFieldId;
    int m_gradientWeightXFieldId;
    int m_gradientWeightYFieldId;
    int m_gradientWeightZFieldId;
    int m_gradientWeightEvaluationFlagFieldId;
    int m_deformationGradientXFieldId;
    int m_deformationGradientYFieldId;
    int m_deformationGradientZFieldId;
    int m_deformationGradientDotXFieldId;
    int m_deformationGradientDotYFieldId;
    int m_deformationGradientDotZFieldId;
    int m_greenLagrangeStrainFieldId;
    int m_leftStretchTensorFieldId;
    int m_rotationTensorFieldId;
    int m_unrotatedCauchyStressFieldId;
    int m_cauchyStressFieldId;
    int m_unrotatedRateOfDeformationFieldId;
    int m_weightedVolumeFieldId;

    int m_bondLevelLeftStretchTensorXXFieldId;
    int m_bondLevelLeftStretchTensorXYFieldId;
    int m_bondLevelLeftStretchTensorXZFieldId;
    int m_bondLevelLeftStretchTensorYXFieldId;
    int m_bondLevelLeftStretchTensorYYFieldId;
    int m_bondLevelLeftStretchTensorYZFieldId;
    int m_bondLevelLeftStretchTensorZXFieldId;
    int m_bondLevelLeftStretchTensorZYFieldId;
    int m_bondLevelLeftStretchTensorZZFieldId;
    int m_bondLevelRotationTensorXXFieldId;
    int m_bondLevelRotationTensorXYFieldId;
    int m_bondLevelRotationTensorXZFieldId;
    int m_bondLevelRotationTensorYXFieldId;
    int m_bondLevelRotationTensorYYFieldId;
    int m_bondLevelRotationTensorYZFieldId;
    int m_bondLevelRotationTensorZXFieldId;
    int m_bondLevelRotationTensorZYFieldId;
    int m_bondLevelRotationTensorZZFieldId;
    int m_bondLevelUnrotatedCauchyStressXXFieldId;
    int m_bondLevelUnrotatedCauchyStressXYFieldId;
    int m_bondLevelUnrotatedCauchyStressXZFieldId;
    int m_bondLevelUnrotatedCauchyStressYXFieldId;
    int m_bondLevelUnrotatedCauchyStressYYFieldId;
    int m_bondLevelUnrotatedCauchyStressYZFieldId;
    int m_bondLevelUnrotatedCauchyStressZXFieldId;
    int m_bondLevelUnrotatedCauchyStressZYFieldId;
    int m_bondLevelUnrotatedCauchyStressZZFieldId;
    int m_bondLevelCauchyStressXXFieldId;
    int m_bondLevelCauchyStressXYFieldId;
    int m_bondLevelCauchyStressXZFieldId;
    int m_bondLevelCauchyStressYXFieldId;
    int m_bondLevelCauchyStressYYFieldId;
    int m_bondLevelCauchyStressYZFieldId;
    int m_bondLevelCauchyStressZXFieldId;
    int m_bondLevelCauchyStressZYFieldId;
    int m_bondLevelCauchyStressZZFieldId;
    int m_bondLevelUnrotatedRateOfDeformationXXFieldId;
    int m_bondLevelUnrotatedRateOfDeformationXYFieldId;
    int m_bondLevelUnrotatedRateOfDeformationXZFieldId;
    int m_bondLevelUnrotatedRateOfDeformationYXFieldId;
    int m_bondLevelUnrotatedRateOfDeformationYYFieldId;
    int m_bondLevelUnrotatedRateOfDeformationYZFieldId;
    int m_bondLevelUnrotatedRateOfDeformationZXFieldId;
    int m_bondLevelUnrotatedRateOfDeformationZYFieldId;
    int m_bondLevelUnrotatedRateOfDeformationZZFieldId;
    int m_bondLevelPiolaStressXXFieldId;
    int m_bondLevelPiolaStressXYFieldId;
    int m_bondLevelPiolaStressXZFieldId;
    int m_bondLevelPiolaStressYXFieldId;
    int m_bondLevelPiolaStressYYFieldId;
    int m_bondLevelPiolaStressYZFieldId;
    int m_bondLevelPiolaStressZXFieldId;
    int m_bondLevelPiolaStressZYFieldId;
    int m_bondLevelPiolaStressZZFieldId;
    int m_bondLevelDeformationGradientInvXXFieldId;
    int m_bondLevelDeformationGradientInvXYFieldId;
    int m_bondLevelDeformationGradientInvXZFieldId;
    int m_bondLevelDeformationGradientInvYXFieldId;
    int m_bondLevelDeformationGradientInvYYFieldId;
    int m_bondLevelDeformationGradientInvYZFieldId;
    int m_bondLevelDeformationGradientInvZXFieldId;
    int m_bondLevelDeformationGradientInvZYFieldId;
    int m_bondLevelDeformationGradientInvZZFieldId;
    int m_bondLevelJacobianDeterminantFieldId;
    int m_stressIntegralFieldId;
  };
}

#endif // PERIDIGM_BONDASSOCIATEDCORRESPONDENCEMATERIAL_HPP
