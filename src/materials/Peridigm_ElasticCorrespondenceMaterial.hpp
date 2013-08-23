//! \file Peridigm_ElasticCorrespondenceMaterial.hpp

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

#ifndef PERIDIGM_ELASTICCORRESPONDENCEMATERIAL_HPP
#define PERIDIGM_ELASTICCORRESPONDENCEMATERIAL_HPP

#include "Peridigm_Material.hpp"

namespace PeridigmNS {

  class ElasticCorrespondenceMaterial : public Material{
  public:

	//! Constructor.
    ElasticCorrespondenceMaterial(const Teuchos::ParameterList & params);

    //! Destructor.
    virtual ~ElasticCorrespondenceMaterial();

    //! Return name of material type
    virtual std::string Name() const { return("Elastic Correspondence"); }

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

    //! Evaluate the internal force.
    virtual void
    computeForce(const double dt,
		 const int numOwnedPoints,
		 const int* ownedIDs,
		 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager) const;

  protected:

    //! Evaluate the hourglass force for suppression of zero-energy modes.
    virtual void
    computeHourglassForce(const double dt,
			  const int numOwnedPoints,
			  const int* ownedIDs,
			  const int* neighborhoodList,
			  PeridigmNS::DataManager& dataManager) const;

    // material parameters
    double m_bulkModulus;
    double m_shearModulus;
    double m_youngsModulus;
    double m_poissonsRatio;
    double m_density;
    double m_horizon;
    double m_hourglassCoefficient;

    // field spec ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_volumeFieldId;
    int m_modelCoordinatesFieldId;
    int m_coordinatesFieldId;
    int m_hourglassForceDensityFieldId;
    int m_forceDensityFieldId;
    int m_bondDamageFieldId;
    int m_deformationGradientXXFieldId;
    int m_deformationGradientXYFieldId;
    int m_deformationGradientXZFieldId;
    int m_deformationGradientYXFieldId;
    int m_deformationGradientYYFieldId;
    int m_deformationGradientYZFieldId;
    int m_deformationGradientZXFieldId;
    int m_deformationGradientZYFieldId;
    int m_deformationGradientZZFieldId;
    int m_shapeTensorInverseXXFieldId;
    int m_shapeTensorInverseXYFieldId;
    int m_shapeTensorInverseXZFieldId;
    int m_shapeTensorInverseYXFieldId;
    int m_shapeTensorInverseYYFieldId;
    int m_shapeTensorInverseYZFieldId;
    int m_shapeTensorInverseZXFieldId;
    int m_shapeTensorInverseZYFieldId;
    int m_shapeTensorInverseZZFieldId;
    int m_strainXXFieldId;
    int m_strainXYFieldId;
    int m_strainXZFieldId;
    int m_strainYXFieldId;
    int m_strainYYFieldId;
    int m_strainYZFieldId;
    int m_strainZXFieldId;
    int m_strainZYFieldId;
    int m_strainZZFieldId;
    int m_stressXXFieldId;
    int m_stressXYFieldId;
    int m_stressXZFieldId;
    int m_stressYXFieldId;
    int m_stressYYFieldId;
    int m_stressYZFieldId;
    int m_stressZXFieldId;
    int m_stressZYFieldId;
    int m_stressZZFieldId;
  };
}

#endif // PERIDIGM_ELASTICCORRESPONDENCEMATERIAL_HPP
