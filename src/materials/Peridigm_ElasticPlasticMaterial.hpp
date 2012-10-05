/*! \file Peridigm_ElasticPlasticMaterial.hpp */

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

#ifndef PERIDIGM_ELASTICPLASTICMATERIAL_HPP_
#define PERIDIGM_ELASTICPLASTICMATERIAL_HPP_

#include "Peridigm_Material.hpp"

namespace PeridigmNS {

  class ElasticPlasticMaterial : public Material {
  public:

	//! Constructor.
	ElasticPlasticMaterial(const Teuchos::ParameterList & params);

	//! Destructor.
	virtual ~ElasticPlasticMaterial();

	//! Return name of material type
	virtual string Name() const {return("Elastic Plastic");}

	//! Returns the density of the material.
	virtual double Density() const { return m_density; }

	//! Returns the bulk modulus of the material.
	virtual double BulkModulus() const { return m_bulkModulus; }
	        
	//! Returns the shear modulus of the material.
	virtual double ShearModulus() const { return m_shearModulus; }

	//! Returns the horizon.
	virtual double Horizon() const { return m_horizon; }

    //! Returns a vector of field specs that specify the variables associated with the material
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > VariableSpecs() const { return m_variableSpecs; }

	//! Initialized data containers and computes weighted volume.
	virtual void
	initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) const;

	//! Evaluate the internal force.
	virtual void
	computeForce(const double dt,
				 const int numOwnedPoints,
				 const int* ownedIDs,
				 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager) const;

	//! Evaluate the jacobian.
	virtual void
	computeJacobian(const double dt,
                    const int numOwnedPoints,
                    const int* ownedIDs,
                    const int* neighborhoodList,
                    PeridigmNS::DataManager& dataManager,
                    PeridigmNS::SerialMatrix& jacobian) const;

	//! Evaluate the jacobian via automatic differentiation.
	virtual void
	computeAutomaticDifferentiationJacobian(const double dt,
                                            const int numOwnedPoints,
                                            const int* ownedIDs,
                                            const int* neighborhoodList,
                                            PeridigmNS::DataManager& dataManager,
                                            PeridigmNS::SerialMatrix& jacobian) const;

  protected:

    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > m_variableSpecs;

	// material parameters
	double m_bulkModulus;
	double m_shearModulus;
	double m_horizon;
	double m_density;
	double m_yieldStress;
	double m_thickness;
    bool m_applyShearCorrectionFactor;
    bool m_disablePlasticity;
    bool m_applyAutomaticDifferentiationJacobian;
    bool m_isPlanarProblem;
  };
}



#endif /* PERIDIGM_ELASTICPLASTICMATERIAL_HPP_ */
