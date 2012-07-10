//! \file Peridigm_LinearElasticIsotropicMaterial.hpp

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

#ifndef PERIDIGM_LINEARELASTICISOTROPICMATERIAL_HPP
#define PERIDIGM_LINEARELASTICISOTROPICMATERIAL_HPP

#include "Peridigm_Material.hpp"

namespace PeridigmNS {

  /*! \brief State-based peridynamic linear elastic isotropic material model.
   *
   * The state-based peridynamic linear elastic isotropic material is an
   * ordinary peridynamic material, meaning that the force resulting from the
   * bond between two nodes acts along the line connecting these nodes.
   * 
   * The magnitude of the pairwise force in a linear peridynamic solid is given by
   *
   * \f$ \underline{t} = \frac{\-3p}{m}\underline{\omega} \, \underline{x} + \frac{15\mu}{m} 
   *    \underline{\omega} \, \underline{e}^{d}, \qquad p = -k \theta \f$,
   *
   * where \f$ p \f$ is the peridynamic pressure, \f$ \mu \f$ and \f$ k \f$ are material 
   * constants (the shear modulus and bulk modulus, respectively), and the following 
   * definitions apply:
   *
   * \f$ \underline{\omega} \f$: Influence function.
   *
   * \f$ \underline{x} \f$: Reference position scalar state field, the distance between two
   * nodes in the reference configuration.
   *
   * \f$ \underline{y} \f$: Deformation scalar state field, the distance between two
   * nodes in the deformed configuration.
   *
   * \f$ \underline{e} \f$: Extension scalar state field. \f$ \underline{e} = 
   *    \underline{y} - \underline{x} \f$.
   *
   * \f$ m \f$: Weighted volume. \f$ m = \sum_{i=0}^{N} \underline{\omega}_{i} \, 
   *    \underline{x}_{i} \, \underline{x}_{i} \, \Delta V_{\mathbf{x}_{i}} \f$.
   *
   * \f$ N \f$: Number of cells in the neighborhood of \f$ x \f$.
   *
   * \f$ \theta \f$:  Dilatation. \f$ \theta =  \sum_{i=0}^{N} \frac{3}{m} \underline{\omega}_{i} 
   *    \, \underline{x}_{i} \, \underline{e}_{i} \, \Delta V_{\mathbf{x}_{i}} \f$.
   *
   * \f$ \underline{e}_{i} \f$:  Isotropic part of the extension. \f$ \underline{e}^{i} = 
   *    \frac{\theta \underline{x}}{3} \f$.
   *
   * \f$ \underline{e}_{d} \f$:  Deviatoric part of the extension. \f$ \underline{e}^{d} = 
   *    \underline{e} - \underline{e}^{i} \f$.
   */
  class LinearElasticIsotropicMaterial : public Material{
  public:

	//! Constructor.
	  LinearElasticIsotropicMaterial(const Teuchos::ParameterList & params);

	//! Destructor.
	virtual ~LinearElasticIsotropicMaterial();

	//! Return name of material type
	virtual string Name() const { return("Linear Elastic"); }

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

	//! Computes the dilatation.
	virtual void
	updateConstitutiveData(const double dt,
						   const int numOwnedPoints,
						   const int* ownedIDs,
						   const int* neighborhoodList,
                           PeridigmNS::DataManager& dataManager) const;

	//! Evaluate the forces on the cells.
	virtual void
	computeForce(const double dt,
				 const int numOwnedPoints,
				 const int* ownedIDs,
				 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager) const;

	//! Evaluate the jacobian
	virtual void
	computeJacobian(const double dt,
                    const int numOwnedPoints,
                    const int* ownedIDs,
                    const int* neighborhoodList,
                    PeridigmNS::DataManager& dataManager,
                    PeridigmNS::SerialMatrix& jacobian) const;

	//! Evaluate the jacobian via automatic differentiation
	virtual void
	computeAutomaticDifferentiationJacobian(const double dt,
                                            const int numOwnedPoints,
                                            const int* ownedIDs,
                                            const int* neighborhoodList,
                                            PeridigmNS::DataManager& dataManager,
                                            PeridigmNS::SerialMatrix& jacobian) const;

  protected:
	
	//! Computes the distance between nodes (a1, a2, a3) and (b1, b2, b3).
	inline double distance(double a1, double a2, double a3,
						   double b1, double b2, double b3) const
	{
	  return ( sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3) ) );
	}

    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > m_variableSpecs;

	// material parameters
	double m_bulkModulus;
	double m_shearModulus;
	double m_density;
    double m_horizon;
    bool m_applyAutomaticDifferentiationJacobian;
  };
}

#endif // PERIDIGM_LINEARELASTICISOTROPICMATERIAL_HPP
