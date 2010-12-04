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

#ifndef PERIDIGM_ISOTROPICELASTICPLASTICMATERIAL_HPP_
#define PERIDIGM_ISOTROPICELASTICPLASTICMATERIAL_HPP_

#include "Peridigm_Material.hpp"
#include "Peridigm_DamageModel.hpp"

namespace PeridigmNS {

  class IsotropicElasticPlasticMaterial : public Material {
  public:

	//! Constructor.
	IsotropicElasticPlasticMaterial(const Teuchos::ParameterList & params);

	//! Destructor.
	virtual ~IsotropicElasticPlasticMaterial();

	//! Return name of material type
	virtual string Name() const {return("Elastic Plastic");}

	//! Returns the density of the material.
	virtual double Density() const { return m_density; }

    //! Returns a vector of field specs that specify the variables associated with the material
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > VariableSpecs() const { return m_variableSpecs; }

	//! Initialized data containers and computes weighted volume.
	virtual void
	initialize(const Epetra_Vector& u,
               const Epetra_Vector& v,
               const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               double* bondState,
               PeridigmNS::DataManager& dataManager,
               Epetra_Vector& force) const;

	//! Computes the dilatation.
	virtual void
	updateConstitutiveData(const Epetra_Vector& u,
						   const Epetra_Vector& v,
						   const double dt,
						   const int numOwnedPoints,
						   const int* ownedIDs,
						   const int* neighborhoodList,
						   double* bondState,
                           PeridigmNS::DataManager& dataManager,
						   Epetra_Vector& force) const;

	//! Evaluate the forces on the cells.
	virtual void
	computeForce(const Epetra_Vector& u,
				 const Epetra_Vector& v,
				 const double dt,
				 const int numOwnedPoints,
				 const int* ownedIDs,
				 const int* neighborhoodList,
				 double* bondState,
                 PeridigmNS::DataManager& dataManager,
				 Epetra_Vector& force) const;

  protected:

    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > m_variableSpecs;

	// material parameters
	double m_bulkModulus;
	double m_shearModulus;
	double m_horizon;
	double m_density;
	double m_yieldStress;

    // damage model
    Teuchos::RCP<DamageModel> m_damageModel;

  };
}



#endif /* PERIDIGM_ISOTROPICELASTICPLASTICMATERIAL_HPP_ */
