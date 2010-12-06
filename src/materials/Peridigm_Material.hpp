/*! \file Peridigm_Material.hpp */

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

#ifndef PERIDIGM_MATERIAL_HPP
#define PERIDIGM_MATERIAL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <map>
#include "Peridigm_DataManager.hpp"

namespace PeridigmNS {

  //! Base class defining the Peridigm material model interface.
  class Material{

  public:
	
	//! Standard constructor.
	Material(const Teuchos::ParameterList & params){}

	//! Destructor.
	virtual ~Material(){}

	//! Return name of material type
	virtual string Name() const = 0;

	//! Returns the density of the material.
	virtual double Density() const = 0;

    //! Returns a vector of field specs that specify the variables associated with the material
    virtual Teuchos::RCP< std::vector<Field_NS::FieldSpec> > VariableSpecs() const = 0;

	//! Initialize the material model.
	virtual void
	initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               PeridigmNS::DataManager& dataManager) const {}

	//! Update the constitutive data based on the current configuration.
	virtual void
	updateConstitutiveData(const double dt,
						   const int numOwnedPoints,
						   const int* ownedIDs,
						   const int* neighborhoodList,
                           PeridigmNS::DataManager& dataManager) const = 0;

	//! Evaluate the forces on the cells
	virtual void
	computeForce(const double dt,
				 const int numOwnedPoints,
				 const int* ownedIDs,
				 const int* neighborhoodList,
                 PeridigmNS::DataManager& dataManager) const = 0;

  private:
	
	//! Default constructor with no arguments, private to prevent use.
	Material(){}
  };
}

#endif // PERIDIGM_MATERIAL_HPP
