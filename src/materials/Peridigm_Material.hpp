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

	//! Returns the number of scalar constitutive variables used by the material model.
	virtual int NumScalarConstitutiveVariables() const = 0;

	//! Returns the number of vector constitutive variables used by the material model.
	virtual int NumVectorConstitutiveVariables() const = 0;

	//! Returns the number of scalar constitutive bond variables used by the material model.
	virtual int NumBondConstitutiveVariables() const = 0;

	//! Returns the name of the scalar constitutive variable at position pos.
	virtual const string & ScalarConstitutiveVariableName(int pos) const = 0;

	//! Returns the name of the scalar constitutive variable at position pos.
	virtual const string & VectorConstitutiveVariableName(int pos) const = 0;

	//! Returns the name of the bond constitutive variable at position pos.
	virtual const string & BondConstitutiveVariableName(int pos) const = 0;

	//! Initialize the material model.
	virtual void
	initialize(const Epetra_Vector& x,
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
               Epetra_Vector& force) const {}

	//! Update the constitutive data based on the current configuration.
	virtual void
	updateConstitutiveData(const Epetra_Vector& x,
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
						   Epetra_Vector& force) const = 0;

	//! Evaluate the forces on the cells
	virtual void
	computeForce(const Epetra_Vector& x,
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
				 Epetra_Vector& force) const = 0;

  private:
	
	//! Default constructor with no arguments, private to prevent use.
	Material(){}
  };
}

#endif // PERIDIGM_MATERIAL_HPP
