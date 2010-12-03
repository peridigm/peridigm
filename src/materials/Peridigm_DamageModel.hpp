/*! \file Peridigm_DamageModel.hpp */

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

#ifndef PERIDIGM_DAMAGEMODEL_HPP
#define PERIDIGM_DAMAGEMODEL_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include "Peridigm_DataManager.hpp"

namespace PeridigmNS {

  //! Base class defining the Peridigm damage model interface.
  class DamageModel{

  public:
	
	//! Standard constructor.
	DamageModel(const Teuchos::ParameterList & params){}

	//! Destructor.
	virtual ~DamageModel(){}

	//! Return name of material type
	virtual string Name() const = 0;

    //! Returns a vector of field specs that specify the variables associated with the damage model
    virtual Teuchos::RCP< std::vector<Field_NS::FieldSpec> > VariableSpecs() const = 0;

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

	//! Initialize the damage model.
	virtual void
	initialize(const Epetra_Vector& u,
               const Epetra_Vector& v,
               const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               double* bondState,
               PeridigmNS::DataManager& dataManager,
               Epetra_MultiVector& vectorConstitutiveData,
               Epetra_Vector& force) const {}

	//! Evaluate the damage
	virtual void
	computeDamage(const Epetra_Vector& u,
                  const Epetra_Vector& v,
                  const double dt,
                  const int numOwnedPoints,
                  const int* ownedIDs,
                  const int* neighborhoodList,
                  double* bondState,
                  PeridigmNS::DataManager& dataManager,
                  Epetra_MultiVector& vectorConstitutiveData,
                  Epetra_Vector& force) const = 0;

  private:
	
	//! Default constructor with no arguments, private to prevent use.
	DamageModel(){}
  };
}

#endif // PERIDIGM_DAMAGEMODEL_HPP
