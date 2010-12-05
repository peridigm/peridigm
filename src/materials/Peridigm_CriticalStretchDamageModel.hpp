/*! \file Peridigm_CriticalStretchDamageModel.hpp */

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

#ifndef PERIDIGM_CRITICALSTRETCHDAMAGEMODEL_HPP
#define PERIDIGM_CRITICALSTRETCHDAMAGEMODEL_HPP

#include "Peridigm_DamageModel.hpp"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>

namespace PeridigmNS {

  //! Base class defining the Peridigm damage model interface.
  class CriticalStretchDamageModel : public DamageModel{

  public:
	
	//! Standard constructor.
	CriticalStretchDamageModel(const Teuchos::ParameterList& params);

	//! Destructor.
	virtual ~CriticalStretchDamageModel();

	//! Return name of material type
	virtual string Name() const { return("Critical Stretch"); }

    //! Returns a vector of field specs that specify the variables associated with the damage model
    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > VariableSpecs() const { return m_variableSpecs; }

	//! Initialize the damage model.
	virtual void
	initialize(const double dt,
               const int numOwnedPoints,
               const int* ownedIDs,
               const int* neighborhoodList,
               double* bondState,
               PeridigmNS::DataManager& dataManager) const ;

	//! Evaluate the damage
	virtual void
	computeDamage(const double dt,
                  const int numOwnedPoints,
                  const int* ownedIDs,
                  const int* neighborhoodList,
                  double* bondState,
                  PeridigmNS::DataManager& dataManager) const ;

  protected:

	//! Computes the distance between nodes (a1, a2, a3) and (b1, b2, b3).
	inline double distance(double a1, double a2, double a3,
						   double b1, double b2, double b3) const
	{
	  return ( sqrt( (a1-b1)*(a1-b1) + (a2-b2)*(a2-b2) + (a3-b3)*(a3-b3) ) );
	}

    Teuchos::RCP< std::vector<Field_NS::FieldSpec> > m_variableSpecs;

    double m_criticalStretch;
  };

}

#endif // PERIDIGM_CRITICALSTRETCHDAMAGEMODEL_HPP
