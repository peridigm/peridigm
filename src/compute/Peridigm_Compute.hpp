/*! \file Peridigm_Compute.hpp */

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

#ifndef PERIDIGM_COMPUTE_HPP
#define PERIDIGM_COMPUTE_HPP

#include <Teuchos_RCP.hpp>

#include "Peridigm_DataManager.hpp"

namespace PeridigmNS {

  //! Base class defining the Peridigm compute class interface.
  class Compute{

  public:
	
  //! Default constructor.
  Compute(){};

  //! Destructor.
  virtual ~Compute(){}

  //! Returns the fieldspecs computed by this class
  virtual std::vector<Field_NS::FieldSpec> getFieldSpecs() const = 0;

  //! Perform computation
  virtual int compute(Teuchos::RCP<PeridigmNS::DataManager> dataManager) const = 0;

  private:

  //! Copy constructor.
  Compute( const Compute& C );

  //! Assignment operator.
  Compute& operator=( const Compute& C );
	
  };
}

#endif // PERIDIGM_MATERIAL_HPP
