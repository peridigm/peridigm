/*! \file Peridigm_ModelEvaluator.hpp */

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

#include <Teuchos_TestForException.hpp>
#include "Peridigm_DiscretizationFactory.hpp"
#include "Peridigm_PdQuickGridDiscretization.hpp"

PeridigmNS::DiscretizationFactory::DiscretizationFactory(const Teuchos::RCP<Teuchos::ParameterList>& discParams_) :
  discParams(discParams_)
{
  // check to see if a test configuration has been specified
  // or if a mesh file has been supplied
  if(!discParams->isParameter("Type")){
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
					   "Discretization not specified, \"Input File\" or \"PdQuickGrid\" required.");
  }

}

Teuchos::RCP<PeridigmNS::AbstractDiscretization>
PeridigmNS::DiscretizationFactory::create(const Teuchos::RCP<const Epetra_Comm>& epetra_comm)
{
  Teuchos::RCP<PeridigmNS::AbstractDiscretization> discretization;

  string type = discParams->get<string>("Type");

  // currently there is only one available discretization, the PdQuickGridDiscretization
  if(type == "Input File"){
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
					   "Software is currently unable to read mesh file, please try again later.");
  }
  else if(type == "PdQuickGrid"){
	discretization = Teuchos::rcp(new PeridigmNS::PdQuickGridDiscretization(epetra_comm, discParams));
  }
  else{
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
					   "Invalid discretization type.  Valid types are \"Input File\" and \"PdQuickGrid\"");
  }
 
  return discretization;
}
