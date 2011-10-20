/*! \file Peridigm_DiscretizationFactory.hpp */

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

#include <Teuchos_Assert.hpp>
#include "Peridigm_DiscretizationFactory.hpp"
#include "Peridigm_PdQuickGridDiscretization.hpp"
#ifdef PERIDIGM_STK
#include "Peridigm_STKDiscretization.hpp"
#endif

PeridigmNS::DiscretizationFactory::DiscretizationFactory(const Teuchos::RCP<Teuchos::ParameterList>& discParams_) :
  discParams(discParams_)
{
  // check to see if a test configuration has been specified
  // or if a mesh file has been supplied
  if(!discParams->isParameter("Type")){
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
					   "Discretization not specified, \"Exodus\" or \"PdQuickGrid\" required.");
  }

}

Teuchos::RCP<PeridigmNS::AbstractDiscretization>
PeridigmNS::DiscretizationFactory::create(const Teuchos::RCP<const Epetra_Comm>& epetra_comm)
{
  Teuchos::RCP<PeridigmNS::AbstractDiscretization> discretization;

  string type = discParams->get<string>("Type");

  if(type == "PdQuickGrid"){
	discretization = Teuchos::rcp(new PeridigmNS::PdQuickGridDiscretization(epetra_comm, discParams));
  }
  else if(type == "Exodus"){
#ifdef PERIDIGM_STK
	discretization = Teuchos::rcp(new PeridigmNS::STKDiscretization(epetra_comm, discParams));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
		       "**** Exodus file support requires Netcdf and Trilinos STK package (rebuild with Netcdf and STK).\n");
#endif
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
		       "**** Invalid discretization type.  Valid types are \"Exodus\" and \"PdQuickGrid\".\n");
  }
 
  return discretization;
}

Teuchos::RCP<PeridigmNS::AbstractDiscretization>
PeridigmNS::DiscretizationFactory::create(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                          const Teuchos::RCP<const QUICKGRID::QuickGridData>& decomp)
{
  Teuchos::RCP<PeridigmNS::AbstractDiscretization> discretization;

  string type = discParams->get<string>("Type");

  if(type == "PdQuickGrid"){
	discretization = Teuchos::rcp(new PeridigmNS::PdQuickGridDiscretization(epetra_comm, decomp));
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
					   "Invalid type for construction of discretization from PdGridData object.  Valid type is \"PdQuickGrid\"");
  }
 
  return discretization;
}
