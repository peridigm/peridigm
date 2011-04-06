/*! \file Peridigm_SolverFactory.hpp */

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

#ifndef PERIDIGM_SOLVERFACTORY_HPP
#define PERIDIGM_SOLVERFACTORY_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Rythmos_IntegrationObserverBase.hpp>

#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include <Epetra_SerialComm.h>
#endif

#include "Peridigm_ModelEvaluator.hpp"

namespace PeridigmNS {

  /*!
   * \brief A factory class to instantiate AbstractSolver objects.
   */
  class SolverFactory {
  public:

    //! Default constructor.
    SolverFactory(const std::string inputFile, const MPI_Comm& appComm);

    //! Destructor
    virtual ~SolverFactory() {}

	/*! \brief Creates an ENAT::RythmosSolver and returns an RCP to it.
	 *  This is the ModelEvaluator that allow external drivers such as
	 *  Dakota to communicate with the code. */
     virtual Teuchos::RCP<EpetraExt::ModelEvaluator> create();

  private:

    //! Private function to set default problem parameter values in lieu of InArgs.
    void setProblemParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> appParams_);
    
    //! Private function to set default solver parameter values in lieu of InArgs.
    void setSolverParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> appParams_);

    //! Private copy constructory to prohibit copying.
    SolverFactory(const SolverFactory&);

    //! Private assignment operator to prohibit copying.
    SolverFactory& operator=(const SolverFactory&);

  protected:

    typedef double Scalar;

    //! Parameter list specifying what solver to create
    Teuchos::RCP<Teuchos::ParameterList> appParams;
    Teuchos::RCP<Epetra_Comm> Comm;
    Teuchos::RCP<EpetraExt::ModelEvaluator> model;
  };

}

#endif // PERIDIGM_SOLVERFACTORY_HPP
