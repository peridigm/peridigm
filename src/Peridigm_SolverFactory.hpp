/*! \file Peridigm_SolverFactory.hpp */

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

namespace Peridigm {

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
