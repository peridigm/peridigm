/*! \file Peridigm_Factory.hpp */

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

#ifndef PERIDIGM_FACTORY_HPP
#define PERIDIGM_FACTORY_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
typedef int MPI_Comm;
#define MPI_COMM_WORLD 1
#include <Epetra_SerialComm.h>
#endif

#include "Peridigm.hpp"

namespace PeridigmNS {

  /*!
   * \brief A factory class to instantiate Peridigm object
   */
  class PeridigmFactory {
  public:

    //! Default constructor.
    PeridigmFactory();

    //! Destructor
    virtual ~PeridigmFactory() {}

    virtual Teuchos::RCP<Peridigm::Peridigm> create(const std::string inputFile, const MPI_Comm& solverComm);

  private:

    //! Private function to set default problem parameter values in lieu of InArgs.
    void setProblemParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_);
    
    //! Private function to set default solver parameter values in lieu of InArgs.
    void setSolverParamDefaults(Teuchos::Ptr<Teuchos::ParameterList> peridigmParams_);

    //! Private copy constructory to prohibit copying.
    PeridigmFactory(const PeridigmFactory&);

    //! Private assignment operator to prohibit copying.
    PeridigmFactory& operator=(const PeridigmFactory&);

  protected:

  };

}

#endif // PERIDIGM_FACTORY_HPP
