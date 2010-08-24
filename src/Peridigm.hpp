/*! \file Peridigm.hpp */

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

#ifndef PERIDIGM_HPP
#define PERIDIGM_HPP

#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include "Peridigm_AbstractDiscretization.hpp"

namespace Peridigm {

  class Peridigm {

  public:

    //! Constructor
    Peridigm(const Teuchos::RCP<const Epetra_Comm>& comm,
             const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Create discretization object
    void createDiscretization();

    //! Destructor
    ~Peridigm(){};

  private:

    //! Parameterlist of entire input deck
    Teuchos::RCP<Teuchos::ParameterList> peridigmParams;

    //! Epetra communicator established by Peridigm_Factory
    Teuchos::RCP<const Epetra_Comm> peridigmComm;

    //! Output stream
    Teuchos::RCP<Teuchos::FancyOStream> out;

    //! Discretization object
    Teuchos::RCP<AbstractDiscretization> peridigmDisc;

  };
}

#endif // PERIDIGM_HPP
