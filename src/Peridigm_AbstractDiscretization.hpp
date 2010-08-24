/*! \file Peridigm_AbstractDiscretization.hpp */

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

#ifndef PERIDIGM_ABSTRACTDISCRETIZATION_HPP
#define PERIDIGM_ABSTRACTDISCRETIZATION_HPP

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "Peridigm_NeighborhoodData.hpp"

namespace Peridigm {

  class AbstractDiscretization {
  public:

    //! Constructor.
    AbstractDiscretization() : horizon(0.0) {}

    //! Destructor.
    virtual ~AbstractDiscretization() {}

    //! One-dimensional map, used for cell volumes and scalar constitutive data.
    virtual Teuchos::RCP<const Epetra_BlockMap> getOneDimensionalMap() const = 0;

    //! One-dimensional overlap map, used for cell volumes and scalar constitutive data (includes ghosts).
    virtual Teuchos::RCP<const Epetra_BlockMap> getOneDimensionalOverlapMap() const = 0;

    /** \brief Bond map, used for constitutive data stored on each bond. This is
     *   a non-overlapping map. */
    virtual Teuchos::RCP<const Epetra_BlockMap> getBondMap() const = 0;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getSolverInitialX() const = 0;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const = 0;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<Peridigm::NeighborhoodData> getNeighborhoodData() const = 0;

    //! Get the number of bonds on this processor
    virtual unsigned int getNumBonds() const = 0;

    //! Get the horizon
    double getHorizon() const { return horizon; }

  protected:

    //! Family horizon
    double horizon;

  private:

    //! Private to prohibit copying.
    AbstractDiscretization(const AbstractDiscretization&);

    //! Private to prohibit copying.
    AbstractDiscretization& operator=(const AbstractDiscretization&);
  };

}

#endif // PERIDIGM_ABSTRACTDISCRETIZATION_HPP
