/*! \file Peridigm_AbstractDiscretization.hpp */

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

#ifndef PERIDIGM_ABSTRACTDISCRETIZATION_HPP
#define PERIDIGM_ABSTRACTDISCRETIZATION_HPP

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>

#include "Peridigm_NeighborhoodData.hpp"

namespace PeridigmNS {

  class AbstractDiscretization {
  public:

    //! Constructor.
    AbstractDiscretization() : horizon(0.0) {}

    //! Destructor.
    virtual ~AbstractDiscretization() {}

    //! Return d-dimensional map
    virtual Teuchos::RCP<const Epetra_BlockMap> getMap(int d) const = 0;

    //! Return d-dimensional overlap map
    virtual Teuchos::RCP<const Epetra_BlockMap> getOverlapMap(int d) const = 0;

    /** \brief Bond map, used for constitutive data stored on each bond. This is
     *   a non-overlapping map. */
    virtual Teuchos::RCP<const Epetra_BlockMap> getBondMap() const = 0;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getInitialX() const = 0;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const = 0;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const = 0;

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
