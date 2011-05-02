/*! \file Peridigm_STKDiscretization.hpp */

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

#ifndef PERIDIGM_STKDISCRETIZATION_HPP
#define PERIDIGM_STKDISCRETIZATION_HPP

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Comm.h>
#include "Peridigm_AbstractDiscretization.hpp"
#include "PdGridData.h"

namespace PeridigmNS {

  //! Discretization class that creates discretizations using STK (reads Exodus II mesh files).
  class STKDiscretization : public PeridigmNS::AbstractDiscretization {

  public:

    //! Constructor
    STKDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetraComm,
                      const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~STKDiscretization();

    //! Return d-dimensional map
    virtual Teuchos::RCP<const Epetra_BlockMap> getMap(int d) const;

    //! Return d-dimensional overlap map (includes ghosts)
    virtual Teuchos::RCP<const Epetra_BlockMap> getOverlapMap(int d) const;

    //! Bond map, used for constitutive data stored on each bond. This is a non-overlapping map.
    virtual Teuchos::RCP<const Epetra_BlockMap> getBondMap() const;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getInitialX() const;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const;

    //! Get the number of bonds on this processor
    unsigned int getNumBonds() const;

  private:

    //! Private to prohibit copying
    STKDiscretization(const STKDiscretization&);

    //! Private to prohibit copying
    STKDiscretization& operator=(const STKDiscretization&);

  protected:

    //! Create maps
    void createMaps(const PdGridData& decomp);

    //! Create vectors
    void createVectors();

    //! Create NeighborhoodData
    void createNeighborhoodData(const PdGridData& decomp);

    //! Epetra communicator
    Teuchos::RCP<const Epetra_Comm> comm;

    //! Maps
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> bondMap;

    //! Vector containing initial positions
    Teuchos::RCP<Epetra_Vector> initialX;

    //! Vector containing cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolume;
	
    //! Struct containing neighborhoods for owned nodes.
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! Returns number of bonds on this processor
    unsigned int numBonds;

    //! Processor ID
    unsigned int myPID;

    //! Number of Processors
    unsigned int numPID;
  };
}

#endif // PERIDIGM_STKDISCRETIZATION_HPP
