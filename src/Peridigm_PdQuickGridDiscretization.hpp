/*! \file Peridigm_PdQuickGridDiscretization.hpp */

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

#ifndef PERIDIGM_TESTDISCRETIZATION_HPP
#define PERIDIGM_TESTDISCRETIZATION_HPP

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Comm.h>
#include "Peridigm_AbstractDiscretization.hpp"
#include "PdGridData.h"

namespace Peridigm {

  //! Discretization class that creates discretizations using PdQuickGrid.
  class PdQuickGridDiscretization : public Peridigm::AbstractDiscretization {
  public:

    //! Constructor
    PdQuickGridDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetraComm,
							  const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~PdQuickGridDiscretization();

    //! Returns the discretization object, switches on types of PdQuickGrids.
    PdGridData getDiscretization(const Teuchos::RCP<Teuchos::ParameterList>& param);

    //! One-dimensional map, used for cell volumes and scalar constitutive data.
    virtual Teuchos::RCP<const Epetra_BlockMap> getOneDimensionalMap() const;

    //! One-dimensional overlap map, used for cell volumes and scalar constitutive data (includes ghosts).
    virtual Teuchos::RCP<const Epetra_BlockMap>  getOneDimensionalOverlapMap() const;

	/** \brief Bond map, used for constitutive data stored on each bond. This is
     *   a non-overlapping map. */
    virtual Teuchos::RCP<const Epetra_BlockMap> 
    getBondMap() const;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector>
    getSolverInitialX() const;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector>
    getCellVolume() const;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<Peridigm::NeighborhoodData> getNeighborhoodData() const;

	//! Get the number of bonds on this processor
	unsigned int getNumBonds() const;

  private:

    //! Private to prohibit copying
    PdQuickGridDiscretization(const PdQuickGridDiscretization&);

    //! Private to prohibit copying
    PdQuickGridDiscretization& operator=(const PdQuickGridDiscretization&);

  protected:

    //! Create maps
    void createMaps(PdGridData& decomp);

    //! Create vectors
    void createVectors();

	//! Create NeighborhoodData
	void createNeighborhoodData(PdGridData& decomp);

    //! Epetra communicator
    Teuchos::RCP<const Epetra_Comm> comm;

    //! Maps
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> bondMap;

    //! Vector containing initial positions and velocities
    Teuchos::RCP<Epetra_Vector> solverInitialX;

    //! Vector containing cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolume;
	
	//! Struct containing neighborhoods for owned nodes.
    Teuchos::RCP<Peridigm::NeighborhoodData> neighborhoodData;

	//! Returns number of bonds on this processor
	unsigned int numBonds;

    //! Processor ID
    unsigned int myPID;

    //! Number of Processors
    unsigned int numPID;
  };
}

#endif // PERIDIGM_TESTDISCRETIZATION_HPP
