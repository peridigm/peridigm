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

#include <Teuchos_RCP.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>

#include "mesh_input/quick_grid/QuickGrid.h"
#include "Peridigm_NeighborhoodData.hpp"

namespace PeridigmNS {

  class AbstractDiscretization {
  public:

    //! Constructor
    AbstractDiscretization() :
      horizon(0.0),
      sphereMeshElementBlocks(Teuchos::rcp(new std::map< std::string, std::vector<int> >())),
      sphereMeshNodeSets(Teuchos::rcp(new std::map< std::string, std::vector<int> >()))
    {}

    //! Destructor
    virtual ~AbstractDiscretization() {}

    //! Return the number of element blocks
    int getNumBlocks() { return (int)sphereMeshElementBlocks->size(); }

    //! Return a sorted list of element block names
    std::vector<std::string> getBlockNames() {
      std::vector<std::string> blockNames;
      std::map< std::string, std::vector<int> >::const_iterator it;
      for(it = sphereMeshElementBlocks->begin() ; it != sphereMeshElementBlocks->end() ; it++)
        blockNames.push_back(it->first);
      std::sort(blockNames.begin(), blockNames.end());
      return blockNames;
    }

    //! Return d-dimensional map
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalMap(int d) const = 0;

    //! Return d-dimensional overlap map
    // \todo Obsolete, remove.
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOverlapMap(int d) const = 0;

    /** \brief Bond map, used for constitutive data stored on each bond. This is
     *   a non-overlapping map. */
    // \todo Obsolete, remove.
    virtual Teuchos::RCP<const Epetra_BlockMap> getBondMap() const = 0;

    //! Return d-dimensional map for the given element block
    virtual Teuchos::RCP<const Epetra_BlockMap> getElementBlockOwnedMap(std::string& blockName, int dimension) const = 0;

    //! Return d-dimensional overlap map for the given element block (owned points + ghosts)
    virtual Teuchos::RCP<const Epetra_BlockMap> getElementBlockOverlapMap(std::string& blockName, int dimension) const = 0;

    //! Return 1-dimensional bond map for the given element block
    virtual Teuchos::RCP<const Epetra_BlockMap> getElementBlockBondMap(std::string& blockName) const = 0;

    //! Get the neighbor list for all locally-owned nodes for the given element block
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getElementBlockNeighborhoodData(std::string& blockName) const = 0;

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

    //! Get the locally-owned IDs for each element block
    virtual Teuchos::RCP< std::map< std::string, std::vector<int> > > getElementBlocks() { return sphereMeshElementBlocks; } ;

    //! Get the locally-owned IDs for each node set
    virtual Teuchos::RCP< std::map< std::string, std::vector<int> > > getNodeSets() { return sphereMeshNodeSets; } ;

    //! Get the owned (non-overlap) map.
    static Epetra_BlockMap getOwnedMap(const Epetra_Comm& comm, const QUICKGRID::Data& gridData, int ndf);

    //! Get the overlap map.
    static Epetra_BlockMap getOverlapMap(const Epetra_Comm& comm, const QUICKGRID::Data& gridData, int ndf);

  protected:

    //! Get the overlap map.
    static Epetra_BlockMap getOverlap(int ndf, int numShared, int*shared, int numOwned, const  int* owned, const Epetra_Comm& comm);

    //! Get the shared global IDs.
    static UTILITIES::Array<int> getSharedGlobalIds(const QUICKGRID::Data& gridData);

    //! Get the local owned IDs.
    static shared_ptr<int> getLocalOwnedIds(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap);

    //! Get the local neighborhood list.
    static shared_ptr<int> getLocalNeighborList(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap);

    //! Horizon
    double horizon;

    //! Map containing element blocks (block name and list of locally-owned element IDs for each block).
    Teuchos::RCP< std::map< std::string, std::vector<int> > > sphereMeshElementBlocks;

    //! Map containing node sets (node set name and list of locally-owned node IDs for each node set).
    Teuchos::RCP< std::map< std::string, std::vector<int> > > sphereMeshNodeSets;

  private:

    //! Private to prohibit copying.
    AbstractDiscretization(const AbstractDiscretization&);

    //! Private to prohibit copying.
    AbstractDiscretization& operator=(const AbstractDiscretization&);
  };

}

#endif // PERIDIGM_ABSTRACTDISCRETIZATION_HPP
