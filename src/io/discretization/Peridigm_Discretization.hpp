/*! \file Peridigm_Discretization.hpp */

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

#ifndef PERIDIGM_DISCRETIZATION_HPP
#define PERIDIGM_DISCRETIZATION_HPP

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include "Peridigm_NeighborhoodData.hpp"
#include "Peridigm_InterfaceData.hpp"
#include "QuickGrid.h"
#include "BondFilter.h"

namespace PeridigmNS {

  class Discretization {
  public:

    //! Constructor
    Discretization() :
      elementBlocks(Teuchos::rcp(new std::map< std::string, std::vector<int> >())),
      nodeSets(Teuchos::rcp(new std::map< std::string, std::vector<int> >()))
    {}

    //! Destructor
    virtual ~Discretization() {}

    //! Return the number of element blocks
    int getNumBlocks() { return (int)elementBlocks->size(); }

    //! Return a sorted list of element block names
    std::vector<std::string> getBlockNames() {
      std::vector<std::string> blockNames;
      std::map< std::string, std::vector<int> >::const_iterator it;
      for(it = elementBlocks->begin() ; it != elementBlocks->end() ; it++)
        blockNames.push_back(it->first);
      std::sort(blockNames.begin(), blockNames.end());
      return blockNames;
    }

    //! Return d-dimensional map
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOwnedMap(int d) const = 0;

    //! Return d-dimensional overlap map
    // \todo Obsolete, remove.
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOverlapMap(int d) const = 0;

    /** \brief Bond map, used for constitutive data stored on each bond. This is
     *   a non-overlapping map. */
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalBondMap() const = 0;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getInitialX() const = 0;

    //! Get the horizon for each node
    virtual Teuchos::RCP<Epetra_Vector> getHorizon() const = 0;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const = 0;

    //! Get a vector containing the block ID of each element
    virtual Teuchos::RCP<Epetra_Vector> getBlockID() const = 0;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const = 0;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::InterfaceData> getInterfaceData() const = 0;

    //! determine if interfaces are constructed
    virtual bool InterfacesAreConstructed() const = 0;

    //! Get the number of bonds on this processor
    virtual unsigned int getNumBonds() const = 0;

    //! Get the number of elems on this processor
    virtual unsigned int getNumElem() const = 0;

    //! Get the max number of bonds per cell
    virtual unsigned int getMaxNumBondsPerElem() const = 0;

    //! Get the minimum element radius in the model (used for example for determining magnitude of finite-difference probe).
    virtual double getMinElementRadius() const = 0;

    //! Get the minimum element radius in the model (used for example for determining magnitude of finite-difference probe).
    virtual double getMaxElementRadius() const = 0;

    //! Get the maximum element dimension (for example the diagonal of a hex element, used for partial volume neighbor search).
    virtual double getMaxElementDimension() const = 0;

    //! Get the locally-owned IDs for each element block
    virtual Teuchos::RCP< std::map< std::string, std::vector<int> > > getElementBlocks() { return elementBlocks; } ;

    //! Get the locally-owned IDs for each node set
    virtual Teuchos::RCP< std::map< std::string, std::vector<int> > > getNodeSets() { return nodeSets; } ;

    //! Get the locally-owned IDs for each node set
    Teuchos::RCP< std::map< std::string, int> > getNodeSetIds() { return nodeSetIds; } ;

    //! Get the node positions in the original Exodus hex/tet mesh.
    virtual void getExodusMeshNodePositions(int globalNodeID, std::vector<double>& nodePositions){
      // The default implementation sets the nodePositions vector to length zero.
      nodePositions.clear();
      return;
    }

    //! Get the owned (non-overlap) map.
    static Epetra_BlockMap getOwnedMap(const Epetra_Comm& comm, const QUICKGRID::Data& gridData, int ndf);

    //! Get the overlap map.
    static Epetra_BlockMap getOverlapMap(const Epetra_Comm& comm, const QUICKGRID::Data& gridData, int ndf);

    void createBondFilters(const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Get the block id for a given block name
    int blockNameToBlockId(std::string blockName) const;

  protected:

    //! Get the overlap map.
    static Epetra_BlockMap getOverlap(int ndf, int numShared, int*shared, int numOwned, const  int* owned, const Epetra_Comm& comm);

    //! Get the shared global IDs.
    static UTILITIES::Array<int> getSharedGlobalIds(const QUICKGRID::Data& gridData);

    //! Get the local owned IDs.
    static std::tr1::shared_ptr<int> getLocalOwnedIds(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap);

    //! Get the local neighborhood list.
    static std::tr1::shared_ptr<int> getLocalNeighborList(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap);

    //! \todo Eliminate old-style elementBlocks data structure.
    //! Map containing element blocks (block name and list of locally-owned element IDs for each block).
    Teuchos::RCP< std::map< std::string, std::vector<int> > > elementBlocks;

    //! Map containing node sets (node set name and list of locally-owned node IDs for each node set).
    Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets;

    //! Map containing the node id for each node set
    Teuchos::RCP< std::map< std::string, int> > nodeSetIds;

    std::vector< std::tr1::shared_ptr<PdBondFilter::BondFilter> > bondFilters;

  private:

    //! Private to prohibit copying.
    Discretization(const Discretization&);

    //! Private to prohibit copying.
    Discretization& operator=(const Discretization&);
  };

}

#endif // PERIDIGM_DISCRETIZATION_HPP
