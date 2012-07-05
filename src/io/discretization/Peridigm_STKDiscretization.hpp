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

// \todo Remove this include after next Trilinos release
#include <Trilinos_version.h>

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Comm.h>

// \todo Remove backwards compatibility after next Trilinos release
#if TRILINOS_MAJOR_MINOR_VERSION > 101002
#include <stk_io/MeshReadWriteUtils.hpp>
#else
#include <stk_io/util/UseCase_mesh.hpp>
#endif

#include "Peridigm_AbstractDiscretization.hpp"
#include "mesh_input/quick_grid/QuickGridData.h"

#include <vector>
#include <map>

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
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOwnedMap(int d) const;

    //! Return d-dimensional overlap map (includes ghosts)
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOverlapMap(int d) const;

    //! Bond map, used for constitutive data stored on each bond. This is a non-overlapping map.
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalBondMap() const;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getInitialX() const;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const;

    //! Get a vector containing the block ID of each element
    virtual Teuchos::RCP<Epetra_Vector> getBlockID() const;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const;

    //! Get the number of bonds on this processor
    virtual unsigned int getNumBonds() const;

    //! Get the horizon
    virtual double getHorizon() const { return horizon; }

    //! Get the minimum element radius in the model (used for example for determining magnitude of finite-difference probe).
    virtual double getMinElementRadius() const { return minElementRadius; }

    //! Get the maximum element dimension (for example the diagonal of a hex element, used for partial volume neighbor search).
    virtual double getMaxElementDimension() const { return maxElementDimension; }

    //! Get the node positions in the original Exodus hex/tet mesh.
    virtual Teuchos::RCP< std::vector<double> > getExodusMeshNodePositions(int globalNodeID);

  private:

    //! Private to prohibit copying
    STKDiscretization(const STKDiscretization&);

    //! Private to prohibit copying
    STKDiscretization& operator=(const STKDiscretization&);

    //! Creates a discretization object using STK functionality.
    QUICKGRID::Data getDecomp(const std::string& meshFileName,
                              double horizon);

  protected:

    template<class T>
    struct NonDeleter{
      void operator()(T* d) {}
    };

    //! Create maps
    void createMaps(const QUICKGRID::Data& decomp);

    //! Create vectors
    void createVectors();

    //! Create NeighborhoodData
    void createNeighborhoodData(const QUICKGRID::Data& decomp);

    //! Compute the scalar triple product
    double scalarTripleProduct(std::vector<double>& a,
                               std::vector<double>& b,
                               std::vector<double>& c) const;

    //! Compute the volume of a hexahedron element
    double hexVolume(std::vector<double*>& nodeCoordinates) const;

    //! Compute the maximum length dimension for a hexahedron element
    double hexMaxElementDimension(std::vector<double*>& nodeCoordinates) const;

    //! Maps
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> bondMap;

    //! Horizon
    double horizon;

    //! Minimum element radius
    double minElementRadius;

    //! Maximum element dimension
    double maxElementDimension;

    //! Search horizon, which may be larger than the horizon if partial volumes are used
    double searchHorizon;

    //! Vector containing initial positions
    Teuchos::RCP<Epetra_Vector> initialX;

    //! Vector containing cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolume;

    //! Vector containing the block ID of each element
    Teuchos::RCP<Epetra_Vector> blockID;

    //! Vector containing node positions in the initial hex/tet mesh
    Teuchos::RCP<Epetra_Vector> exodusMeshNodePositions;

    //! Vector containing element connectivity in the initial hex/tet mesh
    std::map< int, std::vector<int> > exodusMeshElementConnectivity;

    //! Map between block name (from Exodus file) and block number (Peridigm numbering)
    std::map<std::string, int> blockNameToBlockNumberMap;

    //! Struct containing neighborhoods for owned nodes.
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! Returns number of bonds on this processor
    unsigned int numBonds;

    //! Processor ID
    unsigned int myPID;

    //! Number of Processors
    unsigned int numPID;

    //! Mesh meta data
    Teuchos::RCP<stk::mesh::fem::FEMMetaData> metaData;

    // \todo Remove backwards compatibility after next Trilinos release
    #if TRILINOS_MAJOR_MINOR_VERSION > 101002
    //! Mesh bulk data
    Teuchos::RCP<stk::io::MeshData> meshData;
    #else
    //! Mesh bulk data
    Teuchos::RCP<stk::io::util::MeshData> meshData;
    #endif

    //! Epetra communicator
    Teuchos::RCP<const Epetra_Comm> comm;
  };
}

#endif // PERIDIGM_STKDISCRETIZATION_HPP
