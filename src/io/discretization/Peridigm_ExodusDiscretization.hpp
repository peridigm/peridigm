/*! \file Peridigm_ExodusDiscretization.hpp */

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

#ifndef PERIDIGM_EXODUSDISCRETIZATION_HPP
#define PERIDIGM_EXODUSDISCRETIZATION_HPP

//#define USE_STK

#ifdef USE_STK
// \todo Remove backwards compatibility after next Trilinos release
#include <Trilinos_version.h>
#if TRILINOS_MAJOR_MINOR_VERSION > 101002
#include <stk_io/MeshReadWriteUtils.hpp>
#else
#include <stk_io/util/UseCase_mesh.hpp>
#endif
#endif

#include "Peridigm_Discretization.hpp"
#include "Peridigm_InterfaceData.hpp"
#include <vector>
#include <map>

namespace PeridigmNS {

  //! Discretization class that reads an exodus/genesis file (Exodus II file format).
  class ExodusDiscretization : public PeridigmNS::Discretization {

  public:

    enum ExodusElementType { UNKNOWN_ELEMENT, SPHERE_ELEMENT, TET_ELEMENT, HEX_ELEMENT };

    //! Constructor
    ExodusDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetraComm,
                         const Teuchos::RCP<Teuchos::ParameterList>& params);

    //! Destructor
    virtual ~ExodusDiscretization();

    //! Return d-dimensional map
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOwnedMap(int d) const;

    //! Return d-dimensional overlap map (includes ghosts)
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOverlapMap(int d) const;

    //! Bond map, used for constitutive data stored on each bond. This is a non-overlapping map.
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalBondMap() const;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getInitialX() const;

    //! Get the horizon for each point
    virtual Teuchos::RCP<Epetra_Vector> getHorizon() const;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const;

    //! Get a vector containing the block ID of each element
    virtual Teuchos::RCP<Epetra_Vector> getBlockID() const;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::InterfaceData> getInterfaceData() const{
      return interfaceData;
    }

    // ! determine if the interface data is available
    virtual bool InterfacesAreConstructed() const{ return constructInterfaces;}

    //! Get the number of bonds on this processor
    virtual unsigned int getNumBonds() const;

    //! Get the number of elems on this processor
    virtual unsigned int getNumElem() const {return oneDimensionalMap->NumMyElements();}

    //! Get the maximum number of bonds per element on this processor
    virtual unsigned int getMaxNumBondsPerElem() const;

    //! Get the minimum element radius in the model (used for example for determining magnitude of finite-difference probe).
    virtual double getMinElementRadius() const { return minElementRadius; }

    //! Get the maximum element radius in the model
    virtual double getMaxElementRadius() const { return maxElementRadius; }

    //! Get the maximum element dimension (for example the diagonal of a hex element, used for partial volume neighbor search).
    virtual double getMaxElementDimension() const { return maxElementDimension; }

    //! Get the node positions in the original Exodus hex/tet mesh.
    virtual void getExodusMeshNodePositions(int globalNodeID, std::vector<double>& nodePositions);

  private:

    //! Compute the maximum element dimension.
    double computeMaxElementDimension();

    //! Private to prohibit copying
    ExodusDiscretization(const ExodusDiscretization&);

    //! Private to prohibit copying
    ExodusDiscretization& operator=(const ExodusDiscretization&);

    //! Loads mesh data into Epetra_Vectors (initial positions, volumes, block ids) and stores original Exodus node locations and connectivity.
    void loadData(const std::string& meshFileName);

  protected:

    template<class T>
    struct NonDeleter{
      void operator()(T* d) {}
    };

    //! Create vectors
    void createVectors();

    //! Create NeighborhoodData
    void createNeighborhoodData(int neighborListSize, int* neighborList);

    //! Create Interfaces between elements
    void constructInterfaceData();

    //! Filter bonds from neighborhood list
    Teuchos::RCP<PeridigmNS::NeighborhoodData> filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData);

    //! Refine the neighborhood list by eliminating neighbors that have zero intersection with the horizon (used only when computing element-sphere intersections).
    void removeNonintersectingNeighborsFromNeighborList(Teuchos::RCP<Epetra_Vector> x,
                                                        Teuchos::RCP<Epetra_Vector> searchRadii,
                                                        Teuchos::RCP<Epetra_BlockMap> ownedMap,
                                                        Teuchos::RCP<Epetra_BlockMap>& overlapMap,
                                                        int& neighborListSize,
                                                        int*& neighborList);

    //! Perform parallel communication to make exodus mesh data available for ghosted (overlap) element.
    void ghostExodusMeshData();

    //! Error reporting for calls to ExodusII API
    void reportExodusError(int errorCode, const char *methodName, const char *exodusMethodName);

    //! Verbosity flag
    bool verbose;

    //! Maps
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> bondMap;

    //! Minimum element radius
    double minElementRadius;

    //! Maximum element radius
    double maxElementRadius;

    //! Vector containing initial positions
    Teuchos::RCP<Epetra_Vector> initialX;

    //! Vector containing horizons
    Teuchos::RCP<Epetra_Vector> horizonForEachPoint;

    //! Vector containing cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolume;

    //! Vector containing the block ID of each element
    Teuchos::RCP<Epetra_Vector> blockID;

    //! Boolean flag for storing exodus mesh
    bool storeExodusMesh;

    //! Boolean flag for constructing interfaces
    bool constructInterfaces;

    //! Boolean flag indicating that element-horizon intersections should be computed
    bool computeIntersections;

    //! Maximum element dimension of the original exodus mesh
    double maxElementDimension;

    //! Vector containing node positions in the initial hex/tet mesh
    Teuchos::RCP<Epetra_Vector> exodusMeshNodePositions;

    //! Vector containing element connectivity in the initial hex/tet mesh
    Teuchos::RCP<Epetra_Vector> exodusMeshElementConnectivity;

    //! Struct containing neighborhoods for owned nodes.
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! Struct containing the interface varaibles and right and left elements
    Teuchos::RCP<PeridigmNS::InterfaceData> interfaceData;

    //! Returns number of bonds on this processor
    unsigned int numBonds;

    //! Returns the max number of bonds per element on this processor
    unsigned int maxNumBondsPerElem;

    //! Processor ID
    unsigned int myPID;

    //! Number of Processors
    unsigned int numPID;

    //! Discretization parameter controling the formation of bonds
    std::string bondFilterCommand;

#ifdef USE_STK
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
#endif

    //! Epetra communicator
    Teuchos::RCP<const Epetra_Comm> comm;
  };
}

#endif // PERIDIGM_EXODUSDISCRETIZATION_HPP
