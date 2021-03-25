/*! \file Peridigm_AlbanyDiscretization.hpp */

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

#ifndef PERIDIGM_ALBANYDISCRETIZATION_HPP
#define PERIDIGM_ALBANYDISCRETIZATION_HPP

#include "Peridigm_Discretization.hpp"
#include <Epetra_Vector.h>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Comm.h>
#include "mpi.h"

namespace PeridigmNS {

  //! Discretization class that creates discretization from a Albany data.
  class AlbanyDiscretization : public PeridigmNS::Discretization {

  public:

    //! Constructor
    AlbanyDiscretization(const MPI_Comm& mpiComm,
                         const Teuchos::RCP<Teuchos::ParameterList>& params,
                         int numGlobalIds,
                         const int* globalIds,
                         const double* refCoord,
                         const double* volume,
                         const int* blockId,
                         const int numNodeIds = 0,
                         int * nodeGlobalIds = 0,
                         const double* nodeCoord = 0,
                         const int* nodeblockId = 0);

    //! Destructor
    virtual ~AlbanyDiscretization();

    //! Return d-dimensional map
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOwnedMap(int d) const;

    //! Return d-dimensional overlap map (includes ghosts)
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalOverlapMap(int d) const;

    //! Bond map, used for constitutive data stored on each bond. This is a non-overlapping map.
    virtual Teuchos::RCP<const Epetra_BlockMap> getGlobalBondMap() const;

    //! Get initial positions
    virtual Teuchos::RCP<Epetra_Vector> getInitialX() const;

    //! Get the horizon value for each point.
    virtual Teuchos::RCP<Epetra_Vector> getHorizon() const;

    //! Get cell volumes
    virtual Teuchos::RCP<Epetra_Vector> getCellVolume() const;

    //! Get a vector containing the block ID of each element
    virtual Teuchos::RCP<Epetra_Vector> getBlockID() const;

    //! Get the neighbor list for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::NeighborhoodData> getNeighborhoodData() const;

    //! Get the neighbor list for all locally-owned nodes
    Teuchos::RCP<PeridigmNS::NeighborhoodData> getAlbanyPartialStressNeighborhoodData() const{return albanyPartialStressNeighborhoodData;}

    //! Get interface data for all locally-owned nodes
    virtual Teuchos::RCP<PeridigmNS::InterfaceData> getInterfaceData() const{return Teuchos::null;}

    //! determine if interface data was constructed
    virtual bool InterfacesAreConstructed() const{return false;}

    //! Get the number of bonds on this processor
    virtual unsigned int getNumBonds() const;

    //! Get the number of elems on this processor
    virtual unsigned int getNumElem() const {return oneDimensionalMap->NumMyElements();}

    //! Get the maximum number of bonds per element on this processor
    virtual unsigned int getMaxNumBondsPerElem() const;

    //! Get the minimum element radius in the model (used for example for determining magnitude of finite-difference probe).
    virtual double getMinElementRadius() const { return minElementRadius; }

    //! Get the maximum element radius in the model (used for example for determining magnitude of finite-difference probe).
    virtual double getMaxElementRadius() const { return maxElementRadius; }

    //! Get the maximum element dimension (for example the diagonal of a hex element, used for partial volume neighbor search).
    virtual double getMaxElementDimension() const { return maxElementDimension; }

    //! Return the albany interface 1D mapxs
    Teuchos::RCP<Epetra_BlockMap> getAlbanyInterface1DMap()const {return albanyInterface1DMap;}

  private:

    //! Private to prohibit copying
    AlbanyDiscretization(const AlbanyDiscretization&);

    //! Private to prohibit copying
    AlbanyDiscretization& operator=(const AlbanyDiscretization&);

    //! Create a list of on-processor elements for each block
    void createBlockElementLists();

  protected:

    //! Create NeighborhoodData
    void createNeighborhoodData(int neighborListSize, int* neighborList);

    //! Filter bonds from neighborhood list
    Teuchos::RCP<PeridigmNS::NeighborhoodData> filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData);

    //! Maps
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> albanyInterface1DMap;
    Teuchos::RCP<Epetra_BlockMap> oneDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalMap;
    Teuchos::RCP<Epetra_BlockMap> threeDimensionalOverlapMap;
    Teuchos::RCP<Epetra_BlockMap> bondMap;

    //! Minimum element radius
    double minElementRadius;

    //! Maximum element radius
    double maxElementRadius;

    //! Maximum element dimension
    double maxElementDimension;

    //! Vector containing initial positions
    Teuchos::RCP<Epetra_Vector> initialX;

    //! Vector containing horizons
    Teuchos::RCP<Epetra_Vector> horizonForEachPoint;

    //! Vector containing cell volumes
    Teuchos::RCP<Epetra_Vector> cellVolume;

    //! Vector containing the block ID of each element
    Teuchos::RCP<Epetra_Vector> blockID;

    //! Struct containing neighborhoods for owned nodes.
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;

    //! Struct containing neighborhoods for albany partial stress nodes
    Teuchos::RCP<PeridigmNS::NeighborhoodData> albanyPartialStressNeighborhoodData;

    //! Returns number of bonds on this processor
    unsigned int numBonds;

    //! Returns the max number of bonds per element on this processor
    unsigned int maxNumBondsPerElem;

    //! Discretization parameter controling the formation of bonds
    std::string bondFilterCommand;

    //! Epetra communicator
    Teuchos::RCP<const Epetra_Comm> comm;
  };
}

#endif // PERIDIGM_ALBANYDISCRETIZATION_HPP
