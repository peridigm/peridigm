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

#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

#include <set>

#include "BondFilter.h"
#include "Array.h"

#include <Teuchos_RCP.hpp>
//#include <Epetra_BlockMap.h>
#include <Epetra_Vector.h>
#include <vector>
#include <map>


class Epetra_Comm;
struct Zoltan_Struct;
class Epetra_Distributor;

/**
 *
What is a neighborhood and what does it do?
* Is associated with a horizon.
* Contains a list of Global Ids (this is similar or perhaps the same as an Epetra_Map).
* Stores owned coordinates (these coordinates are the basis for the neighborhood calculation for each point).
* Is associated with a set of filters -- filters prevent bonds from crossing particular surfaces.
* Produces Epetra overlap maps (scalar and vector) for a set of points and a corresponding set of horizons.
* Can produce Epetra_MultiVectors (owned and overlap).
* Can perform communications on Epetra_MultiVectors,
 *  for example -- spreading locally owned values to overlap vectors.
*/

/*
 * Use this 'class' to perform parallel search and create a new neighborhood list;
 * This function will construct neighborhood lists including across processor boundaries;
 * @param horizon -- this is the distance that should be used to form the neighborhood list for a given point
 * Use CASE Scenario:
 * 1) Create a mesh
 * 2) Load balance mesh
 * 3) Create NeighborhoodList
 * 4) 4th argument -- BondFilter which defaults to neighborhood search that does not include 'self' point
*/

namespace PDNEIGH {

/**
 * Utilities
 */
using std::shared_ptr;
using std::size_t;
using UTILITIES::Array;

template<class T> struct ArrayDeleter{
  void operator()(T* d) {
    delete [] d;
  }
};

class NeighborhoodList {

class Epetra_MapTag;
friend class Epetra_MapTag;
private:
  enum Epetra_MapType { OWNED=0, OVERLAP=1 };
  class Epetra_MapTag {
  public:
    explicit Epetra_MapTag(Epetra_MapType t, size_t block_size) : type(t), ndf(block_size) {}
    const Epetra_MapType type;
    const size_t ndf;
  };

  struct MapComparator {
    bool operator() (Epetra_MapTag left, Epetra_MapTag right) const {
      if(left.type < right.type)
        return true;
      return left.ndf < right.ndf;
    }
  };

public:
  NeighborhoodList(
      shared_ptr<const Epetra_Comm> comm,
      struct Zoltan_Struct* zz,
      size_t numOwnedPoints,
      shared_ptr<int> ownedGIDs,
      shared_ptr<double> owned_coordinates,
      Teuchos::RCP<Epetra_Vector> horizonList,
      std::vector< shared_ptr<PdBondFilter::BondFilter> > bondFilters = std::vector< shared_ptr<PdBondFilter::BondFilter> >()
      );
  NeighborhoodList(
      shared_ptr<const Epetra_Comm> comm,
      struct Zoltan_Struct* zz,
      size_t numOwnedPoints,
      shared_ptr<int> ownedGIDs,
      shared_ptr<double> owned_coordinates,
      double horizon,
      std::vector< shared_ptr<PdBondFilter::BondFilter> > bondFilters = std::vector< shared_ptr<PdBondFilter::BondFilter> >()
      );
  double get_frameset_buffer_size() const;
  size_t get_num_owned_points() const;
  size_t get_num_shared_points() const;
  int get_num_neigh (int localId) const;
  shared_ptr<int> get_neighborhood_ptr() const;
  shared_ptr<int> get_neighborhood() const;
  shared_ptr<int> get_local_neighborhood() const;
  shared_ptr<int> get_owned_gids() const;
  shared_ptr<int> get_shared_gids() const;
  const int* get_neighborhood (int localId) const;
  const int* get_local_neighborhood (int localId) const;

  int get_size_neighborhood_list() const;
  shared_ptr<double> get_owned_x() const;
  shared_ptr<Epetra_BlockMap> getOwnedMap(int ndf) const;
  shared_ptr<Epetra_BlockMap> getOverlapMap(int ndf) const;
  shared_ptr<const Epetra_Comm> get_Epetra_Comm() const;
  shared_ptr<Epetra_Distributor> create_Epetra_Distributor() const;

private:

  void buildNeighborhoodList(int numOverlapPoints,shared_ptr<double> xOverlapPtr);
  Array<int> createLocalNeighborList(const Epetra_BlockMap& overlapMap);
  Array<int> createSharedGlobalIds() const;
  void createAndAddNeighborhood();
  shared_ptr<Epetra_BlockMap> create_Epetra_BlockMap(Epetra_MapTag key);

private:
  shared_ptr<const Epetra_Comm> epetraComm;
  std::map<Epetra_MapTag, shared_ptr<Epetra_BlockMap>, MapComparator > epetra_block_maps;
  size_t num_owned_points, size_neighborhood_list;
  double frameset_buffer_size;
    Teuchos::RCP<Epetra_Vector> horizons;
  shared_ptr<int> owned_gids;
  shared_ptr<double> owned_x;
  Array<int> neighborhood, local_neighborhood, neighborhood_ptr, num_neighbors, sharedGIDs;
  struct Zoltan_Struct* zoltan;
  std::vector< shared_ptr<PdBondFilter::BondFilter> > filter_ptrs;

};


}

#endif /* NEIGHBORHOOD_H_ */
