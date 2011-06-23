/*
 * Neighborhood.h
 *
 *  Created on: Feb 25, 2011
 *      Author: jamitch
 */

#ifndef NEIGHBORHOOD_H_
#define NEIGHBORHOOD_H_

#include <tr1/memory>
#include <set>

#include "BondFilter.h"
#include "utilities/Array.h"

#include "Epetra_BlockMap.h"
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
* Produces Epetra overlap maps (scalar and vector) for a given horizon and set of points.
* Can produce Epetra_MultiVectors (owned and overlap).
* Can perform communications on Epetra_MultiVectors,
 *  for example -- spreading locally owned values to overlap vectors.
*/

/*
 * Use this 'class' to perform parallel search and create a new neighborhood list;
 * This function will construct neighborhood lists including across processor boundaries;
 * @param horizon -- this is the distance that should be used to form the neighborhood list
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
using std::tr1::shared_ptr;
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
			shared_ptr<Epetra_Comm> comm,
			struct Zoltan_Struct* zz,
			size_t numOwnedPoints,
			shared_ptr<int>& ownedGIDs,
			shared_ptr<double>& owned_coordinates,
			double horizon,
			shared_ptr<PdBondFilter::BondFilter> bondFilterPtr = shared_ptr<PdBondFilter::BondFilter>(new PdBondFilter::BondFilterDefault())
			);
	double get_horizon() const;
	size_t get_num_owned_points() const;
	size_t get_num_shared_points() const;
	int get_num_neigh (int localId) const;
	shared_ptr<int> get_neighborhood_ptr() const;
	shared_ptr<int> get_neighborhood() const;
	shared_ptr<int> get_local_neighborhood() const;
	shared_ptr<int> get_shared_gids() const;
	const int* get_neighborhood (int localId) const;
	const int* get_local_neighborhood (int localId) const;

	int get_size_neighborhood_list() const;
	shared_ptr<double> get_owned_x() const;
	shared_ptr<Epetra_BlockMap> getOwnedMap(int ndf) const;
	shared_ptr<Epetra_BlockMap> getOverlapMap(int ndf) const;
	shared_ptr<Epetra_Comm> get_Epetra_Comm() const;
	shared_ptr<Epetra_Distributor> create_Epetra_Distributor() const;
	/*
	 * This function is primarily intended for internal force operators
	 * that do not include 'x' in the neighborhood H(x) but it is desirable
	 * to have 'x' included in the newly cloned neighborhood -- hence the
	 * default value 'withSelf=true'
	 */
	NeighborhoodList cloneAndShare(double newHorizon, bool withSelf=true);

private:

	void buildNeighborhoodList(int numOverlapPoints,shared_ptr<double> xOverlapPtr);
	Array<int> createLocalNeighborList(const Epetra_BlockMap& overlapMap);
	Array<int> createSharedGlobalIds() const;
	void createAndAddNeighborhood();
	shared_ptr<Epetra_BlockMap> create_Epetra_BlockMap(Epetra_MapTag key);

private:
	shared_ptr<Epetra_Comm> epetraComm;
	std::map<Epetra_MapTag, shared_ptr<Epetra_BlockMap>, MapComparator > epetra_block_maps;
	size_t num_owned_points, size_neighborhood_list;
	double horizon;
	shared_ptr<int> owned_gids;
	shared_ptr<double> owned_x;
	Array<int> neighborhood, local_neighborhood, neighborhood_ptr, num_neighbors, sharedGIDs;
	struct Zoltan_Struct* zoltan;
	shared_ptr<PdBondFilter::BondFilter> filter_ptr;

};


}

#endif /* NEIGHBORHOOD_H_ */
