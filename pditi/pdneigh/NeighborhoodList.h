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
#include "Epetra_BlockMap.h"
#include "PdBondFilter.h"

class Epetra_Comm;
struct Zoltan_Struct;
/**
 *
 *
What does a neighborhood do?
* Is associated with a horizon
* Contains a list of Global Ids (this is similar or perhaps the same as an Epetra_Map)
* Stores owned coordinates (these coordinates are the basis for the neighborhood calculation for each point)
* Is associated with a set of filters -- filters prevent bonds from crossing particular surfaces
* Produces Epetra overlap maps (scalar and vector)for a given horizon and set of points
* Can produce Epetra_MultiVectors (owned and overlap)
* Can perform communications on Epetra_MultiVectors
 * for example -- spreading locally owned values to overlap vectors

/*
 * Use this function to perform parallel search and create a new neighborhood list;
 * This function will construct neighborhood lists including across processor boundaries;
 * The resulting list will then be added to the incoming pdGridData -- note that any pre-existing
 * list on pdGridData will be overwritten/lost.
 * @param horizon -- this is the distance that should be used to form the neighborhood list
 * Use CASE Scenario:
 * 1) Create a mesh
 * 2) Load balance mesh
 * 3) Call this function
 * 4) 4th argument -- BondFilter which defaults to neighborhood search that does not include 'self' point
*/

namespace PDNEIGH {

/**
 * Utilities
 */
using std::tr1::shared_ptr;
using std::size_t;
using std::pair;
template<class T> struct ArrayDeleter{
	void operator()(T* d) {
		delete [] d;
	}
};


class NeighborhoodList {

public:
	NeighborhoodList(
			const Epetra_Comm& comm,
			struct Zoltan_Struct* zz,
			size_t numOwnedPoints,
			shared_ptr<int>& ownedGIDs,
			shared_ptr<double>& owned_coordinates,
			double horizon,
			shared_ptr<PdBondFilter::BondFilter> bondFilterPtr = shared_ptr<PdBondFilter::BondFilter>(new PdBondFilter::BondFilterDefault())
			);
	double get_horizon() const;
	int get_num_owned_points() const;
	int get_num_neigh (int localId) const;
	shared_ptr<int> get_neighborhood_ptr() const;
	shared_ptr<int> get_neighborhood() const;
	shared_ptr<int> get_local_neighborhood() const;
	const int* get_neighborhood (int localId) const;
	const int* get_local_neighborhood (int localId) const;
	int get_size_neighborhood_list() const;
	shared_ptr<double> get_owned_x() const;
	const Epetra_BlockMap getOverlapMap(const Epetra_Comm& comm, int ndf) const;
	const Epetra_BlockMap getOwnedMap(const Epetra_Comm& comm, int ndf) const;
	NeighborhoodList cloneAndShare(double newHorizon);
	shared_ptr< std::set<int> > constructParallelDecompositionFrameSet() const;

private:
	pair<int, shared_ptr<int> > getSharedGlobalIds() const;
	void buildNeighborhoodList(int numOverlapPoints,shared_ptr<double> xOverlapPtr);
	shared_ptr<int> createLocalNeighborList(const Epetra_BlockMap& overlapMap);
	void createAndAddNeighborhood();

private:
	const Epetra_Comm& epetraComm;
	size_t num_owned_points, size_neighborhood_list;
	double horizon;
	shared_ptr<int> owned_gids;
	shared_ptr<double> owned_x;
	shared_ptr<int> neighborhood, local_neighborhood, neighborhood_ptr, num_neighbors;
	struct Zoltan_Struct* zoltan;
	shared_ptr<PdBondFilter::BondFilter> filter_ptr;

};

}

#endif /* NEIGHBORHOOD_H_ */
