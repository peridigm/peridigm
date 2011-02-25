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

class Epetra_Comm;

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
	NeighborhoodList(size_t numOwnedPoints, shared_ptr<int>& ownedGIDs, shared_ptr<double>& owned_coordinates, double horizon);
	const Epetra_BlockMap getOverlapMap(const Epetra_Comm& comm, int ndf) const;
	const Epetra_BlockMap getOwnedMap(const Epetra_Comm& comm, int ndf) const;
	shared_ptr< std::set<int> > constructParallelDecompositionFrameSet() const;

private:
	pair<int, shared_ptr<int> > getSharedGlobalIds() const;


private:
	size_t num_owned_points, size_neighborhood_list;
	double horizon;
	shared_ptr<int> owned_gids;
	shared_ptr<double> owned_x;
	shared_ptr<int> neighborhood, neighborhood_ptr;

};

}

#endif /* NEIGHBORHOOD_H_ */
