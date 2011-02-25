/*
 * Neighborhood.cxx
 *
 *  Created on: Feb 25, 2011
 *      Author: jamitch
 */
#include "NeighborhoodList.h"
#include "PdNeighborhood.h"
#include "Epetra_Comm.h"

namespace PDNEIGH {

/*
 * Prototype for private function
 */
const Epetra_BlockMap getOverlap(int ndf, int numShared, const int* shared, int numOwned,const int* owned, const Epetra_Comm& comm);

NeighborhoodList::NeighborhoodList(std::size_t numOwnedPoints, shared_ptr<int>& ownedGIDs, shared_ptr<double>& owned_coordinates, double horizon)
:
		num_owned_points(numOwnedPoints),
		size_neighborhood_list(1),
		horizon(horizon),
		owned_gids(ownedGIDs),
		owned_x(owned_coordinates),
		neighborhood(new int[1],ArrayDeleter<int>()),
		neighborhood_ptr(new int[numOwnedPoints],ArrayDeleter<int>())
{
}

pair<int, shared_ptr<int> > NeighborhoodList::getSharedGlobalIds() const {
	std::set<int> ownedIds(owned_gids.get(),owned_gids.get()+num_owned_points);
	std::set<int> shared;
	const int *neighPtr = neighborhood_ptr.get();
	const int *neigh = neighborhood.get();
	std::set<int>::const_iterator ownedIdsEnd = ownedIds.end();
	for(int p=0;p<num_owned_points;p++){
		int ptr = neighPtr[p];
		int numNeigh = neigh[ptr];
		for(int n=1;n<=numNeigh;n++){
			int id = neigh[ptr+n];
			/*
			 * look for id in owned points
			 */
			if(ownedIdsEnd == ownedIds.find(id)){
				/*
				 * add this point to shared
				 */
				shared.insert(id);
			}
		}
	}

	// Copy set into shared ptr
	shared_ptr<int> sharedGlobalIds(new int[shared.size()],ArrayDeleter<int>());
	int *sharedPtr = sharedGlobalIds.get();
	std::set<int>::iterator it;
	for ( it=shared.begin() ; it != shared.end(); it++, sharedPtr++ )
		*sharedPtr = *it;

	return pair<int, shared_ptr<int> >(shared.size(),sharedGlobalIds);
}


const Epetra_BlockMap NeighborhoodList::getOwnedMap(const Epetra_Comm& comm, int ndf) const {
	int numShared=0;
	const int *sharedPtr=NULL;
	int numOwned = num_owned_points;
	const int *ownedPtr = owned_gids.get();
	return getOverlap(ndf, numShared,sharedPtr,numOwned,ownedPtr,comm);
}

const Epetra_BlockMap NeighborhoodList::getOverlapMap(const Epetra_Comm& comm, int ndf) const {
	std::pair<int, std::tr1::shared_ptr<int> > sharedPair = this->getSharedGlobalIds();
	shared_ptr<int> sharedPtr = sharedPair.second;
	int numShared = sharedPair.first;
	const int *shared = sharedPtr.get();
	const int *owned = owned_gids.get();
	int numOwned = num_owned_points;
	return getOverlap(ndf,numShared,shared,numOwned,owned,comm);
}


shared_ptr< std::set<int> > NeighborhoodList::constructParallelDecompositionFrameSet() const {
	shared_ptr<double> xPtr = owned_x;
	shared_ptr<int> gIdsPtr = owned_gids;
	int numCells = num_owned_points;
	const PdNeighborhood::Coordinates c(xPtr,numCells);

	int numAxes=3;
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y, PdNeighborhood::Z};
	std::vector<std::tr1::shared_ptr<int> > sortedMaps(numAxes);
	for(int j=0;j<numAxes;j++){

		PdNeighborhood::Coordinates::SortComparator compare = c.getSortComparator(labels[j]);
		sortedMaps[j] = c.getIdentityMap();
		/*
		 * Sort points
		 */
		std::sort(sortedMaps[j].get(),sortedMaps[j].get()+numCells,compare);

	}

	/*
	 * Loop over axes and collect points at min and max ranges
	 * Add Points to frame set
	 */
	shared_ptr< std::set<int> >  frameSetPtr(new std::set<int>);
	PdNeighborhood::CoordinateLabel *label = labels;
	PdNeighborhood::CoordinateLabel *endLabel = labels+numAxes;
	for(;label!=endLabel;label++){

		{
			/*
			 * MINIMUM
			 * Find least upper bound of points for Min+horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			std::tr1::shared_ptr<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();
			/*
			 * First value in map corresponds with minimum value
			 */
			int iMIN = *map;
			double min = x[3*iMIN+axis];
			double value = min + horizon;
			const PdNeighborhood::Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const PdNeighborhood::Coordinates::SearchIterator end=start+numCells;
			PdNeighborhood::Coordinates::SearchIterator lub = std::upper_bound(start,end,value);
			frameSetPtr->insert(lub.mapStart(),lub.mapIterator());

		}

		{
			/*
			 * MAXIMUM
			 * Find greatest lower bound glb for Max-horizon
			 */
			int axis = *label;
			double *x = xPtr.get();
			std::tr1::shared_ptr<int> mapPtr = sortedMaps[axis];
			int *map = mapPtr.get();

			/*
			 * Last value in map corresponds with maximum value
			 */
			int iMAX = *(map+numCells-1);
			double max = x[3*iMAX+axis];
			double value = max - horizon;
			const PdNeighborhood::Coordinates::SearchIterator start=c.begin(*label,mapPtr);
			const PdNeighborhood::Coordinates::SearchIterator end=start+numCells;
			PdNeighborhood::Coordinates::SearchIterator glb = std::upper_bound(start,end,value);
			frameSetPtr->insert(glb.mapIterator(),glb.mapEnd());
		}
	}

	return frameSetPtr;

}


/*
 * PRIVATE FUNCTION
 */
const Epetra_BlockMap getOverlap(int ndf, int numShared, const int* shared, int numOwned, const int* owned, const Epetra_Comm& comm){

	int numPoints = numShared+numOwned;
	shared_ptr<int> ids(new int[numPoints],ArrayDeleter<int>());
	int *ptr = ids.get();

	for(int j=0;j<numOwned;j++,ptr++)
		*ptr=owned[j];

	for(int j=0;j<numShared;j++,ptr++)
		*ptr=shared[j];

	return Epetra_BlockMap(-1,numPoints, ids.get(),ndf, 0,comm);
}


}
