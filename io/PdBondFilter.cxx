/*
 * PdBondFilter.cxx
 *
 *  Created on: Dec 17, 2010
 *      Author: jamitch
 */
#include "PdBondFilter.h"

namespace PdBondFilter {

size_t BondFilterDefault::filterListSize(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) {

	/* THIS RETURNS the length of the neighborhood list which is not the same as numNeighbors;
	 * In general, numNeighbors = sizeList - 1;
	 *
	 * NOTE that search results include "self"; if we want self then we need to add it
	 * to the size of the list;
	 * We already need one extra for "numNeigh" so we must add another entry for self.  On the other
	 * hand, if we don't want self, then we don't have to do anything since
	 * the existing extra entry we can re-use for "numNeigh"
	 */
	return kdTreeList->GetNumberOfIds();
}

void BondFilterDefault::filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

	bool *flagIter = bondFlags;
	for(size_t n=0;n<kdTreeList->GetNumberOfIds();n++,flagIter++){
		/*
		 * All bonds are innocent until proven guilty
		 */
		*flagIter=0;
		size_t uid = kdTreeList->GetId(n);
		if(ptLocalId==uid) *flagIter=1;
	}

}


}
