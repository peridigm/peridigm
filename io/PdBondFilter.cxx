/*
 * PdBondFilter.cxx
 *
 *  Created on: Dec 17, 2010
 *      Author: jamitch
 */
#include "PdBondFilter.h"

namespace PdBondFilter {

size_t BondFilterDefault::filterNumNeighbors(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) {

}

pair<size_t,bool*> BondFilterDefault::filterBonds(vtkIdList* kdTreeList, const double *pt, const double *xOverlap, bool *markNeighborsForExclusion) {
	size_t numNeigh = 0;
	return make_pair<size_t,bool*>(numNeigh,markNeighborsForExclusion);
}


}
