/*
 * PdBondFilter.h
 *
 *  Created on: Dec 17, 2010
 *      Author: jamitch
 */

#ifndef PDBONDFILTER_H_
#define PDBONDFILTER_H_
#include "vtkIdList.h"
#include <cstddef>
#include <utility>
using std::pair;
using std::size_t;

namespace PdBondFilter {


class BondFilter {
public:
	virtual ~BondFilter() {}
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) = 0;
	/*
	 * NOTE: expectation is that bondFlags has been allocated to a sufficient length so that a
	 * single scalar flag can be associated with every point in the neighborhood of 'pt';
	 * bonds are included by default, ie flag=0; if a point is excluded then flag =1 is set
	 */
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* bondFlags) = 0;
};

/**
 * This filter does NOT include the point x in its own neighborhood H(x)
 */
class BondFilterDefault : public BondFilter {
public:
	BondFilterDefault() : BondFilter() {}
	virtual ~BondFilterDefault() {}
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const double *xOverlap);
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion);

};

/**
 * This filter INCLUDES the point x in its own neighborhood H(x); Used for implicit bandwidth
 */
class BondFilterWithSelf : public BondFilter {
public:
	BondFilterWithSelf() : BondFilter() {}
	virtual ~BondFilterWithSelf() {}
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) {}
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion) {}
};

}

#endif /* PDBONDFILTER_H_ */
