/*
 * PdBondFilter.h
 *
 *  Created on: Dec 17, 2010
 *      Author: jamitch
 */

#ifndef PDBONDFILTER_H_
#define PDBONDFILTER_H_

#include <utility>
using std::pair;

namespace PdBondFilter {


class BondFilter {
public:
	virtual ~BondFilter() {}
	virtual size_t filterNumNeighbors(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) = 0;
	virtual pair<size_t,bool*> filterBonds(vtkIdList* kdTreeList, const double *pt, const double *xOverlap, bool* markForExclusion) = 0;
};

/**
 * This filter does NOT include the point x in its own neighborhood H(x)
 */
class BondFilterDefault : public BondFilter {
public:
	BondFilterDefault() : BondFilter() {}
	virtual ~BondFilterDefault() {}
	virtual size_t filterNumNeighbors(vtkIdList* kdTreeList, const double *pt, const double *xOverlap);
	virtual pair<size_t,bool*> filterBonds(vtkIdList* kdTreeList, const double *pt, const double *xOverlap, bool* markForExclusion);

};

/**
 * This filter INCLUDES the point x in its own neighborhood H(x); Used for implicit bandwidth
 */
class BondFilterWithSelf : public BondFilter {
public:
	BondFilterWithSelf() : BondFilter() {}
	virtual ~BondFilterWithSelf() {}
	virtual size_t filterNumNeighbors(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) {}
	virtual pair<size_t,bool*> filterBonds(vtkIdList* kdTreeList, const double *pt, const double *xOverlap, bool* markForExclusion) {}
};

}

#endif /* PDBONDFILTER_H_ */
