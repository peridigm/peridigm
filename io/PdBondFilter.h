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
#include <stdexcept>
#include <string>
#include "VectorUtils.h"

using std::size_t;

namespace PdBondFilter {



class FinitePlane {
public:
	/**
	 * @param normal: unit vector normal to plane
	 * @param lowerLeftCorner: looking down the normal (+dir), this is the lower left hand corner coordinate of the plane
	 * @param bottom_UnitVector: unit vector along bottom edge of plane; tail is at 'lowerLeftCorner' (IMPORTANT)
	 * @param lengthBottom: length of the bottom edge associated with 'bottom_UnitVector'
	 * @param lengthA: length of 2nd edge; assumption is that second edge is perpendicular to both the normal and bottom edge;
	 * The positive direction for edge A is taken along 'bottom_UnitVector cross 'normal'
	 */
	explicit FinitePlane(double normal[3], double lowerLeftCorner[3], double bottom_UnitVector[3], double lengthBottom, double lengthA);
	/**
	Description:
	   Given a line defined by the two points p0,p1; and 'this' plane defined by the
	   normal n and point r0, compute an intersection (this assumes an infinite plane).
	   The parametric coordinate along the line is returned in t, and the coordinates of
	   intersection are returned in x. A zero is returned if the plane and line
	   do not intersect between (0<=t<=1). If the plane and line are parallel,
	   zero is returned and t is set to VTK_LARGE_DOUBLE.
	   */
	int bondIntersectInfinitePlane(const double *p0, const double *p1, double&t, double x[3]);
	/*
	 * Under the assumption that the bond (p1-p0) intersects 'this' infinite plane at
	 * the point 'x', this function determines if the intersection exists
	 * in the "INPUT" finite plane. Note that 'x' is an input parameter
	 * to this function.
	 * Calling Sequence:
	 *  if(0 != bondIntersectInfinitePlane(p0,p1,t,x))
	 *      bondIntersect(x)
	 */
	bool bondIntersect(double x[3], double tolerance=1.0e-15);
private:
	VectorUtilsNS::Vector3D n, r0, ub, ua;
	double a, b;
};


class BondFilter {
public:
	virtual ~BondFilter() {}
	/*
	 * This list length includes an entry for the number of neighbors : "numNeigh"
	 */
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap) = 0;
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
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap);
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion);

};

/**
 * This filter INCLUDES the point x in its own neighborhood H(x); Used for implicit bandwidth
 */
class BondFilterWithSelf : public BondFilter {
public:
	BondFilterWithSelf() : BondFilter() {}
	virtual ~BondFilterWithSelf() {}
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap);
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
};

/**
 * Filter removes bonds from Neighborhood that intersect with "FinitePlane";
 * NOTE: This filter does NOT include the point x in its own neighborhood H(x)
 */
class FinitePlaneFilter: public BondFilter {
public:
	FinitePlaneFilter(const FinitePlane& plane, double tolerance=1.0e-15) : BondFilter(), tolerance(tolerance), plane(plane) {}
	virtual ~FinitePlaneFilter() {}
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap);
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
protected:
	double tolerance;
	FinitePlane plane;
};

/**
 * Filter removes bonds from Neighborhood that intersect with "FinitePlane";
 * NOTE: This filter DOES INCLUDE 'x' in the neighborhood H(x)
 */
class FinitePlaneFilterWithSelf: public FinitePlaneFilter {
public:
//	FinitePlaneFilterWithSelf(const FinitePlane& plane, double tolerance=1.0e-15) : BondFilter(), tolerance(tolerance), plane(plane) {}
	FinitePlaneFilterWithSelf(const FinitePlane& plane, double tolerance=1.0e-15) : FinitePlaneFilter(plane, tolerance) {}
	virtual ~FinitePlaneFilterWithSelf() {}
	/*
	 * This function is implemented in the base class -- turns out that this function really only needs to
	 * compute an upper bound estimate neighborhood list length.  Its only later during the
	 * 'filterBonds' call that the actual list length is finally determined
	 */
//	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap);
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
};

}

#endif /* PDBONDFILTER_H_ */
