/*
 * PdBondFilter.cxx
 *
 *  Created on: Dec 17, 2010
 *      Author: jamitch
 */
#include "PdBondFilter.h"
#include "vtkPlane.h"

namespace PdBondFilter {

static Dot dot;
static Cross cross;
static Minus minus;
static Plus plus;

FinitePlane::FinitePlane(double normal[3], double lowerLeftCorner[3], double edgeA_UnitVector[3], double lengthA, double lengthB)
: n(normal), r0(lowerLeftCorner), ua(edgeA_UnitVector), ub(cross(n,ua)),  a(lengthA), b(lengthB)
{}


int FinitePlane::bondIntersectInfinitePlane(double *p0, double *p1, double&t, double x[3]) {
	return vtkPlane::IntersectWithLine(p0,p1,n.get(),r0.get(),t,x);
}

bool FinitePlane::bondIntersect(double x[3]) {
	bool intersects = false;
	Vector3D r(x);
	Vector3D dr(minus(r,r0));
	double aa=dot(dr,ua);
	double bb=dot(dr,ub);
	if(0<=aa && aa<=a && 0<=bb && bb<=b)
		intersects=true;
	return intersects;
}

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

size_t BondFilterWithSelf::filterListSize(vtkIdList* kdTreeList, const double *pt, const double *xOverlap) {

	/* THIS RETURNS the length of the neighborhood list which is not the same as numNeighbors;
	 * In general, numNeighbors = sizeList - 1;
	 *
	 * NOTE that search results include "self"; if we want self then we need to add it
	 * to the size of the list;
	 * We already need one extra for "numNeigh" so we must add another entry for self.
	 */
	return kdTreeList->GetNumberOfIds()+1;
}

void BondFilterWithSelf::filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

	bool *flagIter = bondFlags;
	for(size_t n=0;n<kdTreeList->GetNumberOfIds();n++,flagIter++){
		/*
		 * All bonds are innocent until proven guilty; in this case they are all 'innocent'
		 */
		*flagIter=0;
	}

}


}
