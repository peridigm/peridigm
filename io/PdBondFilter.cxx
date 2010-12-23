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


int FinitePlane::bondIntersectInfinitePlane(const double *p0, const double *p1, double&t, double x[3]) {
	double *non_const_p0=const_cast<double*>(p0);
	double *non_const_p1=const_cast<double*>(p1);
	return vtkPlane::IntersectWithLine(non_const_p0,non_const_p1,n.get(),r0.get(),t,x);
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

size_t BondFilterDefault::filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap) {

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

size_t BondFilterWithSelf::filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap) {

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

/**
 * THIS FILTER DOES REMOVES 'ptLocalId' -- ie does not include 'self'
 */
size_t FinitePlaneFilter::filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap) {

	/*
	 * Create bond points
	 */
	const double *p0 = pt;
	const double *p1;
	double x[3], t;
	size_t size=kdTreeList->GetNumberOfIds();
	for(size_t p=0;p<size;p++){
		/*
		 * Local id of point within neigbhorhood
		 */
		size_t uid = kdTreeList->GetId(p);

		/*
		 * This filter does not include 'self'
		 */
		if(ptLocalId==uid) {
//			cout << "FinitePlaneFilter::filterListSize DO NOT INCLUDE SELF ID = " << ptLocalId << endl;
			continue;
		}

		/*
		 * Now run plane filter
		 */
		p1 = xOverlap+(3*uid);
		if( 0 != plane.bondIntersectInfinitePlane(p0,p1,t,x) && plane.bondIntersect(x) ){
//			cout << "FinitePlaneFilter::filterListSize DO INCLUDE PT DUE TO INTERSECTION, id  = " << uid << endl;
			size -= 1;
		}
	}

	return size;
}
//
//
//void FinitePlaneFilter::filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {
//
//	/*
//	 * Create bond points
//	 */
//	const double *p0 = pt;
//	const double *p1;
//	bool *flagIter = bondFlags;
//	size_t size=kdTreeList->GetNumberOfIds();
//	for(size_t p=0;p<size;p++){
//		/*
//		 * Local id of point within neigbhorhood
//		 */
//		size_t uid = kdTreeList->GetId(p);
//		/*
//		 * All bonds are innocent until proven guilty
//		 */
//		*flagIter=0;
//		/*
//		 * This filter does not include 'self'
//		 */
//		if(ptLocalId==uid) {
//			*flagIter=1;
//			continue;
//		}
//
//		/*
//		 * Now run plane filter
//		 */
//		p1 = xOverlap+(3*uid);
//
//	}
//
//	return size;
//}

}
