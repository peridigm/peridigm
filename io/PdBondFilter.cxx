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

FinitePlane::FinitePlane(double normal[3], double lowerLeftCorner[3], double bottom_UnitVector[3], double lengthBottom, double lengthA)
: n(normal), r0(lowerLeftCorner), ub(bottom_UnitVector), ua(cross(ub,n)), b(lengthBottom), a(lengthA)
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
//	std::cout << "\tFinitePlane::bondIntersect aa, bb = " << aa << ", " << bb << std::endl;
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
 * THIS FILTER REMOVES 'ptLocalId' -- ie does not include 'self'
 */
size_t FinitePlaneFilter::filterListSize(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap) {

//	cout << "FinitePlaneFilter::filterListSize::Evaluate Neighborhood of Local PT localId = " << ptLocalId << endl;
//	cout << "\tFinitePlaneFilter::filterListSize::Initial num points in neighborhood = " << kdTreeList->GetNumberOfIds() << endl;

	/*
	 * Create bond points
	 */
	const double *p0 = pt;
	const double *p1;
	double x[3], t;
	size_t size=kdTreeList->GetNumberOfIds();
	for(size_t p=0;p<kdTreeList->GetNumberOfIds();p++){
		/*
		 * Local id of point within neigbhorhood
		 */
		size_t uid = kdTreeList->GetId(p);
//		cout << "\tFinitePlaneFilter::filterListSize::Neighborhood PT localId = " << uid << endl;

		/*
		 * This filter does not include 'self'
		 */
		if(ptLocalId==uid) {
//			cout << "\tFinitePlaneFilter::filterListSize DO NOT INCLUDE SELF;  localId =" << uid <<  endl;
			continue;
		}

		/*
		 * Now run plane filter
		 */
		p1 = xOverlap+(3*uid);
		if( 0 != plane.bondIntersectInfinitePlane(p0,p1,t,x) && plane.bondIntersect(x) ){
//			cout << "\tFinitePlaneFilter::filterListSize INTERSECTION; localId, t  = " << uid << ", " << t << endl;
			size -= 1;
		} // else {
//			/*
//			 * ALL DEBUG PRINT STUFF
//			 */
//			cout << "\tplane.bondIntersectInfinitePlane(p0,p1,t,x) = " << plane.bondIntersectInfinitePlane(p0,p1,t,x) << endl;
//			if(plane.bondIntersectInfinitePlane(p0,p1,t,x))
//				cout << "\tplane.bondIntersect(x) = " << plane.bondIntersect(x) << endl;
//		}
	}

	return size;
}


void FinitePlaneFilter::filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

	/*
	 * Create bond points
	 */
	const double *p0 = pt;
	const double *p1;
	double x[3], t;
	bool *flagIter = bondFlags;
	for(size_t p=0;p<kdTreeList->GetNumberOfIds();p++,flagIter++){
		/*
		 * Local id of point within neigbhorhood
		 */
		size_t uid = kdTreeList->GetId(p);
		/*
		 * All bonds are innocent until proven guilty
		 */
		*flagIter=0;
		/*
		 * This filter does not include 'self'
		 */
		if(ptLocalId==uid) {
//			cout << "\tFinitePlaneFilter::filterBonds DO NOT INCLUDE SELF;  localId =" << uid <<  endl;
			*flagIter=1;
			continue;
		}

		/*
		 * Now run plane filter
		 */
		p1 = xOverlap+(3*uid);
		if( 0 != plane.bondIntersectInfinitePlane(p0,p1,t,x) && plane.bondIntersect(x) ){
//			cout << "\tFinitePlaneFilter::filterBonds DO INCLUDE PT DUE TO INTERSECTION; localId, t  = " << uid << ", " << t << endl;
			*flagIter=1;
		}
	}
}

}
