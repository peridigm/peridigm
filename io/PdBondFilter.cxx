/*
 * PdBondFilter.cxx
 *
 *  Created on: Dec 17, 2010
 *      Author: jamitch
 */
#include "PdBondFilter.h"
#include "vtkPlane.h"

using VectorUtilsNS::Vector3D;
static VectorUtilsNS::Dot dot;
static VectorUtilsNS::Cross cross;
static VectorUtilsNS::Minus minus;
static VectorUtilsNS::Plus plus;

namespace PdBondFilter {

FinitePlane::FinitePlane(double normal[3], double lowerLeftCorner[3], double bottom_UnitVector[3], double lengthBottom, double lengthA)
: n(normal), r0(lowerLeftCorner), ub(bottom_UnitVector), ua(cross(ub,n)), b(lengthBottom), a(lengthA)
{}


int FinitePlane::bondIntersectInfinitePlane(const double *p0, const double *p1, double&t, double x[3]) {
	double *non_const_p0=const_cast<double*>(p0);
	double *non_const_p1=const_cast<double*>(p1);
	return vtkPlane::IntersectWithLine(non_const_p0,non_const_p1,n.get(),r0.get(),t,x);
}

bool FinitePlane::bondIntersect(double x[3], double tolerance) {
	double zero=tolerance;
	double one = 1.0+zero;
	bool intersects = false;
	Vector3D r(x);
	Vector3D dr(minus(r,r0));
	double aa=dot(dr,ua);
	double bb=dot(dr,ub);
	if(-zero<aa && aa/a<one && -zero<bb && bb/b<one)
		intersects=true;
	return intersects;
}


void BondFilterDefault::filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

	bool *flagIter = bondFlags;
	for(size_t n=0;n<kdTreeList->GetNumberOfIds();n++,flagIter++){
		/*
		 * All bonds are innocent until proven guilty
		 */
		*flagIter=0;

		/*
		 * If we want to include 'x' then we do not mark
		 */
		if(includeSelf) continue;

		size_t uid = kdTreeList->GetId(n);
		if(ptLocalId==uid) *flagIter=1;
	}

}

shared_ptr<BondFilter> BondFilterDefault::clone(bool withSelf) {
	return shared_ptr<BondFilterDefault>(new BondFilterDefault(withSelf));
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
		 * Local id of point within neighborhood
		 */
		size_t uid = kdTreeList->GetId(p);
		/*
		 * All bonds are innocent until proven guilty
		 */
		*flagIter=0;
		/*
		 * We mark points only if we do not want to include them
		 */
		if(ptLocalId==uid && !includeSelf) {
//			cout << "\tFinitePlaneFilter::filterBonds DO NOT INCLUDE SELF;  localId =" << uid <<  endl;
			*flagIter=1;
			continue;
		}

		/*
		 * Now run plane filter
		 */
		p1 = xOverlap+(3*uid);
		if( 0 != plane.bondIntersectInfinitePlane(p0,p1,t,x) && plane.bondIntersect(x,tolerance) ){
//			cout << "\tFinitePlaneFilter::filterBonds DO INCLUDE PT DUE TO INTERSECTION; localId, t  = " << uid << ", " << t << endl;
			*flagIter=1;
		}
	}
}

shared_ptr<BondFilter> FinitePlaneFilter::clone(bool withSelf) {
	return shared_ptr<FinitePlaneFilter>(new FinitePlaneFilter(plane,withSelf,tolerance));
}

}
