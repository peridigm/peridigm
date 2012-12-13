//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#include "BondFilter.h"
#include "vtkPlane.h"

using UTILITIES::Vector3D;
static UTILITIES::Dot dot;
static UTILITIES::Cross cross;
static UTILITIES::Minus minus;

namespace PdBondFilter {

FinitePlane::FinitePlane(double normal[3], double lowerLeftCorner[3], double bottom_UnitVector[3], double lengthBottom, double lengthA)
: n(normal), r0(lowerLeftCorner), ub(bottom_UnitVector), ua(cross(ub,n)), a(lengthA), b(lengthBottom)
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


void BondFilterDefault::filterBonds(std::vector<int>& treeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

	bool *flagIter = bondFlags;
	for(unsigned int n=0;n<treeList.size();n++,flagIter++){
		/*
		 * All bonds are innocent until proven guilty
		 */
		*flagIter=0;

		/*
		 * If we want to include 'x' then we do not mark
		 */
		if(includeSelf) continue;

		size_t uid = treeList[n];
		if(ptLocalId==uid) *flagIter=1;
	}

}

shared_ptr<BondFilter> BondFilterDefault::clone(bool withSelf) {
	return shared_ptr<BondFilterDefault>(new BondFilterDefault(withSelf));
}

void FinitePlaneFilter::filterBonds(std::vector<int>& treeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

	/*
	 * Create bond points
	 */
	const double *p0 = pt;
	const double *p1;
	double x[3], t;
	bool *flagIter = bondFlags;
	for(unsigned int p=0;p<treeList.size();p++,flagIter++){
		/*
		 * Local id of point within neighborhood
		 */
      size_t uid = treeList[p];
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
