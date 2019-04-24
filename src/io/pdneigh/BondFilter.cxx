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
#include <cmath>
#include <float.h>

#include <iostream>

using UTILITIES::Vector3D;
static UTILITIES::Dot dot;
static UTILITIES::Cross cross;
static UTILITIES::Minus minus;

namespace PdBondFilter {

FinitePlane::FinitePlane(double normal[3], double lowerLeftCorner[3], double bottom_UnitVector[3], double lengthBottom, double lengthA)
: n(normal), r0(lowerLeftCorner), ub(bottom_UnitVector), ua(cross(ub,n)), a(lengthA), b(lengthBottom)
{}

int FinitePlane::bondIntersectInfinitePlane(const double *p0, const double *p1, double&t, double x[3]) {

  double numerator   = (r0[0] - p0[0]) * n[0] + (r0[1] - p0[1]) * n[1] + (r0[2] - p0[2]) * n[2];
  double denominator = (p1[0] - p0[0]) * n[0] + (p1[1] - p0[1]) * n[1] + (p1[2] - p0[2]) * n[2];

  double tolerance = 1.0e-14;

  if(std::abs(denominator) < tolerance){
    // line is parallel to plane
    // it may or may not lie on the plane
    // if it does lie on the plane, then the numerator will be zero
    // in either case, this function will return "no intersection"
    t = DBL_MAX;
  }
  else{
    // the line intersects the plane
    t = numerator/denominator;
  }

  // determine if the line segment intersects the plane
  int intersection = 0;
  if(t >= 0.0 && t <= 1.0)
    intersection = 1;

  // intersection point
  x[0] = p0[0] + t * (p1[0] - p0[0]);
  x[1] = p0[1] + t * (p1[1] - p0[1]);
  x[2] = p0[2] + t * (p1[2] - p0[2]);

  return intersection;
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
		 * If we want to include 'x' then we do not mark
		 */
		if(includeSelf) continue;

		size_t uid = treeList[n];
		if(ptLocalId==uid) *flagIter=1;
	}

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
		 * We mark points only if we do not want to include them
		 */
		if(ptLocalId==uid && !includeSelf) {
			*flagIter=1;
			continue;
		}

		/*
		 * Now run plane filter
		 */
		p1 = xOverlap+(3*uid);
		if( 0 != plane.bondIntersectInfinitePlane(p0,p1,t,x) && plane.bondIntersect(x,tolerance) ){
			*flagIter=1;
		}
	}
}

void DiskFilter::filterBonds(std::vector<int>& treeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

  const double *p0 = pt;
  const double *p1;
  bool *flagIter = bondFlags;
  for(unsigned int p=0;p<treeList.size();p++,flagIter++){

    // Local id of point within neighborhood
    size_t uid = treeList[p];

    // Set flag for bonds that will be excluded from the neighborlist
    p1 = xOverlap+(3*uid);
    if(ptLocalId==uid && !includeSelf) {
      *flagIter=1;
      continue;
    }
    if( bondIntersectsDisk(p0, p1) ) {
      *flagIter=1;
    }
  }
}

bool DiskFilter::bondIntersectsDisk(const double* p0, const double* p1) const {

  double numerator   = (center[0] - p0[0]) * normal[0] + (center[1] - p0[1]) * normal[1] + (center[2] - p0[2]) * normal[2];
  double denominator = (p1[0] - p0[0]) * normal[0] + (p1[1] - p0[1]) * normal[1] + (p1[2] - p0[2]) * normal[2];

  double t;

  if(std::abs(denominator) < tolerance){
    // line is parallel to plane
    // it may or may not lie on the plane
    // if it does lie on the plane, then the numerator will be zero
    // in either case, this function will return "no intersection"
    t = DBL_MAX;
  }
  else{
    // the line intersects the plane
    t = numerator/denominator;
  }

  if (t < 0.0 || t > 1.0)
    return false;

  // intersection point
  double x[3];
  x[0] = p0[0] + t * (p1[0] - p0[0]);
  x[1] = p0[1] + t * (p1[1] - p0[1]);
  x[2] = p0[2] + t * (p1[2] - p0[2]);

  // check if intesection point is within disk
  double distance_squared = (x[0] - center[0])*(x[0] - center[0]) +
                            (x[1] - center[1])*(x[1] - center[1]) +
                            (x[2] - center[2])*(x[2] - center[2]);

  if (distance_squared < radius*radius)
    return true;

  return false;
}

void TriangleFilter::filterBonds(std::vector<int>& treeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool *bondFlags) {

  const double *p0 = pt;
  const double *p1;
  bool *flagIter = bondFlags;
  for(unsigned int p=0;p<treeList.size();p++,flagIter++){

    // Local id of point within neighborhood
    size_t uid = treeList[p];

    // Set flag for bonds that will be excluded from the neighborlist
    p1 = xOverlap+(3*uid);
    if(ptLocalId==uid && !includeSelf) {
      *flagIter=1;
      continue;
    }
    if( bondIntersectsTriangle(p0, p1) ) {
      *flagIter=1;
    }
  }
}

bool TriangleFilter::bondIntersectsTriangle(const double* p0, const double* p1) const {

  double numerator   = (v1_[0] - p0[0]) * normal_[0] + (v1_[1] - p0[1]) * normal_[1] + (v1_[2] - p0[2]) * normal_[2];
  double denominator = (p1[0] - p0[0]) * normal_[0] + (p1[1] - p0[1]) * normal_[1] + (p1[2] - p0[2]) * normal_[2];

  double t;

  if(std::abs(denominator) < tolerance_){
    // line is parallel to plane
    // it may or may not lie on the plane
    // if it does lie on the plane, then the numerator will be zero
    // in either case, this function will return "no intersection"
    t = DBL_MAX;
  }
  else {
    // the line intersects the plane
    t = numerator/denominator;
  }

  if (t < 0.0 || t > 1.0)
    return false;

  // intersection point
  double x[3];
  x[0] = p0[0] + t * (p1[0] - p0[0]);
  x[1] = p0[1] + t * (p1[1] - p0[1]);
  x[2] = p0[2] + t * (p1[2] - p0[2]);

  // determine if the point is within the triangle
  return pointInTriangle(x);
}

bool TriangleFilter::pointInTriangle(const double* x) const {

  // check if intesection point is within triangle by computing the barycentric coordinates
  double a[3], b[3], c[3];
  for (int i=0 ; i<3 ; i++) {
    a[i] = v3_[i] - v1_[i];
    b[i] = v2_[i] - v1_[i];
    c[i] = x[i]   - v1_[i];
  }
  double dot00 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
  double dot01 = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
  double dot02 = a[0]*c[0] + a[1]*c[1] + a[2]*c[2];
  double dot11 = b[0]*b[0] + b[1]*b[1] + b[2]*b[2];
  double dot12 = b[0]*c[0] + b[1]*c[1] + b[2]*c[2];
  double alpha = (dot11 * dot02 - dot01 * dot12) / (dot00 * dot11 - dot01 * dot01);
  double beta  = (dot00 * dot12 - dot01 * dot02) / (dot00 * dot11 - dot01 * dot01);
  bool in_triangle = (alpha > -tolerance_) && (beta > -tolerance_) && (alpha + beta < 1.0 + 2*tolerance_);
  return in_triangle;
}

} // namespace PdBondFilter
