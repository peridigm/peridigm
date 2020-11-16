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

#ifndef BONDFILTER_H_
#define BONDFILTER_H_
#include <cstddef>
#include <utility>
#include <stdexcept>
#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include "Vector3D.h"

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
     zero is returned and t is set to DBL_MAX.
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
  UTILITIES::Vector3D n, r0, ub, ua;
  double a, b;
};

class BondFilter {
public:
  BondFilter(bool withSelf=false) : includeSelf(withSelf) {}
  virtual ~BondFilter()  {}
  /*
   * NOTE: expectation is that bondFlags has been allocated to a sufficient length so that a
   * single scalar flag can be associated with every point in the neighborhood of 'pt';
   * bonds are included by default, ie flag=0; if a point is excluded then flag =1 is set
   */
  virtual void filterBonds(std::vector<int>& treeList, const double *pt, const std::size_t ptLocalId, const double *xOverlap, bool* bondFlags) = 0;
protected:
  bool includeSelf;
};

class BondFilterDefault : public BondFilter {
public:
  BondFilterDefault(bool withSelf=false) : BondFilter(withSelf) {}
  virtual ~BondFilterDefault() {}
  virtual void filterBonds(std::vector<int>& treeList, const double *pt, const std::size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
};

/**
 * Filter removes bonds from Neighborhood that intersect with "FinitePlane";
 * NOTE: This filter does NOT include the point x in its own neighborhood H(x)
 */
class FinitePlaneFilter: public BondFilter {
public:
  FinitePlaneFilter(const FinitePlane& plane) : BondFilter(false), tolerance(1.0e-15),   plane(plane) {}
  FinitePlaneFilter(const FinitePlane& plane, bool withSelf) : BondFilter(withSelf), tolerance(1.0e-15),   plane(plane) {}
  virtual ~FinitePlaneFilter() {}
  virtual void filterBonds(std::vector<int>& treeList, const double *pt, const std::size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
private:
  double tolerance;
  FinitePlane plane;
};

class DiskFilter: public BondFilter {
public:
  DiskFilter(double *c, double* n, double r) : BondFilter(false), tolerance(1.0e-15), radius(r) {
    center[0] = c[0];
    center[1] = c[1];
    center[2] = c[2];
    normal[0] = n[0];
    normal[1] = n[1];
    normal[2] = n[2];
  }
  virtual ~DiskFilter() {}
  virtual void filterBonds(std::vector<int>& treeList, const double *pt, const std::size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
private:
  bool bondIntersectsDisk(const double* p0, const double* p1) const;
  double tolerance;
  double center[3];
  double normal[3];
  double radius;
};

class TriangleFilter: public BondFilter {
public:
  TriangleFilter(double *v1, double* v2, double* v3) : BondFilter(false), tolerance_(1.0e-14) {
    for (int i=0 ; i<3 ; i++) {
      v1_[i] = v1[i];
      v2_[i] = v2[i];
      v3_[i] = v3[i];
    }
    double a[3], b[3];
    a[0] = v2_[0] - v1_[0];
    a[1] = v2_[1] - v1_[1];
    a[2] = v2_[2] - v1_[2];
    b[0] = v3_[0] - v1_[0];
    b[1] = v3_[1] - v1_[1];
    b[2] = v3_[2] - v1_[2];
    normal_[0] = (a[1]*b[2] - a[2]*b[1]);
    normal_[1] = (a[2]*b[0] - a[0]*b[2]);
    normal_[2] = (a[0]*b[1] - a[1]*b[0]);
    double normal_magnitude = std::sqrt(normal_[0]*normal_[0] + normal_[1]*normal_[1] + normal_[2]*normal_[2]);
    normal_[0] /= normal_magnitude;
    normal_[1] /= normal_magnitude;
    normal_[2] /= normal_magnitude;
  }
  virtual ~TriangleFilter() {}
  virtual void filterBonds(std::vector<int>& treeList, const double *pt, const std::size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
private:
  bool bondIntersectsTriangle(const double* p0, const double* p1) const;
  bool pointInTriangle(const double* x) const;
  double v1_[3];
  double v2_[3];
  double v3_[3];
  double normal_[3];
  double tolerance_;
};

}

#endif /* BONDFILTER_H_ */
