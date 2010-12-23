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
#include <functional>

using std::binary_function;

using std::size_t;

namespace PdBondFilter {

class Vector3D {
public:
	Vector3D(){u[0]=0; u[1]=0; u[2]=0;}
	Vector3D(double v[3]) { u[0]=v[0]; u[1]=v[1]; u[2]=v[2]; }
	Vector3D(const Vector3D& rhs) { u[0]=rhs[0]; u[1]=rhs[1]; u[2]=rhs[2]; }
	double* get() { return u; }
	double operator[](int i) const {
		if(0>i || 3<=i){
			std::string message("ERROR\n\tVector3D::operator[](int i) const \'i\' out of range.");
			throw std::domain_error(message);
		}
		return u[i];
	  }
	double & operator[](int i) {
		if(0>i || 3<=i){
			std::string message("ERROR\n\tVector3D::operator[](int i) \'i\' out of range.");
			throw std::domain_error(message);
		}
		return u[i];
	  }
	Vector3D& operator=(const Vector3D& rhs) {
		if(&rhs == this) return *this;
		u[0]=rhs[0]; u[1]=rhs[1]; u[2]=rhs[2];
		return *this;
	}
private:
	double u[3];

};

struct Distance : public binary_function< Vector3D, Vector3D, double > {
	double operator()(const Vector3D& u, const Vector3D& v){
		double dx = v[0]-u[0];
		double dy = v[1]-u[1];
		double dz = v[2]-u[2];
		return sqrt(dx*dx+dy*dy+dz*dz);
	}
};

struct Dot : public binary_function< Vector3D, Vector3D, double >{
	double operator()(const Vector3D& u, const Vector3D& v){
		return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	}
};

struct Cross : public binary_function< Vector3D, Vector3D, Vector3D > {
	/*
	 * u 'cross' v
	 */
	Vector3D operator()(const Vector3D& u, const Vector3D& v){
		double r0=-u[2]* v[1] + u[1] * v[2];
		double r1= u[2]* v[0] - u[0] * v[2];
		double r2=-u[1]* v[0] + u[0] * v[1];
		Vector3D r; r[0]=r0;r[1]=r1;r[2]=r2;
		return r;
	}
};

struct Minus : public binary_function< Vector3D, Vector3D, Vector3D > {
	/*
	 * u - v
	 */
	Vector3D operator()(const Vector3D& u, const Vector3D& v){
		double r0=u[0]-v[0];
		double r1=u[1]-v[1];
		double r2=u[2]-v[2];
		Vector3D r; r[0]=r0;r[1]=r1;r[2]=r2;
		return r;
	}
};

struct Plus : public binary_function< Vector3D, Vector3D, Vector3D > {
	/*
	 * u - v
	 */
	Vector3D operator()(const Vector3D& u, const Vector3D& v){
		double r0=u[0]+v[0];
		double r1=u[1]+v[1];
		double r2=u[2]+v[2];
		Vector3D r; r[0]=r0;r[1]=r1;r[2]=r2;
		return r;
	}
};


class FinitePlane {
public:
	explicit FinitePlane(double normal[3], double lowerLeftCorner[3], double edgeA_UnitVector[3], double lengthA, double lengthB);
	/**
	Description:
	   Given a line defined by the two points p0,p1; and 'this' plane defined by the
	   normal n and point r0, compute an intersection (this assumes an infinite plane).
	   The parametric coordinate along the line is returned in t, and the coordinates of
	   intersection are returned in x. A zero is returned if the plane and line
	   do not intersect between (0<=t<=1). If the plane and line are parallel,
	   zero is returned and t is set to VTK_LARGE_DOUBLE.
	   */
	int bondIntersectInfinitePlane(double *p0, double *p1, double&t, double x[3]);
	/*
	 * Under the assumption that the bond (p1-p0) intersects 'this' infinite plane at
	 * the point 'x', this function determines if the intersection exists
	 * in the "INPUT" finite plane. Note that 'x' is an input parameter
	 * to this function.
	 * Calling Sequence:
	 *  if(0 != bondIntersectInfinitePlane(p0,p1,t,x))
	 *      bondIntersect(x)
	 */
	bool bondIntersect(double x[3]);
private:
	Vector3D n, r0, ua, ub;
	double a, b;
};


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
	virtual size_t filterListSize(vtkIdList* kdTreeList, const double *pt, const double *xOverlap);
	virtual void filterBonds(vtkIdList* kdTreeList, const double *pt, const size_t ptLocalId, const double *xOverlap, bool* markForExclusion);
};

}

#endif /* PDBONDFILTER_H_ */
