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

#ifndef VECTOR_H_
#define VECTOR_H_
#include <functional>
#include <math.h>
#include <string>
#include <stdexcept>

using std::binary_function;

namespace UTILITIES {

class Vector3D {
public:
	Vector3D(){u[0]=0; u[1]=0; u[2]=0;}
	Vector3D(double a, double b, double c) { u[0]=a; u[1]=b; u[2]=c; }
	Vector3D(double v[3]) { u[0]=v[0]; u[1]=v[1]; u[2]=v[2]; }
	Vector3D(const Vector3D& rhs) { u[0]=rhs[0]; u[1]=rhs[1]; u[2]=rhs[2]; }
	double* get() { return u; }

	double norm() const {
		double dx = u[0];
		double dy = u[1];
		double dz = u[2];
		return sqrt(dx*dx+dy*dy+dz*dz);
	}

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

	Vector3D& operator*=(double scale) {
		u[0]*=scale; u[1]*=scale;; u[2]*=scale;;
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
	 * u + v
	 */
	Vector3D operator()(const Vector3D& u, const Vector3D& v){
		double r0=u[0]+v[0];
		double r1=u[1]+v[1];
		double r2=u[2]+v[2];
		Vector3D r; r[0]=r0;r[1]=r1;r[2]=r2;
		return r;
	}
};

struct InsideSphere : public binary_function< Vector3D, Vector3D, bool > {
public:
	InsideSphere(double radius) : radius(radius) {}
	bool operator()(const Vector3D& u, const Vector3D& v) const {
		double dx = v[0]-u[0];
		double dy = v[1]-u[1];
		double dz = v[2]-u[2];
		return dx*dx+dy*dy+dz*dz - radius*radius < 0.0;
	}
	double get_radius() const { return radius; }

private:
	double radius;
};

struct OutsideSphere : public binary_function< Vector3D, Vector3D, bool > {

public:
	OutsideSphere(double radius) : radius(radius) {}
	bool operator()(const Vector3D& u, const Vector3D& v){
		double dx = v[0]-u[0];
		double dy = v[1]-u[1];
		double dz = v[2]-u[2];
		return dx*dx+dy*dy+dz*dz - radius*radius > 0.0;
	}
	double get_radius() const { return radius; }

private:
	double radius;
};

inline double scalar_triple_product(const Vector3D& a, const Vector3D& b, const Vector3D& c){
	double v;
	double a11 = a[0], a12 = b[0], a13 = c[0];
	double a21 = a[1], a22 = b[1], a23 = c[1];
	double a31 = a[2], a32 = b[2], a33 = c[2];
	v=a11*(a22*a33-a23*a32)-a12*(a21*a33-a23*a31)+a13*(a21*a32-a22*a31);
	return v;
}

class PointCenteredBoundingBox {
public:
	PointCenteredBoundingBox (const double* X, double radius): x(X), horizon(radius) {}
	inline double get_xMin() const { return *(x+0)-horizon; }
	inline double get_yMin() const { return *(x+1)-horizon; }
	inline double get_zMin() const { return *(x+2)-horizon; }
	inline double get_xMax() const { return *(x+0)+horizon; }
	inline double get_yMax() const { return *(x+1)+horizon; }
	inline double get_zMax() const { return *(x+2)+horizon; }

private:
	const double *x;
	double horizon;
};

} // UTILITIES

#endif /* VECTOR_H_ */
