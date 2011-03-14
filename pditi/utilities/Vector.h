/*
 * VectorUtils.h
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

#ifndef VECTOR_H_
#define VECTOR_H_
#include <functional>
#include <math.h>

using std::binary_function;

namespace UTILITIES {

class Vector3D {
public:
	Vector3D(){u[0]=0; u[1]=0; u[2]=0;}
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
	InsideSphere(double radius) : r(radius) {}
	bool operator()(const Vector3D& u, const Vector3D& v){
		double dx = v[0]-u[0];
		double dy = v[1]-u[1];
		double dz = v[2]-u[2];
		return dx*dx+dy*dy+dz*dz - r*r < 0.0;
	}

private:
	double r;
};

struct OutsideSphere : public binary_function< Vector3D, Vector3D, bool > {

public:
	OutsideSphere(double radius) : r(radius) {}
	bool operator()(const Vector3D& u, const Vector3D& v){
		double dx = v[0]-u[0];
		double dy = v[1]-u[1];
		double dz = v[2]-u[2];
		return dx*dx+dy*dy+dz*dz - r*r > 0.0;
	}

private:
	double r;
};



} // UTILITIES

#endif /* VECTOR_H_ */
