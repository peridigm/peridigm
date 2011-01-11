/*
 * Q2Cylinder.cpp
 *
 *  Created on: Jan 4, 2011
 *      Author: jamitch
 */

#include "Q2Cylinder.hpp"
#include <Epetra_Vector.h>
#include <string>

using std::string;


namespace InitialConditionsNS {


Q2Cylinder::Q2Cylinder(double _vr0, double _vr1, double _vz0, double _z0, double _a, const VectorUtilsNS::Vector3D& _center):
		vr0(_vr0), vr1(_vr1), vz0(_vz0), z0(_z0), a(_a), center(_center)
{}

void Q2Cylinder::apply(const Epetra_Vector& X, Epetra_Vector& u, Epetra_Vector& v) {
	int numPoints = X.Map().NumMyElements();

	for(int p=0;p<numPoints;p++){
		int ptr = 3*p;
		int iX = ptr;
		int iY = ptr + 1;
		int iZ = ptr + 2;

		double x = X[iX];
		double y = X[iY];
		double z = X[iZ] - z0 - a;

		double vr = vr0 - vr1*(z/a)*(z/a);
		double vz = vz0*(z/a);
		double vtheta = 0.0;
		double theta = atan2(y, x);
		v[iX] = vr*cos(theta) - vtheta*sin(theta);
		v[iY] = vr*sin(theta) + vtheta*cos(theta);
		v[iZ] = vz;
	}

}

} // InitialConditionsNS



