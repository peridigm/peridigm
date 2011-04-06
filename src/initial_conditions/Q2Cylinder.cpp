/*! \file Q2Cylinder.cpp */

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
