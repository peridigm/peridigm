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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "vtkPlane.h"
#include "vtkSmartPointer.h"
#include <iostream>
#include <functional>
#include <valarray>
#include <set>
#include "Vector3D.h"
#include "../BondFilter.h"

using UTILITIES::Vector3D;
using UTILITIES::Dot;
using UTILITIES::Cross;
using UTILITIES::Minus;
using namespace boost::unit_test;
using std::cout;



/*
  TEST CONSTANT arrays
 */
void myConstFunction(const double x[3]) {}
void myNonConstFunction(double x[3]){}

void callFunctions(){
	double x[3]; x[0]=0.0;x[1]=0.0;x[2]=0.0;
	const double *constXPtr = x;
	double *xPtr = x;
	myConstFunction(xPtr);
	myNonConstFunction(xPtr);
	myConstFunction(constXPtr);
	// this will not compile
//	myNonConstFunction(constXPtr);

}
/*
   END TEST CONSTANT arrays
 */


void simplePlaneCase_1(){
	/*
	 * Helper functions
	 */
	Dot dot;
	Cross cross;
	/*
	 * Lower left corner of plane (USER INPUT)
	 */
	Vector3D r0; r0[0]=0;r0[1]=0;r0[2]=0;
	/*
	 * Normal to plane (USER INPUT)
	 */
	Vector3D n; n[0]=1;n[1]=0;n[2]=0;
	/*
	 * Unit vector along edge of plane (USER INPUT)
	 */
	Vector3D ua; ua[0]=0;ua[1]=1;ua[2]=0;
	/*
	 * Unit vector along 'other' edge of plane
	 * ub = n 'cross' ua
	 */
	Vector3D ub = cross(n,ua);
	/*
	 * Assert correct cross product
	 */
	BOOST_CHECK(0.0 == ub[0]);
	BOOST_CHECK(0.0 == ub[1]);
	BOOST_CHECK(1.0 == ub[2]);

	/*
	 * bond b = p1 - p0
	 */
	double p1[3]; p1[0] =  .5; p1[1] = .5; p1[2] = .5;
	double p0[3]; p0[0] = -.5; p0[1] = .5; p0[2] = .5;

	/*
	 * Create plane and set values
	 */
	vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
	double N[3]; N[0]=n[0]; N[1]=n[1];N[2]=n[2];
	double R0[3]; R0[0]=r0[0];R0[1]=r0[1];R0[2]=r0[2];
	plane->SetOrigin(R0);
	plane->SetNormal(N);


	/*
	 * Intersect bond with plane
	 */
	double t;
	double x[3];
	double *p0Ptr=p0;
	double *p1Ptr=p1;
	double *xPtr =x;
	plane->IntersectWithLine(p0Ptr,p1Ptr,plane->GetNormal(),plane->GetOrigin(),t,xPtr);
	/*
	 * Assert that intersection exists and that its at the center of the bond
	 */
	BOOST_CHECK(0.5 == t);
	BOOST_CHECK(0.0==x[0]);
	BOOST_CHECK(0.5==x[1]);
	BOOST_CHECK(0.5==x[2]);

	/*
	 * Does the point of intersection exist within the finite plane?
	 */
	Vector3D r; r[0]=x[0]; r[1]=x[1]; r[2]=x[2];
	Minus minus;
	Vector3D dr = minus(r,r0);
	BOOST_CHECK(0.0==dr[0]);
	BOOST_CHECK(0.5==dr[1]);
	BOOST_CHECK(0.5==dr[2]);

	// Dot 'dr' onto 'ua' and 'ub'
	double alpha = dot(dr,ua);
	double beta = dot(dr,ub);

	// if(alpha <= a &&  beta<= b) then bond intersects input plane OTHERWISE NOT
	BOOST_CHECK(0.5==alpha);
	BOOST_CHECK(0.5==beta);

}

void simplePlaneCase_2(){
	/*
	 * Lower left corner of plane (USER INPUT)
	 */
	double r0[3]; r0[0]=0;r0[1]=1;r0[2]=0;
	/*
	 * Normal to plane (USER INPUT)
	 */
	double n[3]; n[0]=1;n[1]=0;n[2]=0;
	/*
	 * Unit vector along edge of plane (USER INPUT)
	 */
	double ua[3]; ua[0]=0;ua[1]=-1;ua[2]=0;
	/*
	 * Plane dimension along user input edge (USER INPUT)
	 */
	double a = 1.0;
	/*
	 * Plane dimension along edge perpendicular to plane normal and user input edge (USER INPUT)
	 */
	double b = 1.0;

	PdBondFilter::FinitePlane plane(n,r0,ua,a,b);
	double p1[3]; p1[0] =  .5; p1[1] = .5; p1[2] = .5;
	double p0[3]; p0[0] = -.5; p0[1] = .5; p0[2] = .5;
	double *p1Ptr=p1;
	double *p0Ptr=p0;
	double x[3];
	double t;
	BOOST_CHECK(0!=plane.bondIntersectInfinitePlane(p0Ptr,p1Ptr,t,x));
	/*
	 * Assert that intersection exists and that its at the center of the bond
	 */
	BOOST_CHECK(0.5 == t);
	BOOST_CHECK(0.0==x[0]);
	BOOST_CHECK(0.5==x[1]);
	BOOST_CHECK(0.5==x[2]);

	/*
	 * Assert plane function that checks for bond intersection
	 */
	BOOST_CHECK(true==plane.bondIntersect(x));

	/*
	 * Create a 2nd bond; Should find NO intersection
	 */
	p1[0] = -.5; p1[1] = 1.5; p1[2] = .5;
	p0[0] = -.5; p0[1] = .5; p0[2] = .5;
	BOOST_CHECK(0==plane.bondIntersectInfinitePlane(p0Ptr,p1Ptr,t,x));
}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "utPdVtkPlane" );
	proc->add(BOOST_TEST_CASE( &simplePlaneCase_1 ));
	proc->add(BOOST_TEST_CASE( &simplePlaneCase_2 ));
	framework::master_test_suite().add( proc );
	return success;
}


bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(
		int argc,
		char* argv[]
)
{


	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}

