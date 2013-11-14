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

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <float.h>
#include "../../kdtree.h"

using femanica::point;
using femanica::rectangular_range;
using femanica::span_axis;
using std::numeric_limits;
typedef double value_type;
typedef std::size_t ordinal_type;


point<value_type> get_point(){
	point<value_type> a(1);
	return a;
}

TEUCHOS_UNIT_TEST( Point, ConstructorTest) {

	{
		point<value_type> a, b, c(1), d(1,1,1), e(1,2,3), f=get_point();
		point<value_type> g(c);
		TEST_ASSERT(a[femanica::X]==b[femanica::X]);
		TEST_ASSERT(a[femanica::Y]==b[femanica::Y]);
		TEST_ASSERT(a[femanica::Z]==b[femanica::Z]);
		TEST_ASSERT(a[femanica::X]==0);
		TEST_ASSERT(a[femanica::Y]==0);
		TEST_ASSERT(a[femanica::Z]==0);
		TEST_ASSERT(c[femanica::X]==f[femanica::X]);
		TEST_ASSERT(c[femanica::Y]==f[femanica::Y]);
		TEST_ASSERT(c[femanica::Z]==f[femanica::Z]);
		TEST_ASSERT(c[femanica::X]==d[femanica::X]);
		TEST_ASSERT(c[femanica::Y]==d[femanica::Y]);
		TEST_ASSERT(c[femanica::Z]==d[femanica::Z]);
		TEST_ASSERT(g[femanica::X]==f[femanica::X]);
		TEST_ASSERT(g[femanica::Y]==f[femanica::Y]);
		TEST_ASSERT(g[femanica::Z]==f[femanica::Z]);
		TEST_ASSERT(1.0==f[femanica::X]);
		TEST_ASSERT(1.0==f[femanica::Y]);
		TEST_ASSERT(1.0==f[femanica::Z]);
		TEST_ASSERT(1.0==e[femanica::X]);
		TEST_ASSERT(2.0==e[femanica::Y]);
		TEST_ASSERT(3.0==e[femanica::Z]);
		span_axis<value_type> ax={femanica::X,0.0};
		point<value_type> h=intersect(f,ax);
		TEST_ASSERT(h[femanica::X]==ax.cut);
		TEST_ASSERT(h[femanica::Y]==1.0);
		TEST_ASSERT(h[femanica::Z]==1.0);
	}

	{
		span_axis<value_type> ax={femanica::X,0.0};
		span_axis<value_type> v=ax;
		point<value_type> c(1);
		int i=v.axis;
		int j=(v.axis+1)%3;
		int k=(v.axis+2)%3;
		value_type q[3];
		q[i]=v.cut;
		q[j]=c[j];
		q[k]=c[k];
		point<value_type> a(q);
		TEST_ASSERT(a[femanica::X]==0.0);
		TEST_ASSERT(a[femanica::Y]==1.0);
		TEST_ASSERT(a[femanica::Z]==1.0);

	}

	{
		span_axis<value_type> ay={femanica::Y,0.0};
		span_axis<value_type> v=ay;
		point<value_type> c(1);
		int i=v.axis;
		int j=(v.axis+1)%3;
		int k=(v.axis+2)%3;
		value_type q[3];
		q[i]=v.cut;
		q[j]=c[j];
		q[k]=c[k];
		point<value_type> a(q);
		TEST_ASSERT(a[femanica::X]==1.0);
		TEST_ASSERT(a[femanica::Y]==0.0);
		TEST_ASSERT(a[femanica::Z]==1.0);

	}

	{
		span_axis<value_type> az={femanica::Z,0.0};
		span_axis<value_type> v=az;
		point<value_type> c(1);
		int i=v.axis;
		int j=(v.axis+1)%3;
		int k=(v.axis+2)%3;
		value_type q[3];
		q[i]=v.cut;
		q[j]=c[j];
		q[k]=c[k];
		point<value_type> a(q);
		TEST_ASSERT(a[femanica::X]==1.0);
		TEST_ASSERT(a[femanica::Y]==1.0);
		TEST_ASSERT(a[femanica::Z]==0.0);

	}

}

TEUCHOS_UNIT_TEST( Point, ComparatorTest) {

	point<value_type> a(1), b(2), c(3), g(4);
	TEST_ASSERT(a<=b);
	TEST_ASSERT(b<=c);
	TEST_ASSERT(a<=c);
	TEST_ASSERT(!(c<=b));
	point<value_type> d(1,2,3), e(0.5,1.5,2.5);
	TEST_ASSERT(e<=d);
	rectangular_range<value_type> r1(a,c), r2(point<value_type>(1.5),point<value_type>(2.5));
	TEST_ASSERT(r1.contains(r2));
	TEST_ASSERT(!r2.contains(r1));
	TEST_ASSERT(r1.intersects(r2));
	TEST_ASSERT(r2.intersects(r1));
	rectangular_range<value_type> r3(a,b), r4(c,g);
	TEST_ASSERT(!r3.contains(r4));
	TEST_ASSERT(!r4.contains(r3));
	/*
	 * r is the entire 3d domain
	 */
	rectangular_range<value_type> r;
	TEST_ASSERT(r.contains(r1));
	TEST_ASSERT(!r1.contains(r));
	span_axis<value_type> ax={femanica::X,0.0},ay={femanica::Y,0.0};
	rectangular_range<value_type> left=r.left(ax);
	rectangular_range<value_type> lower_left=left.left(ay);
	rectangular_range<value_type> right=r.right(ax);
	rectangular_range<value_type> upper_right=right.right(ay);
	TEST_ASSERT(r.contains(left));
	TEST_ASSERT(!left.contains(r));
	TEST_ASSERT(r.contains(right));
	TEST_ASSERT(!right.contains(r));
	TEST_ASSERT(r.contains(lower_left));
	TEST_ASSERT(left.contains(lower_left));
	TEST_ASSERT(!lower_left.contains(left));
	TEST_ASSERT(!lower_left.contains(r));
	TEST_ASSERT(r.contains(upper_right));
	TEST_ASSERT(right.contains(upper_right));
	TEST_ASSERT(!upper_right.contains(right));
	TEST_ASSERT(!upper_right.contains(r));

	TEST_ASSERT(upper_right.contains(r1));
	TEST_ASSERT(!r1.contains(upper_right));

	point<value_type> y(-1.0),z(-.5),w(0.0);
	rectangular_range<value_type> r5(y,z), r6(y,w);
	TEST_ASSERT(lower_left.contains(r5));
	TEST_ASSERT(left.contains(r5));
	TEST_ASSERT(r.contains(r5));
	TEST_ASSERT(!r5.contains(r));
	TEST_ASSERT(!r5.contains(r6));
	TEST_ASSERT(r6.contains(r5));
	TEST_ASSERT(r5.intersects(r6));
	TEST_ASSERT(r6.intersects(r5));

}

TEUCHOS_UNIT_TEST( Point, _RangeTest) {

	rectangular_range<value_type> r, s;
	point<value_type> rl=r.get_low(), rh=r.get_high();
	point<value_type> sl=s.get_low(), sh=s.get_high();
	TEST_ASSERT(rl[femanica::X]==-numeric_limits<value_type>::max());
	TEST_ASSERT(rl[femanica::Y]==-numeric_limits<value_type>::max());
	TEST_ASSERT(rl[femanica::Z]==-numeric_limits<value_type>::max());
	TEST_ASSERT(rh[femanica::X]==numeric_limits<value_type>::max());
	TEST_ASSERT(rh[femanica::Y]==numeric_limits<value_type>::max());
	TEST_ASSERT(rh[femanica::Z]==numeric_limits<value_type>::max());
	TEST_ASSERT(rl[femanica::X]==sl[femanica::X]);
	TEST_ASSERT(rl[femanica::Y]==sl[femanica::Y]);
	TEST_ASSERT(rl[femanica::Z]==sl[femanica::Z]);
	TEST_ASSERT(rh[femanica::X]==sh[femanica::X]);
	TEST_ASSERT(rh[femanica::Y]==sh[femanica::Y]);
	TEST_ASSERT(rh[femanica::Z]==sh[femanica::Z]);
	span_axis<value_type> ax={femanica::X,0.0};
	rectangular_range<value_type> left=r.left(ax);
	point<value_type> left_low=left.get_low();
	point<value_type> left_high=left.get_high();
	TEST_ASSERT(left_low[femanica::X]==-numeric_limits<value_type>::max());
	TEST_ASSERT(left_low[femanica::Y]==-numeric_limits<value_type>::max());
	TEST_ASSERT(left_low[femanica::Z]==-numeric_limits<value_type>::max());
	TEST_ASSERT(left_high[femanica::X]==ax.cut);
	TEST_ASSERT(left_high[femanica::Y]==numeric_limits<value_type>::max());
	TEST_ASSERT(left_high[femanica::Z]==numeric_limits<value_type>::max());
}




int main(int argc, char* argv[]) {

	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}


