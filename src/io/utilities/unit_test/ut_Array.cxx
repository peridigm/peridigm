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
#include "Array.h"

using UTILITIES::Array;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

void constructor(){

	Array<int> a,b;
	BOOST_CHECK(a.get_size()==0);
	BOOST_CHECK(b.get_size()==0);
	BOOST_CHECK(b.get_shared_ptr()==a.get_shared_ptr());
	BOOST_CHECK(shared_ptr<int>()==a.get_shared_ptr());
	std::size_t size(10);
	Array<double> c(size);
	BOOST_CHECK(c.get_size()==size);

	c.set(3.33);
	Array<double> d(size,c.get_shared_ptr());
	BOOST_CHECK(d.get_size()==size);
	for(std::size_t i=0;i<c.get_size();i++)
		BOOST_CHECK(d.get()[i]==c.get()[i]);

	Array<double> e = Array<double>(size);
	BOOST_CHECK(e.get_size()==size);
	Array<double> f = e;
	BOOST_CHECK(f.get_size()==size);

}

void set(){

	size_t size(10);
	Array<int> a(size);
	BOOST_CHECK(a.get_size()==size);
	a.set(1);
	for(std::size_t i=0;i<a.get_size();i++)
		BOOST_CHECK(1==a.get()[i]);

}

void deep_copy(){

	size_t size(10);
	Array<int> a(size),b;
	BOOST_CHECK(a.get_size()==size);
	BOOST_CHECK(b.get_size()==0);
	a.set(1);
	b.deep_copy(a);
	BOOST_CHECK(b.get_size()==size);
	for(std::size_t i=0;i<a.get_size();i++){
		BOOST_CHECK(1==b.get()[i]);
		BOOST_CHECK(b.get()[i]==a.get()[i]);
	}

}

void scale(){

	size_t size(10);
	Array<int> a(size);
	BOOST_CHECK(a.get_size()==size);
	a.set(1);
	a.scale(10);
	for(std::size_t i=0;i<a.get_size();i++){
		BOOST_CHECK(10==a.get()[i]);
	}

}




bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "ut_Array" );
	proc->add(BOOST_TEST_CASE( &constructor ));
	proc->add(BOOST_TEST_CASE( &set ));
	proc->add(BOOST_TEST_CASE( &deep_copy ));
	proc->add(BOOST_TEST_CASE( &scale ));
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
