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

#include "Array.h"

using UTILITIES::Array;
using std::shared_ptr;


TEUCHOS_UNIT_TEST(Array, ConstructorTest) {


	Array<int> a,b;
	TEST_ASSERT(a.get_size()==0);
	TEST_ASSERT(b.get_size()==0);
	TEST_ASSERT(b.get_shared_ptr()==a.get_shared_ptr());
	TEST_ASSERT(shared_ptr<int>()==a.get_shared_ptr());
	std::size_t size(10);
	Array<double> c(size);
	TEST_ASSERT(c.get_size()==size);

	c.set(3.33);
	Array<double> d(size,c.get_shared_ptr());
	TEST_ASSERT(d.get_size()==size);
	for(std::size_t i=0;i<c.get_size();i++)
		TEST_ASSERT(d.get()[i]==c.get()[i]);

	Array<double> e = Array<double>(size);
	TEST_ASSERT(e.get_size()==size);
	Array<double> f = e;
	TEST_ASSERT(f.get_size()==size);

}

TEUCHOS_UNIT_TEST(Array, SetTest) {

	size_t size(10);
	Array<int> a(size);
        TEST_ASSERT(a.get_size()==size);
	a.set(1);
	for(std::size_t i=0;i<a.get_size();i++)
		TEST_ASSERT(1==a.get()[i]);

}

TEUCHOS_UNIT_TEST(Array, Deep_copyTest) {

	size_t size(10);
	Array<int> a(size),b;
	TEST_ASSERT(a.get_size()==size);
	TEST_ASSERT(b.get_size()==0);
	a.set(1);
	b.deep_copy(a);
	TEST_ASSERT(b.get_size()==size);
	for(std::size_t i=0;i<a.get_size();i++){
		TEST_ASSERT(1==b.get()[i]);
		TEST_ASSERT(b.get()[i]==a.get()[i]);
	}

}


TEUCHOS_UNIT_TEST(Array, ScaleTest) {

	size_t size(10);
	Array<int> a(size);
	TEST_ASSERT(a.get_size()==size);
	a.set(1);
	a.scale(10);
	for(std::size_t i=0;i<a.get_size();i++){
		TEST_ASSERT(10==a.get()[i]);
	}

}





int main
(
		int argc,
		char* argv[]
)
{

	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
