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
#include <math.h>
#include <iostream>

using std::cout;
using std::endl;

using namespace boost::unit_test;




void bit_ids(){
	enum Relation {
				POINT=0,
				BOND,
				RELATION_UNDEFINED
	};

	enum Type {
		        VOLUME=0,
				GID,
				PROC_NUM,
				WEIGHTED_VOLUME,
				DILATATION,
				DAMAGE,
				E_DP,
				PLASTIC_CONSISTENCY,
				SHEAR_CORRECTION_FACTOR,
				NUM_NEIGHBORS,
				COORDINATES,
				DISPLACEMENT,
				CURRENT_COORDINATES,
				VELOCITY,
				ACCELERATION,
				BC_MASK,
				FORCE,
				FORCE_DENSITY,
				CONTACT_FORCE,
				CONTACT_FORCE_DENSITY,
				BOND_DAMAGE,
				TYPE_UNDEFINED
	};

	enum Length {
				SCALAR=1,
				VECTOR2D=2,
				VECTOR3D=3,
				LENGTH_UNDEFINED
	};

	Relation point = POINT;
	Type  gid = GID;
	Length length = VECTOR3D;

	cout << "point = " << point << "; GID = " << gid << "; length = " << length << endl;
	cout << "gid << 8 = " << (gid<<8) << endl;
	cout << "length << 16 = " << (length << 16) << endl;
	cout << "id = point | (gid<<8) | (length << 16) = " << (point | (gid<<8) | (length << 16)) << endl;

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "ut_Bits" );
	proc->add(BOOST_TEST_CASE( &bit_ids ));
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
