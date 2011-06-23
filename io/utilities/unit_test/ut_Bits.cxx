/*
 * ut_Bits.cxx
 *
 *  Created on: Apr 28, 2011
 *      Author: jamitch
 */

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
	unsigned int id = point | (gid<<8) | (length << 16);
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
