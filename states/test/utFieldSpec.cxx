/*
 * utFieldData.cxx
 *
 *  Created on: Oct 20, 2009
 *      Author: jamitch
 */

#define BOOST_TEST_MODULE
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "FieldSpec.h"

using namespace PdStates;
BOOST_AUTO_TEST_CASE( FieldSpecOperatorEqual )
{

	BOOST_CHECK ( PdStates::COORD1D == PdStates::COORD1D );
	BOOST_CHECK ( PdStates::COORD1D != PdStates::COORD2D );
	BOOST_CHECK ( PdStates::COORD1D != PdStates::COORD3D );
	BOOST_CHECK ( PdStates::COORD2D != PdStates::COORD1D );
	BOOST_CHECK ( PdStates::COORD2D == PdStates::COORD2D );
	BOOST_CHECK ( PdStates::COORD2D != PdStates::COORD3D );
	BOOST_CHECK ( PdStates::COORD3D != PdStates::COORD1D );
	BOOST_CHECK ( PdStates::COORD3D != PdStates::COORD2D );
	BOOST_CHECK ( PdStates::COORD3D == PdStates::COORD3D );

}

BOOST_AUTO_TEST_CASE( FieldSpecGetLength )
{

	BOOST_CHECK ( COORD1D.getLength() == FieldSpec::SCALAR  );
	BOOST_CHECK ( COORD2D.getLength() == FieldSpec::VECTOR2D );
	BOOST_CHECK ( COORD3D.getLength() == FieldSpec::VECTOR3D );

	BOOST_CHECK ( DISPL1D.getLength() == FieldSpec::SCALAR  );
	BOOST_CHECK ( DISPL2D.getLength() == FieldSpec::VECTOR2D );
	BOOST_CHECK ( DISPL3D.getLength() == FieldSpec::VECTOR3D );

	BOOST_CHECK ( VELOC1D.getLength() == FieldSpec::SCALAR  );
	BOOST_CHECK ( VELOC2D.getLength() == FieldSpec::VECTOR2D );
	BOOST_CHECK ( VELOC3D.getLength() == FieldSpec::VECTOR3D );

}

