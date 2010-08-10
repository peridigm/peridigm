/*
 * utPeridigm_decompositionStates.cpp
 *
 *  Created on: Apr 28, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "Peridigm_DecompositionStates.hpp"
#include "PdutMpiFixture.h"
#include <iostream>

using namespace boost::unit_test;
using namespace Pdut;
using namespace std;
using namespace Peridigm;

static int myRank;
static int numProcs;

void addGetMethods(){

	Peridigm::DecompositionStates d;
	BOOST_CHECK(3==d.getNumScalarStateVariables());
	BOOST_CHECK(1==d.getNumVectorStateVariables());
	BOOST_CHECK(0==d.getNumScalarStateBondVariables());

	d.addScalarStateVariable("scalarWhatever");
	d.addVectorStateVariable("vectorWhatever");
	d.addScalarStateBondVariable("scalarBondVarWhatever");
	BOOST_CHECK(4==d.getNumScalarStateVariables());
	BOOST_CHECK(2==d.getNumVectorStateVariables());
	BOOST_CHECK(1==d.getNumScalarStateBondVariables());

	BOOST_CHECK("Weighted Volume"==d.getScalarStateName(0));
	BOOST_CHECK("Dilatation"==d.getScalarStateName(1));
	BOOST_CHECK("Damage"==d.getScalarStateName(2));
	BOOST_CHECK("scalarWhatever"==d.getScalarStateName(3));

	BOOST_CHECK("Current Position"==d.getVectorStateName(0));
	BOOST_CHECK("vectorWhatever"==d.getVectorStateName(1));
	BOOST_CHECK("scalarBondVarWhatever"==d.getScalarStateBondVarName(0));

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPeridigm_decompositionStates" );
	proc->add(BOOST_TEST_CASE( &addGetMethods ));
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

	// Initialize MPI and timer
	PdutMpiFixture myMpi = PdutMpiFixture(argc,argv);

	// These are static (file scope) variables
	myRank = myMpi.rank;
	numProcs = myMpi.numProcs;

	/**
	 * This test only make sense for numProcs == 1
	 */
	if(1 != numProcs){
		std::cerr << "Unit test runtime ERROR: utPeridigm_decompositionStates is intended for \"serial\" run only and makes sense on 1 processor" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 1 ./utPeridigm_decompositionStates" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
