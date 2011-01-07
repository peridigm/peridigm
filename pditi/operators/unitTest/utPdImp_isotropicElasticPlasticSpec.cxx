/*
 * utIsotropicElasticPlasticSpec.cxx
 *
 *  Created on: Jul 20, 2010
 *      Author: jamitch
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "../PdImpMpiFixture.h"
#include "../PdImpMaterials.h"
#include <utility>
#include <map>
#include <vector>
#include <cstdlib>
#include <time.h>

static int myRank;
static int numProcs;
/*
 * Young's Modulus (MPa)
 */
static double E = 68.9e3;

/*
 * Poisson's ratio
 */
static double nu = .33;

/*
 * Yield "Stress" estimate for perfect plasticity (MPa)
 * 6061-T6 data
 */
static double Y = 351.79;

/*
 * Density of aluminum g/mm^3
 */
static double rho = 2.7e-3;

/*
 * Bulk Modulus
 */
static double K = E / 3 / (1-2.0 * nu);

/*
 * Shear Modulus
 */
static double mu = E / 2 / (1+nu);


using namespace boost::unit_test;
using namespace PdImp;

void assertIsotropicElasticPlasticSpec() {
	/*
	 * Yield stress must be > 0
	 */
	BOOST_CHECK_THROW(IsotropicElasticPlasticSpec::yieldStress(0.0),std::domain_error);
	BOOST_CHECK_THROW(IsotropicElasticPlasticSpec::yieldStress(-1.0),std::domain_error);
	BOOST_CHECK_NO_THROW(IsotropicElasticPlasticSpec::yieldStress(Y));
	YieldStress yield = IsotropicElasticPlasticSpec::yieldStress(Y);
	/*
	 * Horizon must be > 0
	 */
	BOOST_CHECK_THROW(IsotropicElasticPlasticSpec::materialHorizon(-1.0),std::domain_error);
	BOOST_CHECK_THROW(IsotropicElasticPlasticSpec::materialHorizon(0.0),std::domain_error);
	BOOST_CHECK_NO_THROW(IsotropicElasticPlasticSpec::materialHorizon(.49));
	MaterialHorizon horizon = IsotropicElasticPlasticSpec::materialHorizon(.49);

	YoungsModulus youngsModulus = IsotropicHookeSpec::youngsModulus(E);
	PoissonsRatio poissonsRatio = IsotropicHookeSpec::poissonsRatio(nu);
	IsotropicElasticPlasticSpec elasticPlasticSpec(yield,horizon,IsotropicHookeSpec(youngsModulus,poissonsRatio));
	IsotropicHookeSpec hookeSpec = elasticPlasticSpec.getHookeSpec();
	double tolerance(1.0e-15);
	BOOST_CHECK_CLOSE(E,hookeSpec.youngsModulus(),tolerance);
	BOOST_CHECK_CLOSE(nu,hookeSpec.poissonsRatio(),tolerance);
	BOOST_CHECK_CLOSE(Y,yield.getValue(),tolerance);
	BOOST_CHECK_CLOSE(.49,horizon.getValue(),tolerance);


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPdImp_isotropicElasticPlasticSpec" );
	proc->add(BOOST_TEST_CASE( &assertIsotropicElasticPlasticSpec ));
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
	PdImpRunTime::PimpMpiFixture pimpMPI = PdImpRunTime::PimpMpiFixture::getPimpMPI(argc,argv);
	const Epetra_Comm& comm = pimpMPI.getEpetra_Comm();

	// These are static (file scope) variables
	myRank = comm.MyPID();
	numProcs = comm.NumProc();

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
