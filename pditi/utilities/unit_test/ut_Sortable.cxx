/*
 * ut_Sortable.cxx
 *
 *  Created on: Mar 12, 2011
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <math.h>
#include "Sortable.h"

using UTILITIES::Array;
using UTILITIES::Sortable;
using UTILITIES::CartesianComponent;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

const int NUMBER_OF_POINTS_SORTED = 10000;
const int N = NUMBER_OF_POINTS_SORTED;



void sortPoints(){
	/*
	 * Random coordinates between 0 and PI
	 */
	Array<double> coordinates(3*N);
	{
		/*
		 * Initialize random number generator
		 */
		srand ( time(NULL) );
		double *X = coordinates.get();
		double pi = M_PI;
		for(int p=0;p<N;p++,X+=3){
			*(X+0)= pi*(rand()%N)/N;
			*(X+1)= pi*(rand()%N)/N;
			*(X+2)= pi*(rand()%N)/N;
		}
	}
	/*
	 * sort points
	 */
	int numAxes=3;
	Array<CartesianComponent> components(numAxes);
	components[0] = UTILITIES::X;
	components[1] = UTILITIES::Y;
	components[2] = UTILITIES::Z;

	Sortable points(N, coordinates.get_shared_ptr());
	Array< Array<int> > sorted_maps(numAxes);
	for(size_t c=0;c<components.get_size();c++){

		Sortable::Comparator compare = points.getComparator(components[c]);
		Array<int> mapX = points.getIdentityMap();
		sorted_maps[c] = mapX;

		/*
		 * Note that sort iterator acts as comparator
		 */
		std::sort(mapX.get(),mapX.get()+N,compare);

	}

	for(size_t c=0;c<components.get_size();c++){

		/*
		 * Now check that points in map are in order
		 */
		{
			double *X = coordinates.get();
			int *ids = sorted_maps[c].get();
			int *end = ids+(N-1);
			for(;ids!=end;ids++){
				double x1 = X[3*(*ids)+c];
				double x2 = X[3*(*(ids+1))+c];
				BOOST_CHECK(x1<=x2);
			}
		}

	}


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "ut_Sortable" );
	proc->add(BOOST_TEST_CASE( &sortPoints ));
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
