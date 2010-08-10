/*
 * utPdNeighborhoodSort.cxx
 *
 *  Created on: Mar 10, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "PdNeighborhood.h"
#include "PdQuickGrid.h"
#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <iterator>
#include <time.h>


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using namespace std;

const int NUMBER_OF_POINTS_SORTED = 100000;
const int N = NUMBER_OF_POINTS_SORTED;

void sortPoints(){
	/*
	 * Random coordinates between 0 and PI
	 */
	shared_ptr<double> xPtr(new double[3*N],Deleter<double>());
	{
		/*
		 * Initialize random number generator
		 */
		srand ( time(NULL) );
		double *X = xPtr.get();
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
	PdNeighborhood::CoordinateLabel labels[] = {PdNeighborhood::X,PdNeighborhood::Y,PdNeighborhood::Z};
	for(int label=0;label<3;label++){
		Coordinates c(xPtr,N);
		Coordinates::SortComparator compare = c.getSortComparator(labels[label]);
		std::tr1::shared_ptr<int> mapX = c.getIdentityMap();

		/*
		 * Note that sort iterator acts as comparator
		 */
		std::sort(mapX.get(),mapX.get()+N,compare);

		/*
		 * Now check that points in map are in order
		 */
		{
			double *X = xPtr.get();
			int *ids = mapX.get();
			int *end = ids+(N-1);
			for(;ids!=end;ids++){
				double x1 = X[3*(*ids)+label];
				double x2 = X[3*(*(ids+1))+label];
				BOOST_CHECK(x1<=x2);
			}
		}

	}


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPdNeighborhoodSort" );
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

