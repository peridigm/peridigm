/*! \file utPdNeighborhoodSort.cxx */

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

