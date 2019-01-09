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

#include <math.h>
#include "Sortable.h"

using UTILITIES::Array;
using UTILITIES::Sortable;
using UTILITIES::CartesianComponent;
using std::shared_ptr;


const int NUMBER_OF_POINTS_SORTED = 10000;
const int N = NUMBER_OF_POINTS_SORTED;


TEUCHOS_UNIT_TEST(Array, SortPointsTest) {
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
				TEST_ASSERT(x1<=x2);
			}
		}

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
