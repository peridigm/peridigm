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
#include "Teuchos_UnitTestRepository.hpp"
#include <cstddef>
#include <cstdlib>
#include <vector>
#include "kdtree.h"
#include <iostream>
#include <cmath>

using namespace::std;

typedef double value_type;
typedef std::size_t ordinal_type;
typedef femanica::kdtree<value_type,ordinal_type> my_tree;

TEUCHOS_UNIT_TEST(KdTree_spherical_search, No_PointsTest) {

	ordinal_type num_points=0;
	vector<value_type> points(3*num_points);
	TEST_ASSERT(points.size()==num_points);
	my_tree tree=my_tree::get_tree(points.data(),num_points);

	/*
	 * search
	 */
	std::vector<ordinal_type> neighbors;
	value_type search_point[]={0.0,0.0,0.0};
	value_type R(1.0);
	tree.all_neighbors_within_radius(search_point,R,neighbors);
	TEST_ASSERT(neighbors.size()==0);

}

TEUCHOS_UNIT_TEST(KdTree_spherical_search, One_PointTest) {

	ordinal_type num_points=1;
	vector<value_type> points(3);
    points[0] = 0.0 ; points[1] = 0.0 ; points[2] = 0.0;
	my_tree tree=my_tree::get_tree(points.data(),num_points);
	/*
	 * search
	 */
	std::vector<ordinal_type> neighbors;
	value_type search_point[]={0.0,0.0,0.0};
	value_type R(1.0);
	tree.all_neighbors_within_radius(search_point,R,neighbors);
	TEST_ASSERT(neighbors.size()==1);
	TEST_ASSERT(neighbors[0]==0);

}

TEUCHOS_UNIT_TEST(KdTree_spherical_search, Two_PointsTest) {

	ordinal_type num_points=2;
	vector<value_type> points(6);
    points[0] = 2.0 ; points[1] = 0.0 ; points[2] = 0.0;
    points[3] = 3.0 ; points[4] = 0.0 ; points[5] = 0.0;
	my_tree tree=my_tree::get_tree(points.data(),num_points);

	{
		/*
		 * search point not in tree and not close enough to tree
		 */
		std::vector<ordinal_type> neighbors;
		value_type search_point[]={0.0,0.0,0.0};
		value_type R(0.25);
		tree.all_neighbors_within_radius(search_point,R,neighbors);
		TEST_ASSERT(neighbors.size()==0);

	}

	{
		/*
		 * search point not in tree and but with horizon that finds first point
		 */
		std::vector<ordinal_type> neighbors;
		value_type search_point[]={0.0,0.0,0.0};
		value_type R(2.0);
		tree.all_neighbors_within_radius(search_point,R,neighbors);
		TEST_ASSERT(neighbors.size()==1);
		TEST_ASSERT(neighbors[0]==0);
	}

	{
		/*
		 * search point not in tree and but with horizon that finds both points
		 */
		std::vector<ordinal_type> neighbors;
		value_type search_point[]={0.0,0.0,0.0};
		value_type R(3.0);
		tree.all_neighbors_within_radius(search_point,R,neighbors);
		TEST_ASSERT(neighbors.size()==2);
		TEST_ASSERT(neighbors[0]==0 || neighbors[0]==1);
		TEST_ASSERT(neighbors[1]==0 || neighbors[1]==1);
	}

	{
		/*
		 * search first point
		 */
		std::vector<ordinal_type> neighbors;
		value_type search_point[]={2.0,0.0,0.0};
		value_type R(0.25);
		tree.all_neighbors_within_radius(search_point,R,neighbors);
		TEST_ASSERT(neighbors.size()==1);
		TEST_ASSERT(neighbors[0]==0);
	}
	{
		/*
		 * search first point; change radius
		 */
		std::vector<ordinal_type> neighbors;
		value_type search_point[]={2.0,0.0,0.0};
		value_type R(0.5);
		tree.all_neighbors_within_radius(search_point,R,neighbors);
		TEST_ASSERT(neighbors.size()==1);
		TEST_ASSERT(neighbors[0]==0);
	}
	{
		/*
		 * search first point; change radius
		 */
		std::vector<ordinal_type> neighbors;
		value_type search_point[]={2.0,0.0,0.0};
		value_type R(1.0);
		tree.all_neighbors_within_radius(search_point,R,neighbors);
		TEST_ASSERT(neighbors.size()==2);
		TEST_ASSERT(neighbors[0]==0 || neighbors[0]==1);
		TEST_ASSERT(neighbors[1]==0 || neighbors[1]==1);
		TEST_ASSERT(neighbors[0]==0 || neighbors[1]==0);
		TEST_ASSERT(neighbors[0]==1 || neighbors[1]==1);

	}
}

TEUCHOS_UNIT_TEST(KdTree_spherical_search, Three_PointsTest) {

	{
        size_t num_cells(2);
        const size_t dimension(3);
        vector<value_type> points;
        points.push_back(-0.5) ; points.push_back(0.0) ;  points.push_back(0.0) ;
        points.push_back(0.5)  ; points.push_back(0.0) ;  points.push_back(0.0) ;

		TEST_ASSERT(num_cells*dimension==points.size());

		// add point at origin
		points.push_back(0.0); points.push_back(0.0); points.push_back(0.0);
		num_cells+=1;
		TEST_ASSERT(num_cells*dimension==points.size());

		my_tree tree=my_tree::get_tree(points.data(),num_cells);

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.25,0.0,0.0};
			value_type R(0.2);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * find origin
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(0.25);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==1);
			/*
			 * origin was last point added
			 */
			TEST_ASSERT(2==neighbors[0]);
		}

		{
			/*
			 * find origin and point on right
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.25,0.0,0.0};
			value_type R(0.26);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==2);
			/*
			 * origin was last point added; right point is 2nd point
			 */
			TEST_ASSERT(2==neighbors[0] || 1==neighbors[0]);
			TEST_ASSERT(2==neighbors[1] || 1==neighbors[1]);
		}

		{
			/*
			 * find origin and point on left
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={-0.25,0.0,0.0};
			value_type R(0.26);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==2);
			/*
			 * origin was last point added; left point is first point
			 */
			TEST_ASSERT(2==neighbors[0] || 0==neighbors[0]);
			TEST_ASSERT(2==neighbors[1] || 0==neighbors[1]);
		}

		{
			/*
			 * find all three points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(0.51);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==3);
			/*
			 * origin was last point added; left point is first point
			 */
		        TEST_ASSERT(2==neighbors[0] || 0==neighbors[0] || 1==neighbors[0]);
			TEST_ASSERT(2==neighbors[1] || 0==neighbors[1] || 1==neighbors[1]);
			TEST_ASSERT(2==neighbors[2] || 0==neighbors[2] || 1==neighbors[2]);
		}

	}

	{
        size_t num_cells(2);
        const size_t dimension(3);
        vector<value_type> points;
        points.push_back(-0.5) ; points.push_back(0.0) ;  points.push_back(0.0) ;
        points.push_back(0.5)  ; points.push_back(0.0) ;  points.push_back(0.0) ;

        TEST_ASSERT(num_cells*dimension==points.size());

		// add point at 0,.5,0
		points.push_back(0.0); points.push_back(0.5); points.push_back(0.0);
		num_cells+=1;
		TEST_ASSERT(num_cells*dimension==points.size());

		my_tree tree=my_tree::get_tree(points.data(),num_cells);

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(0.2);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={1.0,0.0,0.0};
			value_type R(0.2);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,1.0,0.0};
			value_type R(0.2);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,1.0};
			value_type R(0.2);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * try finding all points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(0.51);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==3);
		}


	}


}

TEUCHOS_UNIT_TEST(KdTree_spherical_search, Four_PointsTest) {

	{
        size_t num_cells(4);
        const size_t dimension(3);
        double x_length = 2.0;
        vector<value_type> points;
        points.push_back(-0.5) ; points.push_back(-0.5) ;  points.push_back(0.0) ;
        points.push_back(0.5)  ; points.push_back(-0.5) ;  points.push_back(0.0) ;
        points.push_back(-0.5) ; points.push_back(0.5)  ;  points.push_back(0.0) ;
        points.push_back(0.5)  ; points.push_back(0.5)  ;  points.push_back(0.0) ;

        TEST_ASSERT(num_cells*dimension==points.size());

		my_tree tree=my_tree::get_tree(points.data(),num_cells);

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(0.2);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * try finding all points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(x_length/sqrt(2.0));
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==4);
		}

		{
			/*
			 * try finding three points using lower left
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={-0.5,-0.5,0.0};
			value_type R(1.0);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==3);
		}
		{
			/*
			 * try finding three points using lower right
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.5,-0.5,0.0};
			value_type R(1.0);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==3);
		}

		{
			/*
			 * try finding three points using upper right
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.5,0.5,0.0};
			value_type R(1.0);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==3);
		}
		{
			/*
			 * try finding three points using upper left
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={-0.5,0.5,0.0};
			value_type R(1.0);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==3);
		}

	}
}

TEUCHOS_UNIT_TEST(KdTree_spherical_search, n_xn_PointsTest) {

	{
		/*
		 * this routine will work for arbitrary n that are 'even'
		 */
        size_t n = 1000;
        size_t nx = n;
        size_t ny = n;
        size_t nz = 1;
        size_t num_cells(nx*ny*nz);
        const size_t dimension(3);
        double x_length = 2.0;
        double cell_size = x_length/n;
        vector<value_type> points;
        for(int i=0 ; i<static_cast<int>(n) ; ++i){
          for(int j=0 ; j<static_cast<int>(n) ; ++j){
            points.push_back(-1.0 + cell_size/2.0 + i*cell_size);
            points.push_back(-1.0 + cell_size/2.0 + j*cell_size);
            points.push_back(0.0);
          }
        }

        TEST_ASSERT(num_cells*dimension==points.size());

		my_tree tree=my_tree::get_tree(points.data(),num_cells);

		{
			/*
			 * try finding no points
			 */
			std::vector<ordinal_type> neighbors;
			value_type search_point[]={0.0,0.0,0.0};
			value_type R(.9*cell_size/2.0);
			tree.all_neighbors_within_radius(search_point,R,neighbors);
			TEST_ASSERT(neighbors.size()==0);
		}

		{
			/*
			 * incrementally increase search radius
			 */
			size_t num_neighbors[]={0,4,24,52,96};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={0.0,0.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}

		{
			/*
			 * search center right
			 */
			size_t num_neighbors[]={0,4/2,24/2,52/2,96/2};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={1.0,0.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}

		{
			/*
			 * search top center
			 */
			size_t num_neighbors[]={0,4/2,24/2,52/2,96/2};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={0.0,1.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}

		{
			/*
			 * search center left
			 */
			size_t num_neighbors[]={0,4/2,24/2,52/2,96/2};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={-1.0,0.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}


		{
			/*
			 * search center bottom
			 */
			size_t num_neighbors[]={0,4/2,24/2,52/2,96/2};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={0.0,-1.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}

		{
			/*
			 * search upper right
			 */
			size_t num_neighbors[]={0,4/4,24/4,52/4,96/4};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={1.0,1.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}


		{
			/*
			 * search upper left
			 */
			size_t num_neighbors[]={0,4/4,24/4,52/4,96/4};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={-1.0,1.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}

		{
			/*
			 * search lower left
			 */
			size_t num_neighbors[]={0,4/4,24/4,52/4,96/4};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={-1.0,-1.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}

		{
			/*
			 * search lower right
			 */
			size_t num_neighbors[]={0,4/4,24/4,52/4,96/4};
			std::vector<ordinal_type> neighbors;
			for(int i=0;i<5;i++){
				value_type search_point[]={1.0,-1.0,0.0};
				value_type R(sqrt(2.0)*i*cell_size);
				tree.all_neighbors_within_radius(search_point,R,neighbors);
				TEST_ASSERT(neighbors.size()==num_neighbors[i]);
			}

		}


	}

}


int main(int argc, char* argv[]) {

	// Initialize UTF
	return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}


