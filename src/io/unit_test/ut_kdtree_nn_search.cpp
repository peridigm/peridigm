/*! \file ut_kdtree_nn_search.cpp */


#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cwchar>
#include <vector>
#include "kdtree.h"
#include "../quick_grid/QuickGrid.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ostream>
#include <Epetra_SerialComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"


using std::size_t;
using std::vector;

using QUICKGRID::Spec1D;
using femanica::kdtree;

typedef double value_type;
typedef unsigned int ordinal_type;
typedef kdtree<value_type,ordinal_type> my_tree;

void fill_vector_interlaced
	(
			ordinal_type num_points,
			value_type* x,
			value_type* y,
			value_type* z,
			vector<value_type>& fill_me
	){
	fill_me.clear();
	for(ordinal_type n=0;n<num_points;n++){
		fill_me.push_back(x[n]);
		fill_me.push_back(y[n]);
		fill_me.push_back(z[n]);
	}
}


value_type get_random_point(value_type center, value_type width){
	value_type tolerance=1.0e-15;
	value_type eps=width*tolerance;
	value_type min=center-width/2.0+eps;
	value_type max=center+width/2.0-eps;
	return ((value_type)(rand())/((value_type)(RAND_MAX+1.0)))*(max-min)+min;
}



TEUCHOS_UNIT_TEST(KdTree_nn_Search, No_Points) {

	ordinal_type num_points=0;
	value_type *x=0;
	value_type *y=0;
	value_type *z=0;
	vector<value_type> points(3*num_points);
	fill_vector_interlaced(num_points,x,y,z,points);
	TEST_ASSERT(points.size()==num_points);
	my_tree tree=my_tree::get_tree(points.data(),num_points);

	/*
	 * Nearest neighbor search
	 */
	value_type search_point[]={2.0,2.0,0.0};
	vector<value_type> P(search_point,search_point+3);
	TEST_THROW(tree.nearest_neighbor_search(P),std::invalid_argument);

}


TEUCHOS_UNIT_TEST(KdTree_nn_Search, One_Point) {
	ordinal_type num_points=1;
	value_type x[]={0.0};
	value_type y[]={0.0};
	value_type z[]={0.0};
	vector<value_type> points(3*num_points);
	fill_vector_interlaced(num_points,x,y,z,points);
	my_tree tree=my_tree::get_tree(points.data(),num_points);
	/*
	 * Nearest neighbor search
	 */
	value_type search_point[]={2.0,2.0,0.0};
	vector<value_type> P(search_point,search_point+3);
	ordinal_type best=tree.nearest_neighbor_search(P);
	TEST_ASSERT(0==best);

}

TEUCHOS_UNIT_TEST(KdTree_nn_Search, Two_Points){

	ordinal_type num_points=2;
	value_type x[]={2.0,5.0};
	value_type y[]={3.0,4.0};
	value_type z[]={0.0,0.0};
	vector<value_type> points(3*num_points);
	fill_vector_interlaced(num_points,x,y,z,points);
	my_tree tree=my_tree::get_tree(points.data(),num_points);

	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,2.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={5.0,5.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(1==best);
	}

}

TEUCHOS_UNIT_TEST(KdTree_nn_Search, Three_Points){

	ordinal_type num_points=3;
	value_type x[]={2.0,5.0,9.0};
	value_type y[]={3.0,4.0,6.0};
	value_type z[]={0.0,0.0,0.0};
	vector<value_type> points(3*num_points);
	fill_vector_interlaced(num_points,x,y,z,points);
	my_tree tree=my_tree::get_tree(points.data(),num_points);

	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,2.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{

		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={5.0,5.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(1==best);
	}
	{

		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,3.1,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{

		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,4.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={3.6,4.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(1==best);
	}

	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={8.0,7.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(2==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={10.0,10.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(2==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={0.0,0.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
}


TEUCHOS_UNIT_TEST(KdTree_nn_Search, Four_Points){

	ordinal_type num_points=4;
	value_type x[]={2.0,5.0,9.0,3.0};
	value_type y[]={3.0,4.0,6.0,9.0};
	value_type z[]={0.0,0.0,0.0,0.0};
	vector<value_type> points(3*num_points);
	fill_vector_interlaced(num_points,x,y,z,points);
	my_tree tree=my_tree::get_tree(points.data(),num_points);

	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,2.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{

		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={5.0,5.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(1==best);
	}
	{

		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,3.1,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{

		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={2.0,4.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={3.6,4.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(1==best);
	}

	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={8.0,7.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(2==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={10.0,10.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(2==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={0.0,0.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(0==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={3.0,8.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(3==best);
	}
	{
		/*
		 * Nearest neighbor search
		 */
		value_type search_point[]={4.0,9.0,0.0};
		vector<value_type> P(search_point,search_point+3);
		ordinal_type best=tree.nearest_neighbor_search(P);
		TEST_ASSERT(3==best);
	}
	{

		/* initialize random seed: */
		srand (time(NULL));
		value_type x0=4.0,y0=9.0,z0=0.0;
		value_type wx=2.0,wy=2.0,wz=2.0;

		ordinal_type num_points=100;
		for(ordinal_type p=0;p<num_points;p++){
			value_type x=get_random_point(x0,wx);
			value_type y=get_random_point(y0,wy);
			value_type z=get_random_point(z0,wz);
			value_type search_point[]={x,y,z};
			vector<value_type> P(search_point,search_point+3);
			ordinal_type best=tree.nearest_neighbor_search(P);
			//std::cout << "x,y,z" << std::scientific << std::setw(15) << x << ", " << y << ", " << z << std::endl;
			TEST_ASSERT(3==best);
		}

	}
}

vector<value_type> getDiscretization(const Spec1D& spec){
	size_t numCells = spec.getNumCells();
	vector<double> x(numCells);
	double x0=spec.getX0();
	double cellSize=spec.getCellSize();
	double p = x0+cellSize/2.0;
	for(size_t i=0;i<numCells;p+=cellSize,i++)
		x[i]=p;
	return x;
}

vector<value_type> get_discretization(const Spec1D& xSpec, const Spec1D& ySpec, const Spec1D& zSpec){
	// Set points and cells
	// note number of points is same as number of cells
	size_t nx = xSpec.getNumCells();
	size_t ny = ySpec.getNumCells();
	size_t nz = zSpec.getNumCells();
	size_t numCells = nx*ny*nz;

	vector<value_type> x = ::getDiscretization(xSpec);
	vector<value_type> y = ::getDiscretization(ySpec);
	vector<value_type> z = ::getDiscretization(zSpec);

	size_t dimension=3;
	vector<value_type> X(numCells*dimension);
	ordinal_type c=0;
	value_type point[3]={0.0,0.0,0.0};
	for(size_t k=0;k<nz;k++){
		point[2]=z[k];
		for(size_t j=0;j<ny;j++){
			point[1]=y[j];
			for(size_t i=0;i<nx;i++){
				point[0]=x[i];
				for(size_t p=0;p<3;p++,c++)
					X[c] = point[p];
			}
		}
	}

	return X;
}


TEUCHOS_UNIT_TEST(KdTree_nn_Search, N_X_N_Points){
	/*
	 * this routine will work for arbitrary n that are 'even'
	 */
	size_t n=300;
	size_t nx=n,ny=n,nz=1;
	value_type x0=-1.0;
	value_type x_length=2.0;
	Spec1D x_spec(nx,x0,x_length);
	Spec1D y_spec(ny,x0,x_length);
	Spec1D z_spec(nz,x0,x_length);
	const size_t dimension(3);
	size_t num_cells=x_spec.getNumCells()*z_spec.getNumCells()*y_spec.getNumCells();
	vector<value_type> points=get_discretization(x_spec,y_spec,z_spec);
	TEST_ASSERT(num_cells*dimension==points.size());
	value_type cell_size=x_spec.getCellSize();

	TEST_ASSERT(num_cells*dimension==points.size());
	std::cout << "START kdtree construction with " << num_cells << " points." << std::endl;
	my_tree tree=my_tree::get_tree(points.data(),num_cells);
	std::cout << "\tFINISHED kdtree construction." << std::endl;

	/* initialize random seed: */
	srand (time(NULL));

	ordinal_type num_search_points=1;
	std::cout << "Start nearest neighbor search of " << num_search_points*num_cells << " randomly positioned points." << std::endl;
	std::cout << "\tRandom points are positioned near each point in kdtree so that answer is known apriori." << std::endl;
	// search points randomly generated about each point
	for(size_t i=0;i<num_cells;i++){
		value_type x0=points[3*i],y0=points[3*i+1],z0=points[3*i+2];
		value_type wx=cell_size,wy=cell_size,wz=cell_size;
		//std::cout << "search 100 random points about\n";
		//std::cout << "x0,y0,z0" << std::scientific << std::setw(15) << x0 << ", " << y0 << ", " << z0  << std::endl;
		//std::cout << "wx,wy,wz" << std::scientific << std::setw(15) << wx << ", " << wy << ", " << wz << std::endl;
		for(ordinal_type p=0;p<num_search_points;p++){
			value_type x=get_random_point(x0,wx);
			value_type y=get_random_point(y0,wy);
			value_type z=get_random_point(z0,wz);
			value_type search_point[]={x,y,z};
			//std::cout << "\tRandom point: x,y,z" << std::scientific << std::setw(15) << x << ", " << y << ", " << z  << std::endl;
			vector<value_type> P(search_point,search_point+3);
			ordinal_type best=tree.nearest_neighbor_search(P);
			TEST_ASSERT(best==i);
		}
	}
}


int main
(int argc, char* argv[])
{
  
    // Run the tests
    
    
    return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
 
}
