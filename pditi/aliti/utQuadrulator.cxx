/*
 * utQuadrulator.cxx
 *
 *  Created on: Nov 13, 2010
 *      Author: awesome
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "quadrulator.h"
#include <cstdio>
#include <iostream>
#include "PdVTK.h"
#include <vtkCellType.h>

using std::size_t;
using namespace boost::unit_test;

static int numProcs;
static int myRank;

void getGrid(){
	double x0 = 0.0, y0 = 0.0;
	size_t nx = 3, ny = 3;
	size_t hx = 1.0, hy = 1.0;
	Quadrulator reader (x0, y0, nx, hx, ny, hy);

	ArrayRCP<size_t> links = reader.getVertexLinks();
	BOOST_CHECK(16 == reader.getNumVertices());
	BOOST_CHECK(4  == reader.getNpe());
	BOOST_CHECK(9  == reader.getNem());
	size_t nem=reader.getNem();
	size_t numVertices = reader.getNumVertices();
	size_t npe = reader.getNpe();

	// Check vertex connectivity
	size_t l[] =
	{
			0,1,5,4,
			1,2,6,5,
			2,3,7,6,
			4,5,9,8,
			5,6,10,9,
			6,7,11,10,
			8,9,13,12,
			9,10,14,13,
			10,11,15,14
	};
	for(size_t e=0;e<nem;e++){
		for(size_t n=0;n<npe;n++){
			size_t left = links[npe*e+n];
			size_t right = l[npe*e+n];
			BOOST_CHECK(left==right);
		}
	}

	// now assert the local element coordinates; these are coordinates
	//   associated with each element based on the nodal connectivity
	ArrayRCP<double> eX = reader.getElementCoordinates();
	int el = 0;
	BOOST_CHECK(0.0 == eX[12*el]);
	BOOST_CHECK(0.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(1.0 == eX[12*el+3]);
	BOOST_CHECK(0.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(1.0 == eX[12*el+6]);
	BOOST_CHECK(1.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(0.0 == eX[12*el+9]);
	BOOST_CHECK(1.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 1;
	BOOST_CHECK(1.0 == eX[12*el]);
	BOOST_CHECK(0.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(2.0 == eX[12*el+3]);
	BOOST_CHECK(0.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(2.0 == eX[12*el+6]);
	BOOST_CHECK(1.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(1.0 == eX[12*el+9]);
	BOOST_CHECK(1.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 2;
	BOOST_CHECK(2.0 == eX[12*el]);
	BOOST_CHECK(0.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(3.0 == eX[12*el+3]);
	BOOST_CHECK(0.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(3.0 == eX[12*el+6]);
	BOOST_CHECK(1.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(2.0 == eX[12*el+9]);
	BOOST_CHECK(1.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 3;
	BOOST_CHECK(0.0 == eX[12*el]);
	BOOST_CHECK(1.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(1.0 == eX[12*el+3]);
	BOOST_CHECK(1.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(1.0 == eX[12*el+6]);
	BOOST_CHECK(2.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(0.0 == eX[12*el+9]);
	BOOST_CHECK(2.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 4;
	BOOST_CHECK(1.0 == eX[12*el]);
	BOOST_CHECK(1.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(2.0 == eX[12*el+3]);
	BOOST_CHECK(1.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(2.0 == eX[12*el+6]);
	BOOST_CHECK(2.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(1.0 == eX[12*el+9]);
	BOOST_CHECK(2.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 5;
	BOOST_CHECK(2.0 == eX[12*el]);
	BOOST_CHECK(1.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(3.0 == eX[12*el+3]);
	BOOST_CHECK(1.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(3.0 == eX[12*el+6]);
	BOOST_CHECK(2.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(2.0 == eX[12*el+9]);
	BOOST_CHECK(2.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 6;
	BOOST_CHECK(0.0 == eX[12*el]);
	BOOST_CHECK(2.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(1.0 == eX[12*el+3]);
	BOOST_CHECK(2.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(1.0 == eX[12*el+6]);
	BOOST_CHECK(3.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(0.0 == eX[12*el+9]);
	BOOST_CHECK(3.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 7;
	BOOST_CHECK(1.0 == eX[12*el]);
	BOOST_CHECK(2.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(2.0 == eX[12*el+3]);
	BOOST_CHECK(2.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(2.0 == eX[12*el+6]);
	BOOST_CHECK(3.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(1.0 == eX[12*el+9]);
	BOOST_CHECK(3.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);
	el = 8;
	BOOST_CHECK(2.0 == eX[12*el]);
	BOOST_CHECK(2.0 == eX[12*el+1]);
	BOOST_CHECK(0.0 == eX[12*el+2]);
	BOOST_CHECK(3.0 == eX[12*el+3]);
	BOOST_CHECK(2.0 == eX[12*el+4]);
	BOOST_CHECK(0.0 == eX[12*el+5]);
	BOOST_CHECK(3.0 == eX[12*el+6]);
	BOOST_CHECK(3.0 == eX[12*el+7]);
	BOOST_CHECK(0.0 == eX[12*el+8]);
	BOOST_CHECK(2.0 == eX[12*el+9]);
	BOOST_CHECK(3.0 == eX[12*el+10]);
	BOOST_CHECK(0.0 == eX[12*el+11]);

	vtkSmartPointer<vtkPoints> vtkX = PdVTK::createVTK_Points(reader.getCoordinates().get(),reader.getNumVertices());
	vtkSmartPointer<vtkCellArray> vtkQuadCells = PdVTK::createVTK_quadCells(links.get(),nem);
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(vtkX,vtkQuadCells,VTK_QUAD);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utQuadrulator.pvtu",numProcs,myRank,PdVTK::vtkASCII);
	PdVTK::write(writer,grid);
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utQuadrulator" );
	proc->add(BOOST_TEST_CASE( &getGrid ));
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

	numProcs = 1;
	myRank = 0;

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
