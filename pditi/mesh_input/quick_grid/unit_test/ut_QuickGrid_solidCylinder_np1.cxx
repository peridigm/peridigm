/*
 * utPdQuickGrid_solidCylinder_np1.cxx
 *
 *  Created on: May 18, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "vtkXMLStructuredGridWriter.h"
#include "vtkXMLStructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "../QuickGrid.h"
#include "vtk/PdVTK.h"
#include "vtk/Field.h"
#include <tr1/memory>
#include <valarray>
#include <iostream>
#include <cmath>



using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;

static size_t myRank;
static size_t numProcs;
const double cylinderRadius = 1.0;
const size_t numRings = 2;
const double zStart = 0.0;
const double cylinderLength = 1.0;


QUICKGRID::TensorProductSolidCylinder getMeshGenerator(){

	QUICKGRID::TensorProductSolidCylinder meshGen(numProcs,cylinderRadius,numRings,zStart,cylinderLength);
	return meshGen;
}

void runTest() {
	QUICKGRID::TensorProductSolidCylinder meshGen = getMeshGenerator();
	BOOST_CHECK(57==meshGen.getNumGlobalCells());
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, meshGen);
	int numPoints = decomp.numPoints;
	BOOST_CHECK(57==numPoints);

	/*
	 * Assert that there are 3 core points
	 */
	double *x = decomp.myX.get();
	double dz = cylinderLength/3;
	double z0 = zStart + dz/2;
	const double tolerance = 1.0e-15;
	BOOST_CHECK_CLOSE(*(x+0),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x+1),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x+2),z0,tolerance);
	BOOST_CHECK_CLOSE(*(x+19*3+0),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x+19*3+1),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x+19*3+2),z0+dz,tolerance);
	BOOST_CHECK_CLOSE(*(x+2*19*3+0),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x+2*19*3+1),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x+2*19*3+2),z0+2.0*dz,tolerance);

	/*
	 * Check volume of core of 3 core points
	 */
	double rI = cylinderRadius/numRings/sqrt(M_PI);
	double vol = M_PI*rI*rI*dz;
	double *v = decomp.cellVolume.get();
	BOOST_CHECK_CLOSE(*(v+0),vol,tolerance);
	BOOST_CHECK_CLOSE(*(v+19),vol,tolerance);
	BOOST_CHECK_CLOSE(*(v+2*19),vol,tolerance);

	/*
	 * Check global assigned "global ids"
	 */
	int *ids = decomp.myGlobalIDs.get();
	for(int i=0;ids!=decomp.myGlobalIDs.get()+numPoints;ids++,i++){
		BOOST_CHECK(*ids==i);
	}

	/*
	 * Write problem set up parameters to file
	 */
	Field_NS::Field<double> volField = Field_NS::getVOLUME(decomp.cellVolume,numPoints);
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);
	PdVTK::writeField<double>(grid,volField);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("solidCylinderMesh.pvtu", numProcs, myRank);
	PdVTK::write(writer,grid);


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "ut_QuickGrid_solidCylinder_np1" );
	proc->add(BOOST_TEST_CASE( &runTest ));
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

	// These are static (file scope) variables
	myRank = 0;
	numProcs = 1;

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
