/*
 * utPdQuickGrid_solidCylinder_np2.cxx
 *
 *  Created on: May 18, 2010
 *      Author: jamitch
 */


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "PdVTK.h"
#include "vtkXMLStructuredGridWriter.h"
#include "vtkXMLStructuredGridReader.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"
#include "vtkSmartPointer.h"
#include "vtkStructuredGrid.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdGridData.h"
#include "PdutMpiFixture.h"
#include "PdNeighborhood.h"
#include "zoltan.h"
#include "PdZoltan.h"
#include "Field.h"
#include "mpi.h"
#include <tr1/memory>
#include <valarray>
#include <iostream>
#include <cmath>
#include <map>


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

using namespace Pdut;
using std::tr1::shared_ptr;
using std::cout;
using std::set;
using std::map;

static int myRank;
static int numProcs;
const double cylinderRadius = 1.0;
const int numRings = 2;
const double zStart = 0.0;
const double cylinderLength = 1.0;
const double horizon = .6;


TensorProductSolidCylinder getMeshGenerator(){

	TensorProductSolidCylinder meshGen(numProcs,cylinderRadius,numRings,zStart,cylinderLength);
	return meshGen;
}

void runTest() {
	TensorProductSolidCylinder meshGen = getMeshGenerator();
	BOOST_CHECK(57==meshGen.getNumGlobalCells());
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, meshGen);
	decomp = getLoadBalancedDiscretization(decomp);
	decomp = createAndAddNeighborhood(decomp,horizon);
	BOOST_CHECK(57==decomp.globalNumPoints);

	/*
	 * Write problem set up parameters to file
	 */
	int numPoints = decomp.numPoints;
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

	test_suite* proc = BOOST_TEST_SUITE( "utPdQuickGrid_solidCylinder_np2" );
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

	// Initialize MPI and timer
	PdutMpiFixture myMpi = PdutMpiFixture(argc,argv);

	// These are static (file scope) variables
	myRank = myMpi.rank;
	numProcs = myMpi.numProcs;
	/**
	 * This test only make sense for numProcs == 2
	 */
	if(2 != numProcs){
		std::cerr << "Unit test runtime ERROR: utPdQuickGrid_solidCylinder_np2 only makes sense on 2 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 2 ./utPdQuickGrid_solidCylinder_np2" << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
