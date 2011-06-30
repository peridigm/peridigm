/*
 * ut_twoPointReloadBalance.cxx
 *
 *  Created on: Jun 29, 2011
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
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include <tr1/memory>
#include <iostream>
#include <cmath>
#include <map>
#include <set>

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "Array.h"
#include "../PdZoltan.h"
#include "../BondFilter.h"
#include "../NeighborhoodList.h"
#include "quick_grid/QuickGrid.h"

#include "vtk/Field.h"
#include "PdutMpiFixture.h"


using std::tr1::shared_ptr;
using namespace boost::unit_test;

using namespace Pdut;
using std::cout;
using std::set;
using std::map;

static size_t myRank;
static size_t numProcs;
const size_t nx = 2;
const size_t ny = 1;
const size_t nz = 1;
const double xStart =  -2.0;
const double yStart =  -0.5;
const double zStart =  -0.5;
const double xLength = 4.0;
const double yLength = 1.0;
const double zLength = 1.0;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;
const double horizon = 2.1;

QUICKGRID::QuickGridData getGrid() {
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,QUICKGRID::SphericalNorm);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	QUICKGRID::print_meta_data(decomp);
	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}

void twoPointReloadBalance() {
	QUICKGRID::QuickGridData decomp_1 = getGrid();

	/*
	 * Reload balance
	 */
	QUICKGRID::QuickGridData decomp_2 = PDNEIGH::getLoadBalancedDiscretization(decomp_1);

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "ut_twoPointReloadBalance" );
	proc->add(BOOST_TEST_CASE( &twoPointReloadBalance ));
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

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}

