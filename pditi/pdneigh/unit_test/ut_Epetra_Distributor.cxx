/*
 * ut_Epetra_Distributor.cxx
 *
 *  Created on: Mar 19, 2011
 *      Author: wow
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
#include "../PdZoltan.h"
#include "quick_grid/QuickGrid.h"
#include "../NeighborhoodList.h"
#include "../BondFilter.h"
#include "../OverlapDistributor.h"
#include "vtk/PdVTK.h"
#include "PdutMpiFixture.h"
#include "zoltan.h"
#include "vtk/Field.h"
#include "Array.h"
#include "mpi.h"
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



using namespace Field_NS;
using namespace PdBondFilter;
using namespace PDNEIGH;
using UTILITIES::Array;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

using namespace Pdut;
using std::tr1::shared_ptr;
using std::cout;
using std::endl;
using std::set;
using std::map;


using namespace QUICKGRID;
using namespace PDNEIGH;
using namespace Field_NS;

static size_t myRank;
static size_t numProcs;
const int nx = 4;
const int ny = nx;
const int nz = 1;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = xLength;
const double zStart = 0.0;
const double zLength = xLength/nx;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;
const double horizon = 1.1 * std::sqrt(xSpec.getCellSize()*xSpec.getCellSize()+ySpec.getCellSize()*ySpec.getCellSize());

void printNeighborhood(int numNeigh, int* neigh){
	for(int i=0;i<numNeigh;i++,neigh++){
		cout << ", " << *neigh;
	}
	cout << endl;
}

QUICKGRID::QuickGridData getGrid() {
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);

	// This load-balances
	decomp = getLoadBalancedDiscretization(decomp);
	return decomp;
}

PDNEIGH::NeighborhoodList createNeighborhood(QUICKGRID::QuickGridData decomp) {
	shared_ptr<BondFilter> bondFilterPtr(new PdBondFilter::BondFilterDefault(true));
	shared_ptr<Epetra_Comm> comm(new Epetra_MpiComm(MPI_COMM_WORLD));
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon,bondFilterPtr);
	return list;
}

void createDistributor() {
	QUICKGRID::QuickGridData decomp = getGrid();
	BOOST_CHECK(4==decomp.numPoints);
	PDNEIGH::NeighborhoodList list = createNeighborhood(decomp);
	Array<int> sharedGIDs(list.get_num_shared_points(),list.get_shared_gids());
	BOOST_CHECK(5==sharedGIDs.get_size());

	Field<char> ownedBc(BC_MASK,decomp.numPoints);
	ownedBc.set(2*myRank);
	Field<char> overlapField = PDNEIGH::createOverlapField<char>(list,ownedBc);
	size_t num_overlap_points = overlapField.get_num_points();
	BOOST_CHECK(9==num_overlap_points);
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "createDistributor" );
	proc->add(BOOST_TEST_CASE( &createDistributor ));
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
	 * This test only make sense for numProcs == 4
	 */
	if(4 != numProcs){
		std::cerr << "Unit test runtime ERROR: ut_Epetra_Distributor only makes sense on 4 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 4 ./ut_Epetra_Distributor." << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
