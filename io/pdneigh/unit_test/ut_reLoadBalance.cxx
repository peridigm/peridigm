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

#include "utilities/Array.h"
#include "../PdZoltan.h"
#include "../BondFilter.h"
#include "../NeighborhoodList.h"
#include "mesh_input/quick_grid/QuickGrid.h"

#include "mesh_output/vtk/Field.h"
#include "utilities/PdutMpiFixture.h"


using std::tr1::shared_ptr;
using namespace boost::unit_test;

using namespace Pdut;
using std::cout;
using std::set;
using std::map;

static size_t myRank;
static size_t numProcs;

void assert_grid(QUICKGRID::QuickGridData decomp, const std::string& label="");
using QUICKGRID::QuickGridMeshGenerationIterator;

QUICKGRID::QuickGridData getGrid(const string json_filename) {
	shared_ptr<QuickGridMeshGenerationIterator> g = QUICKGRID::getMeshGenerator(numProcs,json_filename);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, *g);

	// This load-balances
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
	return decomp;
}

void assert_grid(QUICKGRID::QuickGridData decomp,const std::string& label) {
	std::cout << label << "\n";
	int dimension_answer=3;
	int globalNumPoints_answer=2;
	int numPoints_answer=2;
	int sizeNeighborhoodList_answer=4;
	int gids_answer[] = {0,1};
	int neighborhood_answer[] = {1,1,1,0};
	int neighborhoodPtr_answer[] = {0,2};

	BOOST_CHECK(dimension_answer==decomp.dimension);
	BOOST_CHECK(globalNumPoints_answer==decomp.globalNumPoints);
	BOOST_CHECK(numPoints_answer==decomp.numPoints);
	BOOST_CHECK(sizeNeighborhoodList_answer==decomp.sizeNeighborhoodList);
	//	BOOST_CHECK(0==decomp.numExport);
	//	BOOST_CHECK(decomp.unPack);
	{
		const int* gids=decomp.myGlobalIDs.get();
		for(int n=0;n<numPoints_answer;n++,gids++)
			BOOST_CHECK(gids_answer[n]==*gids);
	}
	{
		const int* neighborhood = decomp.neighborhood.get();
		for(int n=0;n<sizeNeighborhoodList_answer;n++,neighborhood++)
			BOOST_CHECK(neighborhood_answer[n]==*neighborhood);
	}
	{
		const int* neighborhoodPtr = decomp.neighborhoodPtr.get();
		for(int n=0;n<numPoints_answer;n++,neighborhoodPtr++)
			BOOST_CHECK(neighborhoodPtr_answer[n]==*neighborhoodPtr);
	}

}

void twoPointReloadBalance() {
	const string json_file="input_files/ut_twoPointReLoadBalance.json";
	QUICKGRID::QuickGridData decomp_1 = getGrid(json_file);
	assert_grid(decomp_1,"\ttwoPointReloadBalance: UNBALANCED MESH assert_grid()\n");

	/*
	 * Reload balance
	 */
	QUICKGRID::QuickGridData decomp_2 = PDNEIGH::getLoadBalancedDiscretization(decomp_1);
	assert_grid(decomp_1,"\ttwoPointReloadBalance: LOADBALANCED MESH assert_grid()\n");

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

