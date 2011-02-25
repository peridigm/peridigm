/*
 * ut_frameset_2x2x1_.cxx
 *
 *  Created on: Feb 25, 2011
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
#include "../NeighborhoodList.h"
#include "zoltan.h"
#include "PdZoltan.h"
#include "Field.h"
#include "PdBondFilter.h"
#include <Teuchos_RCP.hpp>
#include "mpi.h"
#include <tr1/memory>
#include <valarray>
#include <iostream>
#include <cmath>
#include <map>
#include <set>


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

using Teuchos::RCP;
using PdBondFilter::BondFilter;


using namespace Pdut;
using std::tr1::shared_ptr;
using std::cout;
using std::set;
using std::map;


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using namespace Field_NS;
using std::tr1::shared_ptr;
using namespace boost::unit_test;


static int numProcs;
static int myRank;

const int nx = 2;
const int ny = nx;
const double lX = 2.0;
const double lY = lX;
const double lZ = 1.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const int nz = 1;
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=1.1;
using std::cout;
using std::endl;



PdGridData getGrid() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	// This load-balances
	decomp = getLoadBalancedDiscretization(decomp);
	return decomp;
}


shared_ptr< std::set<int> > constructFrame(PDNEIGH::NeighborhoodList& list) {
	shared_ptr< std::set<int> > frameSet = list.constructParallelDecompositionFrameSet();
	/*
	 * There is only 1 point on each processor and by default it must show up in the frameset
	 */
	BOOST_CHECK(1==frameSet->size());
	return frameSet;
}

void createNeighborhood() {
	PdGridData decomp = getGrid();
	BOOST_CHECK(1==decomp.numPoints);
	BOOST_CHECK(4==decomp.globalNumPoints);
	PDNEIGH::NeighborhoodList list(decomp.numPoints,decomp.myGlobalIDs,decomp.myX,2.0*horizon);
	shared_ptr< std::set<int> > frameSet = constructFrame(list);
	RCP<BondFilter> bondFilterPtr(new PdBondFilter::BondFilterWithSelf());
	decomp = createAndAddNeighborhood(decomp,2*horizon,bondFilterPtr);

	BOOST_CHECK(1==decomp.numPoints);
	BOOST_CHECK(4==decomp.globalNumPoints);
	BOOST_CHECK(5==decomp.sizeNeighborhoodList);
	/*
	 * Neighborhood of every point should have every other point
	 */
	int n[] = {0,1,2,3};
	set<int> neighSet(n,n+4);
	int *neighborhood = decomp.neighborhood.get();
	int numNeigh = *neighborhood;
	BOOST_CHECK(4==numNeigh);
	neighborhood++;
	for(int i=0;i<numNeigh;i++,neighborhood++){
		BOOST_CHECK(neighSet.end()!=neighSet.find(*neighborhood));
	}

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "createNeighborhood" );
	proc->add(BOOST_TEST_CASE( &createNeighborhood ));
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
		std::cerr << "Unit test runtime ERROR: ut_frameset_2x2x1_np4 only makes sense on 4 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 4 ./ut_frameset_2x2x1_np4." << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
