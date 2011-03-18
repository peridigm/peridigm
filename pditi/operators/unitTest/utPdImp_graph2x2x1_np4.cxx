/*
 * utPimp_utPimp_graph2x2x1_np4.cxx
 *
 *  Created on: Jun 3, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "Sortable.h"
#include "Array.h"
#include "quick_grid/QuickGrid.h"
#include "NeighborhoodList.h"
#include "PdZoltan.h"
#include "vtk/Field.h"
#include "vtk/PdVTK.h"
#include "../PdImpMaterials.h"
#include "../PdITI_Operator.h"
#include "../PdITI_Utilities.h"
#include "../DirichletBcSpec.h"
#include "../BodyLoadSpec.h"
#include "../StageFunction.h"
#include "../Loader.h"
#include "../StageComponentDirichletBc.h"
#include "../ComponentDirichletBcSpec.h"
#include "../IsotropicElasticConstitutiveModel.h"
#include "../ConstitutiveModel.h"
#include "PdutMpiFixture.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <set>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include <utility>
#include <map>
#include <vector>
#include <cstdlib>
#include <time.h>

using UTILITIES::CartesianComponent;
using namespace Pdut;
using namespace Field_NS;
using std::tr1::shared_ptr;
using namespace boost::unit_test;


static size_t myRank;
static size_t numProcs;

const size_t nx = 2;
const size_t ny = nx;
const double lX = 2.0;
const double lY = lX;
const double lZ = 1.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const size_t nz = 1;
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const size_t numCells = nx*ny*nz;
const double horizon=1.1;
const PdImp::BulkModulus _K(130000.0);
const PdImp::PoissonsRatio _MU(0.0);
const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
const double g = 9.807e-3;
const double rho = 7800e-6;
const int vectorNDF=3;

using PdVTK::writeField;
using PdImp::ComponentDirichletBcSpec;
using PdImp::StageFunction;
using PdImp::StageComponentDirichletBc;
using std::cout;
using std::endl;

//shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm);
shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian);

void testGraph() {
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,QUICKGRID::SphericalNorm);
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

	/*
	 * Each processor should own 1 point
	 */
	BOOST_CHECK(decomp.numPoints==1);

	/*
	 * Create Pimp Operator
	 */
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);
	PdITI::PdITI_Operator op(comm,list,decomp.cellVolume);

	/*
	 * There should be 3 points total in the overlap vectors
	 * Each point only has 3 neighbors -- it does not have the diagonal neighbor
	 */
	Epetra_BlockMap overlapMap = list.getOverlapMap(comm,3);
	BOOST_CHECK(3==overlapMap.NumMyElements());

	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	Field<double> uOverlapField(DISPL3D,overlapMap.NumMyElements());
	uOverlapField.set(0.0);
	shared_ptr<RowStiffnessOperator> jacobian = op.getJacobian(uOverlapField);

	/*
	 * Get points for bc's
	 */
	CartesianComponent axis = UTILITIES::Z;
	UTILITIES::Array<int> bcIds = UTILITIES::getPointsAxisAlignedMaximum(axis,decomp.myX,decomp.numPoints,horizon);
	std::sort(bcIds.get(),bcIds.end());

	/**
	 * Create boundary conditions spec
	 */
	ComponentDirichletBcSpec allFixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIds);
	StageFunction constStageFunction(0.0,0.0);
	StageComponentDirichletBc bcOp(allFixedSpec,constStageFunction);

	/*
	 * Create graph
	 */
	shared_ptr<Epetra_CrsGraph> graphPtr = getGraph(jacobian);

}


//shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm) {
//	PdImp::PdImpOperator *op = new PdImp::PdImpOperator(comm,decomp);
//	return shared_ptr<PdImp::PdImpOperator>(op);
//}

shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian){
	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();
	const Epetra_BlockMap& colMap = jacobian->getColMap();

	/*
	 * Epetra Graph
	 */
	UTILITIES::Array<int> numCols = jacobian->getNumColumnsPerRow();
	Epetra_CrsGraph *graph = new Epetra_CrsGraph(Copy,rowMap,colMap,numCols.get());
	for(int row=0;row<jacobian->getNumRows();row++){
		UTILITIES::Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.get_size();
		/*
		 * Each row should have 4 entries
		 */
		BOOST_CHECK(4==numCol);
		int *cols = rowLIDs.get();
		BOOST_CHECK(0==graph->InsertMyIndices(row,numCol,cols));
	}
	BOOST_CHECK(0==graph->FillComplete());
	return shared_ptr<Epetra_CrsGraph>(graph);
}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "testGraph" );
	proc->add(BOOST_TEST_CASE( &testGraph ));
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
		std::cerr << "Unit test runtime ERROR: utPimp_graph2x2x1_np4 only makes sense on 4 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 4 ./utPimp_graph2x2x1_np4." << std::endl;
		myMpi.PdutMpiFixture::~PdutMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
