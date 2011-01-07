/*
 * utPimp_constantBodyLoad.cxx
 *
 *  Created on: Apr 27, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "../PdImpMpiFixture.h"
#include "PdNeighborhood.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdNeighborhood.h"
#include "Field.h"
#include "../PdImpMaterials.h"
#include "../PdImpOperator.h"
#include "../PdImpOperatorUtilities.h"
#include "../DirichletBcSpec.h"
#include "../BodyLoadSpec.h"
#include "../StageFunction.h"
#include "../Loader.h"
#include "../StageComponentDirichletBc.h"
#include "../ComponentDirichletBcSpec.h"
#include "PdVTK.h"
#include <set>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <utility>
#include <map>
#include <vector>
#include <cstdlib>
#include <set>
#include <time.h>
#include <tr1/memory>

using namespace PdQuickGrid;
using namespace PdNeighborhood;
using namespace Field_NS;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::vector;


static int numProcs;
static int myRank;

const int nx = 5;
const int ny = nx;
const double lX = 1.0;
const double lY = lX;
const double lZ = 10.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const int nz = (int)(lZ * nx / lX);
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=1.5*sqrt(pow(lX/nx,2)+pow(lY/ny,2)+pow(lZ/nz,2));
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

void constantBodyLoad() {


	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	int numPoints = decomp.numPoints;
	BOOST_CHECK(numCells==numPoints);
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);

	/*
	 * Create force field to assemble body load into
	 */
	FieldSpec::FieldSpec fNP1Spec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fNP1");
	FieldSpec::FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fN(fNSpec,numPoints);
	Field_NS::Field<double> fNP1(fNP1Spec,numPoints);

	/*
	 * Get points for bc's
	 */
	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::Z;
	Pd_shared_ptr_Array<int> bcIds = PdNeighborhood::getPointsAxisAlignedMaximum(axis,decomp.myX,numPoints,horizon);
	/**
	 * Create array of boundary conditions
	 */
	ComponentDirichletBcSpec allFixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIds);
	StageFunction constStageFunction(0.0,0.0);
	StageComponentDirichletBc bcOp(allFixedSpec,constStageFunction);


	/*
	 * Create body load
	 */
	Pd_shared_ptr_Array<int> localIds(decomp.numPoints);
	{
		/*
		 * Create list of local ids
		 */
		int *ids=localIds.get();
		const int *end=localIds.end();
		for(int i=0;ids!=end;i++,ids++)
			*ids = i;
	}
	double u[] = {0,0,-1};
	double uMag = g*rho;
	PdImp::BodyLoadSpec bodyLoadSpec(u,localIds);
	StageFunction bodyLoadFunction(0,uMag);
	shared_ptr<PdImp::Loader> bodyLoad = bodyLoadSpec.getStageLoader(bodyLoadFunction);
	bodyLoad->computeOwnedExternalForce(1.0,fN);
	bodyLoad->computeOwnedExternalForce(1.0,fNP1);
	bcOp.applyHomogeneousForm(fNP1);

	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);
	PdVTK::writeField(grid,fN);
	PdVTK::writeField(grid,fNP1);

	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("constantBodyLoad.pvtu", comm.NumProc(), comm.MyPID());
	PdVTK::write(writer,grid);

}



bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPimp_constantBodyLoad" );
	proc->add(BOOST_TEST_CASE( &constantBodyLoad ));
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
	PdImpRunTime::PimpMpiFixture pimpMPI = PdImpRunTime::PimpMpiFixture::getPimpMPI(argc,argv);
	const Epetra_Comm& comm = pimpMPI.getEpetra_Comm();

	// These are static (file scope) variables
	myRank = comm.MyPID();
	numProcs = comm.NumProc();
	/**
	 * This test only make sense for numProcs == 1
	 */
	if(1 != numProcs){
		std::cerr << "Unit test runtime ERROR: utPimp_constantBodyLoad is intended for \"serial\" run only and makes sense on 1 processor" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 1 ./utPimp_constantBodyLoad" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
