/*
 * utPimp_componentDirichletBcSpec.cxx
 *
 *  Created on: Apr 26, 2010
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
#include "../IsotropicElasticConstitutiveModel.h"
#include "../ConstitutiveModel.h"
#include "../DirichletBcSpec.h"
#include "../ComponentDirichletBcSpec.h"
#include "../StageComponentDirichletBc.h"
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
using PdImp::IsotropicElasticConstitutiveModel;
using PdImp::ConstitutiveModel;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::vector;
using std::cout;
using std::endl;

static int numProcs;
static int myRank;

const int nx = 4;
const int ny = 4;
const int nz = 8;
const double lX = 1.0;
const double lY = 1.0;
const double lZ = 1.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=1.01*sqrt(pow(lX/nx,2)+pow(lY/ny,2)+pow(lZ/nz,2));
const PdImp::BulkModulus _K(130000.0);
const PdImp::PoissonsRatio _MU(0.0);
const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
const int vectorNDF=3;

void assertBcIds(const Pd_shared_ptr_Array<int>& bcIds);
using PdImp::DirichletBcSpec;
using PdImp::ComponentDirichletBcSpec;
using PdImp::StageComponentDirichletBc;
using PdImp::StageFunction;

void allComponentsFixed_1() {

	/**
	 * Create array of boundary conditions
	 */
	std::vector< DirichletBcSpec::ComponentLabel > comps(3);
	comps[0] = DirichletBcSpec::X;
	comps[1] = DirichletBcSpec::Y;
	comps[2] = DirichletBcSpec::Z;
	Pd_shared_ptr_Array<int> ids(0);
	PdImp::ComponentDirichletBcSpec allFixedSpec(comps,ids);
	vector<DirichletBcSpec::ComponentLabel> c = allFixedSpec.getComponents();
	BOOST_CHECK(c.size()==3);
	BOOST_CHECK(c[0]==DirichletBcSpec::X);
	BOOST_CHECK(c[1]==DirichletBcSpec::Y);
	BOOST_CHECK(c[2]==DirichletBcSpec::Z);
	vector< vector<double> > dirs = allFixedSpec.getUnitDirections();
	BOOST_CHECK(dirs.size()==3);
	vector<double> u1 = dirs[0];
	vector<double> u2 = dirs[1];
	vector<double> u3 = dirs[2];

	/*
	 * check u1
	 */
	double tolerance=1.0e-15;
	BOOST_CHECK(u1.size()==3);
	BOOST_CHECK_CLOSE(u1[0],1.0,tolerance);
	BOOST_CHECK_CLOSE(u1[1],0.0,tolerance);
	BOOST_CHECK_CLOSE(u1[2],0.0,tolerance);

	/*
	 * check u2
	 */
	BOOST_CHECK(u2.size()==3);
	BOOST_CHECK_CLOSE(u2[0],0.0,tolerance);
	BOOST_CHECK_CLOSE(u2[1],1.0,tolerance);
	BOOST_CHECK_CLOSE(u2[2],0.0,tolerance);

	/*
	 * check u3
	 */
	BOOST_CHECK(u3.size()==3);
	BOOST_CHECK_CLOSE(u3[0],0.0,tolerance);
	BOOST_CHECK_CLOSE(u3[1],0.0,tolerance);
	BOOST_CHECK_CLOSE(u3[2],1.0,tolerance);


}

void allComponentsFixed_2() {

	Pd_shared_ptr_Array<int> ids(0);
	ComponentDirichletBcSpec allFixedSpec = ComponentDirichletBcSpec::getAllComponents(ids);

	vector<DirichletBcSpec::ComponentLabel> c = allFixedSpec.getComponents();
	BOOST_CHECK(c.size()==3);
	BOOST_CHECK(c[0]==DirichletBcSpec::X);
	BOOST_CHECK(c[1]==DirichletBcSpec::Y);
	BOOST_CHECK(c[2]==DirichletBcSpec::Z);

	vector< vector<double> > dirs = allFixedSpec.getUnitDirections();
	BOOST_CHECK(dirs.size()==3);
	vector<double> u1 = dirs[0];
	vector<double> u2 = dirs[1];
	vector<double> u3 = dirs[2];

	/*
	 * check u1
	 */
	double tolerance = 1.0e-15;
	BOOST_CHECK(u1.size()==3);
	BOOST_CHECK_CLOSE(u1[0],1.0,tolerance);
	BOOST_CHECK_CLOSE(u1[1],0.0,tolerance);
	BOOST_CHECK_CLOSE(u1[2],0.0,tolerance);

	/*
	 * check u2
	 */
	BOOST_CHECK(u2.size()==3);
	BOOST_CHECK_CLOSE(u2[0],0.0,tolerance);
	BOOST_CHECK_CLOSE(u2[1],1.0,tolerance);
	BOOST_CHECK_CLOSE(u2[2],0.0,tolerance);

	/*
	 * check u3
	 */
	BOOST_CHECK(u3.size()==3);
	BOOST_CHECK_CLOSE(u3[0],0.0,tolerance);
	BOOST_CHECK_CLOSE(u3[1],0.0,tolerance);
	BOOST_CHECK_CLOSE(u3[2],1.0,tolerance);


}

void applyHomogeneousForm() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	int numPoints = decomp.numPoints;
	BOOST_CHECK(numCells==numPoints);
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);


	Field<double> uOwnedField = PdImp::getPureShearXY(Field_NS::Field<double>(Field_NS::COORD3D,decomp.myX,numPoints));

	/*
	 * Compute Tangent
	 * For each row, probe force operator and compare
	 */
	FieldSpec::FieldSpec fNP1Spec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fNP1");
	FieldSpec::FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::TemporalField<double> force = Field_NS::TemporalField<double>(Field_NS::FORCE3D,numPoints);
	Field_NS::Field<double> fN(fNSpec,numPoints);
	Field_NS::Field<double> fNP1(fNP1Spec,numPoints);
	PdImp::PdImpOperator op(comm,decomp);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);

	op.computeInternalForce(uOwnedField,fN);
	op.computeInternalForce(uOwnedField,fNP1);

	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::Y;
	Pd_shared_ptr_Array<int> bcIds = PdNeighborhood::getPointsAxisAlignedMinimum(axis,decomp.myX,numPoints,horizon);
	std::sort(bcIds.get(),bcIds.get()+bcIds.getSize());
	assertBcIds(bcIds);
	ComponentDirichletBcSpec allFixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIds);
	StageFunction constStageFunction(0.0,0.0);
	StageComponentDirichletBc bcOp(allFixedSpec,constStageFunction);
	bcOp.applyHomogeneousForm(fNP1);

	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);
	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,fN);
	PdVTK::writeField(grid,fNP1);

	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("applyHomogeneousForm.pvtu", comm.NumProc(), comm.MyPID());
	PdVTK::write(writer,grid);

}

/*
 * Note that for this to work, incoming bcIds must be sorted
 */
void assertBcIds(const Pd_shared_ptr_Array<int>& bcIds){

	int delta=nx*ny;
	int numPointsPerXY=8;

	BOOST_CHECK(64==bcIds.getSize());
	const int *ids = bcIds.get();
	int c=0;
	for(int j=0;j<nz;j++,c+=delta){
		for(int p=0;p<numPointsPerXY;p++,ids++)
		BOOST_CHECK(c+p==*ids);
	}

}


bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "ComponentDirichletBcSpec" );
	proc->add(BOOST_TEST_CASE( &allComponentsFixed_1 ));
	proc->add(BOOST_TEST_CASE( &allComponentsFixed_2 ));
	proc->add(BOOST_TEST_CASE( &applyHomogeneousForm ));
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
		std::cerr << "Unit test runtime ERROR: utPimp_componentDirichletBcSpec is intended for \"serial\" run only and makes sense on 1 processor" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 1 ./utPimp_componentDirichletBcSpec" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
