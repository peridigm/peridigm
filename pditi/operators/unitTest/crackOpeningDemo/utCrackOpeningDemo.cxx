/*
 * utCrackOpeningDemo.cxx
 *
 *  Created on: Feb 11, 2011
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdBondFilter.h"
#include "PdVTK.h"
#include "Field.h"
#include "../utPdITI.h"
#include "../../PdImpMpiFixture.h"
#include "../../PdImpMaterials.h"
#include "../../PdImpOperator.h"
#include "../../IsotropicElasticConstitutiveModel.h"
#include "../../DirichletBcSpec.h"
#include "../../StageFunction.h"
#include "../../StageComponentDirichletBc.h"
#include "../../ComponentDirichletBcSpec.h"
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include <iostream>


using namespace PdQuickGrid;
using namespace PdBondFilter;
using namespace PdVTK;
using namespace Field_NS;
using namespace PdImp;
using PdITI::IsotropicElasticConstitutiveModel;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;
using std::endl;

static int numProcs;
static int myRank;
/*
 * This should be even so that the crack plane lies between to rows of points
 */
const int nx = 6;
const int ny = 6;
const double xStart = -2.5;
const double xLength = 5.0;
const double yStart = -2.5;
const double yLength = 5.0;
const double zStart = -0.5;
const double zLength = 1.0;
static double xMax=xStart+xLength;
static double xMin=xStart;
const int nz = nx* (zLength / xLength);
const int numCells = nx*ny*nz;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const double dx = xSpec.getCellSize();
const double dy = ySpec.getCellSize();
const double dz = zSpec.getCellSize();
const double _cellVolume = dx*dy*dz;

/*
 * Horizon
 */
const double horizon=1.1*sqrt( (3.0*dx)*(3.0*dx) );

/*
 * Function prototypes in this file
 */
FinitePlane getYZ_CrackPlane();

/*
 * Young's Modulus (MPa)
 */
static double E = 68.9e3;

/*
 * Poisson's ratio
 */
static double nu = 0.0;

/*
 * Density of aluminum g/mm^3
 */
static double rho = 2.7e-3;




/*
 * This demonstrates how the first and second coordinate
 * along an axis are computed
 */
//const double x0 = xStart+xSpec.getCellSize()/2.0;
//const double x1 = x0 + xSpec.getCellSize();
//const double y0 = yStart+ySpec.getCellSize()/2.0;
//const double y1 = y0 + ySpec.getCellSize();
//const double z0 = zStart+zSpec.getCellSize()/2.0;
//const double z1 = z0 + zSpec.getCellSize();


PdGridData getGrid() {

	if(0==myRank){
		cout << "Creating and load balancing mesh..." << endl;
	}

	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData gridData =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	gridData=getLoadBalancedDiscretization(gridData);
	FinitePlane crackPlane=getYZ_CrackPlane();
	RCP<BondFilter> filterPtr=RCP<BondFilter>(new FinitePlaneFilter(crackPlane));
	gridData = createAndAddNeighborhood(gridData,horizon,filterPtr);

	/*
	 * Lower left hand corner of crack plane when viewing down
	 * normal in the +dir
	 */
	const double x0 = xStart+xLength/2;
	const double y0 = yStart;
	const double z0 = zStart;

	if(0==myRank){
		cout << "\t\tDONE." << endl;
		cout << "Total number of points in mesh = " << gridData.globalNumPoints << endl;
		cout << "nx,ny,nz = " << nx << ", " << ny << ", "<< nz << ", "<< endl;
		cout << "x0,y0,z0 = " << x0 << ", " << y0 << ", "<< z0 << ", "<< endl;
	}

	/*
	 * Write file for debugging
	 */
//	const FieldSpec myRankSpec(FieldSpec::DEFAULT_FIELDTYPE,FieldSpec::SCALAR,"MyRank");
//	Field<double> X(COORD3D,gridData.myX,gridData.numPoints);
//	Field<int> rankField(myRankSpec,gridData.numPoints);
//	rankField.setValue(myRank);
//
//	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(gridData.myX.get(), gridData.numPoints);
//	PdVTK::writeField(grid,X);
//	PdVTK::writeField(grid,rankField);
//	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer= PdVTK::getWriter("utCrackOpeningDemo.pvtu", numProcs, myRank, PdVTK::vtkBINARY);
//	PdVTK::write(writer,grid);

	return gridData;
}



FinitePlane getYZ_CrackPlane() {

	/*
	 * Crack runs along y-axis
	 * Crack length along y-axis is yLength/2
	 * Crack runs from bottom to top of plate in z-dir; a=zLength
	 */

	/*
	 * Lower left hand corner of crack plane when viewing down
	 * normal in the +dir
	 */
	const double x0 = xStart+xLength/2;
	const double y0 = yStart;
	const double z0 = zStart;

	/*
	 * normal is along -x-dir
	 */
	double n[3]; n[0]=-1.0;n[1]=0.0;n[2]=0.0;
	/*
	 * lower left corner of plane
	 */
	double r0[3]; r0[0]=x0; r0[1]=y0; r0[2]=z0;
	/*
	 * vector along bottom edge is in the +y-dir
	 */
	double ub[3]; ub[0]=0; ub[1]=1.0;ub[2]=0.0;
	/*
	 * b is length of crack along bottom edge
	 * a is length of crack along z-dir
	 */
	double b=yLength/2.0, a=zLength;
	return FinitePlane(n,r0,ub,b,a);
}


void crackOpeningDemo(){
	/*
	 * Get mesh and decomposition
	 */
	PdGridData gridData = getGrid();
	/*
	 * Communicator
	 */
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);

	/*
	 * Material Properties
	 */
	IsotropicHookeSpec isotropicSpec = utPdITI::getMaterialSpec(E,nu);

	/*
	 * Create PdITI Operator
	 */
	shared_ptr<PdImp::PdImpOperator> op = utPdITI::getPimpOperator(gridData,comm);
	shared_ptr<PdITI::ConstitutiveModel> fIntOperator(new PdITI::IsotropicElasticConstitutiveModel(isotropicSpec));
	op->addConstitutiveModel(fIntOperator);

	/*
	 * Get points for bc's
	 */
	/*
	 * Note that we are looking for the first three planes of points on ends
	 */
	double searchRadius=3.1*dx;
	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::X;
	Pd_shared_ptr_Array<int> bcIdsFixed = PdNeighborhood::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMin);
	std::sort(bcIdsFixed.get(),bcIdsFixed.get()+bcIdsFixed.getSize());
	Pd_shared_ptr_Array<int> bcIdsApplied = PdNeighborhood::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMax);
	std::sort(bcIdsApplied.get(),bcIdsApplied.get()+bcIdsApplied.getSize());

	/**
	 * Create boundary conditions spec
	 */
	vector<shared_ptr<StageComponentDirichletBc> > bcs(2);
	ComponentDirichletBcSpec fixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIdsFixed);
	StageFunction constStageFunction(0.0,0.0);
	shared_ptr<StageComponentDirichletBc> bcFixed(new StageComponentDirichletBc(fixedSpec,constStageFunction));
	bcs[0] = bcFixed;
	std::vector< DirichletBcSpec::ComponentLabel > c(1);
	c[0] = DirichletBcSpec::X;
	ComponentDirichletBcSpec appliedSpec(c,bcIdsApplied);
	StageFunction dispStageFunction(1.0e-3,1.0e-3);
	shared_ptr<StageComponentDirichletBc> bcApplied(new StageComponentDirichletBc(appliedSpec,dispStageFunction));
	bcs[1] = bcApplied;

	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	Field<double> uOwnedField(DISPL3D,gridData.numPoints);
	uOwnedField.setValue(0.0);
	bcApplied->applyKinematics(1.0,uOwnedField);
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op->getRowStiffnessOperator(uOwnedField,horizon);

	/*
	 * Create graph
	 */
	shared_ptr<Epetra_CrsGraph> graphPtr = utPdITI::getGraph(jacobian);

	/*
	 * Create Epetra_RowMatrix
	 */
	shared_ptr<Epetra_RowMatrix> mPtr = utPdITI::getOperator(bcs,graphPtr,jacobian);

	/*
	 * Create force field
	 * IN this case,
	 * 1) Compute internal force with displacement vector that has kinematics applied
	 * 2) Negate internal force since we are using it as a residual on the RHS (THIS IS NOT DONE -- SIGNS need to be investigated)
	 * 3) Apply kinematics to this vector so that solution properly includes the applied kinematics
	 */
	Field_NS::FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fN(fNSpec,gridData.numPoints);
	fN.setValue(0.0);
	op->computeInternalForce(uOwnedField,fN);
	bcApplied->applyKinematics(1.0,fN);

	Epetra_LinearProblem linProblem;
	linProblem.SetOperator(mPtr.get());
	linProblem.AssertSymmetric();

	const Epetra_BlockMap& rangeMap  = mPtr->OperatorRangeMap();
	const Epetra_BlockMap& domainMap  = mPtr->OperatorDomainMap();

	double *f = fN.getArray().get();
	Epetra_Vector rhs(View,rangeMap,f);
	Epetra_Vector lhs(View,domainMap,uOwnedField.getArray().get());

	linProblem.SetRHS(&rhs);
	linProblem.SetLHS(&lhs);
	BOOST_CHECK(0==linProblem.CheckInput());

	AztecOO solver(linProblem);
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
	BOOST_CHECK(0==solver.CheckInput());
	solver.Iterate(500,1e-6);
	/*
	 * Write problem set up parameters to file
	 */
	Field_NS::FieldSpec deltaSpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D,"delta");
	Field_NS::Field<double> delta(deltaSpec,gridData.numPoints);
	delta.setValue(0.0);
	bcApplied->applyKinematics(1.0,delta);

	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(gridData.myX,gridData.numPoints);

	/*
	 * Add processor rank to output
	 */
	Field_NS::Field<int> myRank(Field_NS::PROC_NUM,gridData.numPoints);
	myRank.setValue(comm.MyPID());

	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,delta);
	writeField<int>(grid,myRank);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utPimp_pullCylinder_npX.pvtu", comm.NumProc(), comm.MyPID());
	PdVTK::write(writer,grid);

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "utCrackOpeningDemo" );
	proc->add(BOOST_TEST_CASE( &crackOpeningDemo ));
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

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
