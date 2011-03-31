/*
 * utCrackOpeningDemo.cxx
 *
 *  Created on: Feb 11, 2011
 *      Author: jamitch
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_ALTERNATIVE_INIT_API
//#include <boost/test/unit_test.hpp>
//#include <boost/test/parameterized_test.hpp>
#include "Sortable.h"
#include "Array.h"
#include "quick_grid/QuickGrid.h"
#include "NeighborhoodList.h"
#include "OverlapDistributor.h"
#include "BondFilter.h"
#include "PdZoltan.h"
#include "vtk/Field.h"
#include "vtk/PdVTK.h"
#include "../../PdImpMaterials.h"
#include "../../PdITI_Operator.h"
#include "../../PdITI_Utilities.h"
#include "../../DirichletBcSpec.h"
#include "../../BodyLoadSpec.h"
#include "../../StageFunction.h"
#include "../../Loader.h"
#include "../../StageComponentDirichletBc.h"
#include "../../ComponentDirichletBcSpec.h"
#include "../../IsotropicElasticConstitutiveModel.h"
#include "../../ConstitutiveModel.h"
#include "../utPdITI.h"
#include "PdutMpiFixture.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include <set>
#include <Teuchos_RCPDecl.hpp>
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
using UTILITIES::Array;
using PdITI::ConstitutiveModel;
using PdITI::IsotropicElasticConstitutiveModel;
using PdImp::StageComponentDirichletBc;
using PdImp::StageDirichletBc;
using PdImp::DirichletBcSpec;
using PdImp::ComponentDirichletBcSpec;
using PdImp::StageFunction;
using namespace PdBondFilter;
using namespace Pdut;
using namespace Field_NS;
using std::tr1::shared_ptr;
using Teuchos::RCP;

const int vectorNDF=3;
static size_t myRank;
static size_t numProcs;

/*
 * This should be even so that the crack plane lies between to rows of points
 */
const size_t nx = 30;
const size_t ny = 30;
const double xStart = -2.5;
const double xLength = 5.0;
const double yStart = -2.5;
const double yLength = 5.0;
const double zStart = -0.5;
const double zLength = 1.0;
static double xMax=xStart+xLength;
static double xMin=xStart;
const size_t nz = nx* (zLength / xLength);
const size_t numCells = nx*ny*nz;
const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
const double dx = xSpec.getCellSize();
const double dy = ySpec.getCellSize();
const double dz = zSpec.getCellSize();
const double _cellVolume = dx*dy*dz;

/*
 * Young's Modulus (MPa)
 */
const double E = 68.9e3;

/*
 * Poisson's ratio
 */
const double nu = 0.0;

/*
 * Density of aluminum g/mm^3
 */
const double rho = 2.7e-3;

/*
 * Horizon
 */
const double horizon=1.1*sqrt( (3.0*dx)*(3.0*dx) );

/*
 * Function prototypes in this file
 */
FinitePlane getYZ_CrackPlane();

QUICKGRID::QuickGridData getGrid() {

	/*
	 * This demonstrates how the first and second coordinate
	 * along an axis are computed
	 * const double x0 = xStart+xSpec.getCellSize()/2.0;
	 * const double x1 = x0 + xSpec.getCellSize();
	 * const double y0 = yStart+ySpec.getCellSize()/2.0;
	 * const double y1 = y0 + ySpec.getCellSize();
	 * const double z0 = zStart+zSpec.getCellSize()/2.0;
	 * const double z1 = z0 + zSpec.getCellSize();
	*/

	if(0==myRank){
		cout << "Creating and load balancing mesh..." << endl;
	}

	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	QUICKGRID::QuickGridData gridData =  QUICKGRID::getDiscretization(myRank, cellPerProcIter);
	gridData=PDNEIGH::getLoadBalancedDiscretization(gridData);

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
	QUICKGRID::QuickGridData gridData = getGrid();
	/*
	 * Communicator
	 */
	shared_ptr<Epetra_Comm> comm = shared_ptr<Epetra_Comm>(new Epetra_MpiComm(MPI_COMM_WORLD));

	FinitePlane crackPlane=getYZ_CrackPlane();
//	shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new BondFilterDefault());
	shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	PDNEIGH::NeighborhoodList list(comm,gridData.zoltanPtr.get(),gridData.numPoints,gridData.myGlobalIDs,gridData.myX,horizon,filterPtr);

	/*
	 * Material Properties
	 */
	IsotropicHookeSpec isotropicSpec = utPdITI::getMaterialSpec(E,nu);

	/*
	 * Create PdITI Operator
	 */
	PdITI::PdITI_Operator op(comm,list,gridData.cellVolume);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);
	PDNEIGH::NeighborhoodList row_matrix_list = op.get_row_matrix_neighborhood();


	/*
	 * Get points for bc's
	 */
	/*
	 * Note that we are looking for a discrete number of points at end;
	 * Set the scale factor to just larger than an integer where
	 * the integer corresponds with the number of points to be included in
	 * the boundary conditions
	 */
	double scaleFactor=2.1;
	double searchRadius=scaleFactor*dx;
	CartesianComponent axis = UTILITIES::X;
	Array<int> bcIdsFixed = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMin);
	Array<int> bcIdsApplied = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,gridData.myX,gridData.numPoints,searchRadius,xMax);

	/**
	 * Create boundary conditions spec
	 */
	vector<shared_ptr<StageComponentDirichletBc> > bcs;
	ComponentDirichletBcSpec fixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIdsFixed);
	StageFunction constStageFunction(0.0,0.0);
	shared_ptr<StageComponentDirichletBc> bcFixed(new StageComponentDirichletBc(fixedSpec,constStageFunction));
	bcs.push_back(bcFixed);
	std::vector< DirichletBcSpec::ComponentLabel > c(1);
	c[0] = DirichletBcSpec::X;
	ComponentDirichletBcSpec appliedSpec(c,bcIdsApplied);
	StageFunction dispStageFunction(1.0e-3,1.0e-3);
	shared_ptr<StageComponentDirichletBc> bcApplied(new StageComponentDirichletBc(appliedSpec,dispStageFunction));
	bcs.push_back(bcApplied);
	Field<char> bcMaskFieldOwned(BC_MASK,gridData.numPoints);
	bcMaskFieldOwned.set(0);
	for(int b=0;b<bcs.size();b++)
		bcs[b]->imprint_bc(bcMaskFieldOwned);
	Field<char> bcMaskFieldOverlap = PDNEIGH::createOverlapField(row_matrix_list,bcMaskFieldOwned);

	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	Field<double> uOwnedField(DISPL3D,gridData.numPoints);
	uOwnedField.set(0.0);
	for(int b=0;b<bcs.size();b++)
		bcs[b]->applyKinematics(1.0,uOwnedField);
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getJacobian(uOwnedField);
	Epetra_BlockMap ownedMap = jacobian->getRowMap();


	/*
	 * Create Epetra_RowMatrix
	 */
	shared_ptr<Epetra_RowMatrix> mPtr = utPdITI::getOperator(bcMaskFieldOverlap,jacobian);

	/*
	 * TODO
	 * Investigate SIGN of applied loading
	 * 1) Compute internal force with displacement vector that has kinematics applied
	 * 2) Apply kinematics to this vector so that solution properly includes the applied kinematics
	 */
	Field_NS::FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fNOwnedField(fNSpec,gridData.numPoints);
	fNOwnedField.set(0.0);
	Field_NS::FieldSpec deltaSpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D,"delta");
	Field_NS::Field<double> delta(deltaSpec,gridData.numPoints);
	delta.set(0.0);
	op.computeInternalForce(uOwnedField,fNOwnedField);
	for(int b=0;b<bcs.size();b++){
		bcs[b]->applyKinematics(1.0,fNOwnedField);
		bcs[b]->applyKinematics(1.0,delta);
	}


	/*
	 * AZTEC Setup (THIS WORKS FINE)
	 */
//	Epetra_LinearProblem linProblem;
//	linProblem.SetOperator(mPtr.get());
//	linProblem.AssertSymmetric();
//
//	/*
//	 * Domain and Range Map are the same
//	 */
//	const Epetra_BlockMap& rangeMap  = mPtr->OperatorRangeMap();
//	const Epetra_BlockMap& domainMap = mPtr->OperatorDomainMap();
//	Epetra_Vector rhs(View,rangeMap,fNOwnedField.get());
//	Epetra_Vector lhs(View,domainMap,uOwnedField.get());
//	linProblem.SetRHS(&rhs);
//	linProblem.SetLHS(&lhs);
//	if(0 != linProblem.CheckInput()){
//		cout << "0 != linProblem.CheckInput()" << endl;
//		std::exit(1);
//	}
//
//	AztecOO solver(linProblem);
//	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
//	if(0 != solver.CheckInput()){
//		cout << "0 != solver.CheckInput()" << endl;
//		std::exit(1);
//	}
//	solver.Iterate(500,1e-6);
	/*
	 * END AZTEC Setup (THIS WORKS FINE)
	 */

	/*
	 * BELOS Setup (THIS WORKS FINE)
	 */
	ParameterList belosList;
	int blocksize=1;
	int maxiters=500;
	double tol=1.0e-6;
	int frequency=1;
	belosList.set( "Block Size", blocksize );              // Blocksize to be used by iterative solver
	belosList.set( "Maximum Iterations", maxiters );       // Maximum number of iterations allowed
	belosList.set( "Convergence Tolerance", tol );         // Relative convergence tolerance requested
	belosList.set( "Verbosity", Belos::Errors + Belos::Warnings +
	Belos::TimingDetails + Belos::StatusTestDetails );
	belosList.set( "Output Frequency", frequency );
	belosList.set("Output Style", Belos::Brief);
	//
	// Construct an unpreconditioned linear problem instance.
	//
	typedef Epetra_MultiVector                MV;
	typedef Epetra_Operator                   OP;
	RCP<Epetra_Vector> B = rcp(new Epetra_Vector(View,ownedMap,fNOwnedField.get()));
	RCP<Epetra_Vector> X = rcp(new Epetra_Vector(View,ownedMap,uOwnedField.get()));
	RCP<OP> A = rcp(mPtr.get(),false);
	Belos::LinearProblem<double,MV,OP> problem( A, X, B );
	bool set = problem.setProblem();
	if (set == false) {
		std::cout << std::endl << "ERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
		std::exit(1);
	}
	RCP< Belos::SolverManager<double,MV,OP> > newSolver
	= rcp( new Belos::BlockCGSolMgr<double,MV,OP>(rcp(&problem,false), rcp(&belosList,false)) );
	//
	// Perform solve
	//
	Belos::ReturnType ret = newSolver->solve();
	/*
	 * END BELOS Setup (THIS WORKS FINE)
	 */


	/*
	 * Write problem set up parameters to file
	 */
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(gridData.myX,gridData.numPoints);

	/*
	 * Add processor rank to output
	 */
	Field_NS::Field<int> myRank(Field_NS::PROC_NUM,gridData.numPoints);
	myRank.set(comm->MyPID());

	PdVTK::writeField(grid,fNOwnedField);
	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,delta);
	PdVTK::writeField<int>(grid,myRank);
	PdVTK::writeField<char>(grid,bcMaskFieldOwned);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utCrackOpeningDemo_npX.pvtu", comm->NumProc(), comm->MyPID());
	PdVTK::write(writer,grid);

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

	crackOpeningDemo();

}
