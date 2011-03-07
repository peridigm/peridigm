/*
 * utPimp_pullCylinder_npX.cxx
 *
 *  Created on: Jun 24, 2010
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
#include "PdZoltan.h"
#include "Field.h"
#include "utPdITI.h"
#include "../../pdneigh/NeighborhoodList.h"
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
#include "PdVTK.h"
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


using namespace PdQuickGrid;
using namespace PdNeighborhood;
using namespace Field_NS;
using PdITI::IsotropicElasticConstitutiveModel;
using PdITI::ConstitutiveModel;
using PdVTK::writeField;
using PdImp::DirichletBcSpec;
using PdImp::ComponentDirichletBcSpec;
using PdImp::StageFunction;
using PdImp::StageComponentDirichletBc;

using std::tr1::shared_ptr;
using namespace boost::unit_test;

using std::vector;
using std::set;

static int numProcs;
static int myRank;
static double zMax;
static double zMin;

const double cylinderRadius = 1.0;
const int numRings = 5;
const double zStart = 0.0;
const double cylinderLength = 15.0;
static double horizon;
const PdImp::BulkModulus _K(130000.0);
const PdImp::PoissonsRatio _MU(0.0);
const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
const double g = 9.807e-3;
const double rho = 7800e-6;
const int vectorNDF=3;


TensorProductSolidCylinder getMeshGenerator(){
	TensorProductSolidCylinder meshGen(numProcs,cylinderRadius,numRings,zStart,cylinderLength);
	const std::vector<PdQPointSet1d>& specs = meshGen.getTensorProductSpecs();
	PdQPointSet1d zSpec = specs[2];
	double dZ = zSpec.getCellSize();
	/*
	 * Assign a couple of static variables
	 */
	zMax = (zStart+cylinderLength)-dZ/2.0;
	zMin = zStart + dZ/2.0;
	horizon = 1.5*sqrt(3.0)*dZ;
	return meshGen;
}

void utPimp_pullCylinder() {
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	TensorProductSolidCylinder meshGen = getMeshGenerator();
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, meshGen);
	decomp = getLoadBalancedDiscretization(decomp);
//	decomp = createAndAddNeighborhood(decomp,horizon);

	/*
	 * Create Pimp Operator
	 */
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);
	PdITI::PdITI_Operator op(comm,list,decomp.cellVolume);
//	shared_ptr<PdImp::PdImpOperator> op = utPdITI::getPimpOperator(decomp,comm);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);


	/*
	 * Get points for bc's
	 */
	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::Z;
	Pd_shared_ptr_Array<int> bcIdsFixed = PdNeighborhood::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,decomp.myX,decomp.numPoints,horizon,zMin);
	std::sort(bcIdsFixed.get(),bcIdsFixed.get()+bcIdsFixed.getSize());
	Pd_shared_ptr_Array<int> bcIdsApplied = PdNeighborhood::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,decomp.myX,decomp.numPoints,horizon,zMax);
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
	c[0] = DirichletBcSpec::Z;
	ComponentDirichletBcSpec appliedSpec(c,bcIdsApplied);
	StageFunction dispStageFunction(1.0e-3,1.0e-3);
	shared_ptr<StageComponentDirichletBc> bcApplied(new StageComponentDirichletBc(appliedSpec,dispStageFunction));
	bcs[1] = bcApplied;

	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	Field<double> uOwnedField(DISPL3D,decomp.numPoints);
	uOwnedField.setValue(0.0);
	bcApplied->applyKinematics(1.0,uOwnedField);
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getJacobian(uOwnedField);

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
	Field_NS::Field<double> fN(fNSpec,decomp.numPoints);
	fN.setValue(0.0);
	op.computeInternalForce(uOwnedField,fN);
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
	Field_NS::Field<double> delta(deltaSpec,decomp.numPoints);
	delta.setValue(0.0);
	bcApplied->applyKinematics(1.0,delta);

	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,decomp.numPoints);
	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,delta);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utPimp_pullCylinder_npX.pvtu", comm.NumProc(), comm.MyPID());
	PdVTK::write(writer,grid);

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPimp_pullCylinder" );
	proc->add(BOOST_TEST_CASE( &utPimp_pullCylinder ));
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
