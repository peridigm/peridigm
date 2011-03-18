/*
 * utPimp_smallCylinderPull_npX.cxx
 *
 *  Created on: Jul 1, 2010
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
#include "utPdITI.h"
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

using namespace Pdut;
using UTILITIES::CartesianComponent;
using namespace Field_NS;
using UTILITIES::Array;
using PdITI::IsotropicElasticConstitutiveModel;
using PdITI::ConstitutiveModel;
using PdImp::StageComponentDirichletBc;
using PdImp::StageFunction;
using PdImp::ComponentDirichletBcSpec;
using PdImp::DirichletBcSpec;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::vector;
using std::cout;
using std::endl;

static size_t numProcs;
static size_t myRank;

const double cylinderRadius = 1.0;
const int numRings = 2;
const double zStart = 0.0;
const double cylinderLength = 1.0;
const double horizon = .6;
static double zMax=0.0;
static double zMin=0.0;

const PdImp::BulkModulus _K(130000.0);
const PdImp::PoissonsRatio _MU(0.0);
const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);

//shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm);

QUICKGRID::TensorProductSolidCylinder getMeshGenerator(){

	QUICKGRID::TensorProductSolidCylinder meshGen(numProcs,cylinderRadius,numRings,zStart,cylinderLength);
	return meshGen;
}

QUICKGRID::QuickGridData getGrid() {
	QUICKGRID::TensorProductSolidCylinder meshGen = getMeshGenerator();
	QUICKGRID::QuickGridData decomp =  QUICKGRID::getDiscretization(myRank, meshGen);
	decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

	const std::vector<QUICKGRID::Spec1D>& specs = meshGen.getTensorProductSpecs();
	QUICKGRID::Spec1D zSpec = specs[2];
	double dZ = zSpec.getCellSize();
	/*
	 * Assign a couple of static variables
	 */
	zMax = (zStart+cylinderLength)-dZ/2.0;
	zMin = zStart + dZ/2.0;

	return decomp;
}

void computeInternalForce(){
	QUICKGRID::QuickGridData decomp = getGrid();

	/*
	 * Get points for bc's
	 */
	CartesianComponent axis = UTILITIES::Z;
	Array<int> bcIdsFixed = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,decomp.myX,decomp.numPoints,horizon,zMax);
	std::sort(bcIdsFixed.get(),bcIdsFixed.end());
	Array<int> bcIdsApplied = UTILITIES::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,decomp.myX,decomp.numPoints,horizon,zMin);
	std::sort(bcIdsApplied.get(),bcIdsApplied.end());

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
	StageFunction dispStageFunction(-1.0e-3,-1.0e-3);
	shared_ptr<StageComponentDirichletBc> bcApplied(new StageComponentDirichletBc(appliedSpec,dispStageFunction));
	bcs[1] = bcApplied;

	Field<double> uOwnedField(DISPL3D,decomp.numPoints);
	uOwnedField.set(0.0);
	bcApplied->applyKinematics(1.0,uOwnedField);


	/*
	 * Create PdITI Operator
	 */
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);
	PdITI::PdITI_Operator op(comm,list,decomp.cellVolume);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);

	FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field<double> fN(fNSpec,decomp.numPoints);
	fN.set(0.0);
	op.computeInternalForce(uOwnedField,fN);

	/*
	 * Write problem set up parameters to file
	 */
	int numPoints = decomp.numPoints;
	Field_NS::Field<double> volField = Field_NS::getVOLUME(decomp.cellVolume,numPoints);
	Field<int> fieldRank(Field_NS::PROC_NUM,decomp.numPoints);
	fieldRank.set(myRank);
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);
	PdVTK::writeField<double>(grid,volField);
	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,fN);
	PdVTK::writeField(grid,fieldRank);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utPimp_smallCylinderPull_npX.pvtu", numProcs, myRank);
	PdVTK::write(writer,grid);

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPimp_smallCylinderPull_npX" );
	proc->add(BOOST_TEST_CASE( &computeInternalForce ));
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
