/*
 * utPdIti_twoPointPlasticity
 *
 *  Created on: Sept 11, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include "PdQuickGrid.h"
#include <math.h>
#include "PdMaterialUtilities.h"
#include "Field.h"
#include "PdVTK.h"
#include "../../PdImpMpiFixture.h"
#include "../../PdImpMaterials.h"
#include "../../PdImpOperator.h"
#include "../../IsotropicElasticPlasticModel.h"
#include "../../ConstitutiveModel.h"
#include "../../DirichletBcSpec.h"
#include "../../StageFunction.h"
#include "../../StageComponentDirichletBc.h"
#include "../../ComponentDirichletBcSpec.h"
#include <fstream>

using namespace boost::unit_test;
using namespace std;
using namespace PdMaterialUtilities;
using namespace Field_NS;
using PdImp::ConstitutiveModel;
using PdImp::IsotropicElasticPlasticModel;
using std::tr1::shared_ptr;
using namespace PdImp;
using PdImp::DirichletBcSpec;
using PdImp::ComponentDirichletBcSpec;
using PdImp::StageFunction;
using PdImp::StageComponentDirichletBc;
using namespace PdVTK;


static int myRank;
static int numProcs;


static double horizon=sqrt(2);

/*
 * Young's Modulus (MPa)
 */
static double E = 68.9e3;

/*
 * Poisson's ratio
 */
static double nu = .33;

/*
 * Yield "Stress" estimate for perfect plasticity (MPa)
 * 6061-T6 data
 */
static double Y = 351.79;

/*
 * yield strain ~.0051 -- 1/2 the engineering strain -- not the same as the
 * shear strain used here to load
 */
static double epsYield = .0051;

/*
 * Density of aluminum g/mm^3
 */
static double rho = 2.7e-3;

/*
 * Bulk Modulus
 */
static double K = E / 3 / (1-2.0 * nu);

/*
 * Shear Modulus
 */
static double mu = E / 2 / (1+nu);

IsotropicElasticPlasticSpec getMaterialSpec() {
	YieldStress yieldStress = IsotropicElasticPlasticSpec::yieldStress(Y);
	MaterialHorizon materialHorizon = IsotropicElasticPlasticSpec::materialHorizon(horizon);
	YoungsModulus youngsModulus = IsotropicHookeSpec::youngsModulus(E);
	PoissonsRatio poissonsRatio = IsotropicHookeSpec::poissonsRatio(nu);
	IsotropicElasticPlasticSpec elasticPlasticSpec(yieldStress,materialHorizon,IsotropicHookeSpec(youngsModulus,poissonsRatio));
	return elasticPlasticSpec;
}

PdGridData getTwoPointGridData(){
	int numCells = 2;
	int dimension = 3;
	PdGridData pdGridData = PdQuickGrid::allocatePdGridData(numCells,dimension);
	/*
	 * Create points
	 */
	double w=1.0;
	pdGridData.globalNumPoints=2;
	pdGridData.numPoints=2;
	pdGridData.numExport=0;
	double *x = pdGridData.myX.get();
	*(x+0)=0.0;
	*(x+1)=0.0;
	*(x+2)=0.0;
	*(x+3)=w;
	*(x+4)=w;
	*(x+5)=0.0;

	/*
	 * Global ids
	 */
	 int *ids = pdGridData.myGlobalIDs.get();
	 *(ids+0)=0;
	 *(ids+1)=1;

	 /*
	  * Cell volumes
	  */
	 double *v = pdGridData.cellVolume.get();
	 *(v+0)=1.0;
	 *(v+1)=1.0;

	 /*
	  * Create neighborhood
	  */
	 pdGridData.sizeNeighborhoodList=4;
	 shared_ptr<int> neighborhood = shared_ptr<int>(new int[pdGridData.sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	 pdGridData.neighborhood = neighborhood;
	 int *neigh = neighborhood.get();
	 *(neigh+0)=1;
	 *(neigh+1)=1;
	 *(neigh+2)=1;
	 *(neigh+3)=0;

	 shared_ptr<int> neighborhoodPtr = shared_ptr<int>(new int[pdGridData.numPoints],PdQuickGrid::Deleter<int>());
	 pdGridData.neighborhoodPtr=neighborhoodPtr;
	 int *nPtr = neighborhoodPtr.get();
	 *(nPtr+0)=0;
	 *(nPtr+1)=2;

	 return pdGridData;
}

shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm) {
	PdImp::PdImpOperator *op = new PdImp::PdImpOperator(comm,decomp);
	return shared_ptr<PdImp::PdImpOperator>(op);
}


void runPureShear() {
	IsotropicElasticPlasticSpec matSpec = getMaterialSpec();
	PdGridData pdGridData = getTwoPointGridData();
	int numPoints = pdGridData.numPoints;
	BOOST_CHECK(2 == numPoints);
	BOOST_CHECK(4 == pdGridData.sizeNeighborhoodList);

	/*
	 * Communicator
	 */
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);

	/*
	 * Create PdImp Operator
	 */
	shared_ptr<PdImp::PdImpOperator> op = getPimpOperator(pdGridData,comm);
	shared_ptr<ConstitutiveModel> fIntOperator(shared_ptr<ConstitutiveModel>(new IsotropicElasticPlasticModel(matSpec)));
	op->addConstitutiveModel(fIntOperator);

	/*
	 * Point '0' is fixed
	 * Point '1' has applied displacement
	 */
	Pd_shared_ptr_Array<int> bcIdsFixed(1);
	Pd_shared_ptr_Array<int> bcIdsApplied(1);
	{
		int *id0 = bcIdsFixed.get();
		*id0 = 0;
		int *id1 = bcIdsApplied.get();
		*id1 = 1;
	}
	/*
	 * Dirichlet BC specs
	 */
	// Fixed point
	ComponentDirichletBcSpec fixedSpec = ComponentDirichletBcSpec::getAllComponents(bcIdsFixed);
	StageFunction constStageFunction;
	shared_ptr<StageComponentDirichletBc> bcFixed(new StageComponentDirichletBc(fixedSpec,constStageFunction));

	/*
	 * Constrain y,z components at applied end
	 */
	vector<DirichletBcSpec::ComponentLabel> yzLabels(2);
	yzLabels[0] = DirichletBcSpec::Y;
	yzLabels[1] = DirichletBcSpec::Z;
	ComponentDirichletBcSpec yzFixedSpec(yzLabels,bcIdsApplied);
	shared_ptr<StageComponentDirichletBc> bcYZ_Fixed(new StageComponentDirichletBc(yzFixedSpec,constStageFunction));

	/*
	 * Applied displacement
	 */
	vector<DirichletBcSpec::ComponentLabel> xLabel(1);
	xLabel[0] = DirichletBcSpec::X;
	ComponentDirichletBcSpec xAppliedSpec(xLabel,bcIdsApplied);

	/*
	 * Displacement and Internal Force Vectors
	 */
	FieldSpec uSpec(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR3D,"u");
	Field<double> uOwnedField(uSpec,pdGridData.numPoints);

	FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fNField(fNSpec,pdGridData.numPoints);
	const FieldSpec velocitySpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D, "velocity");
	Field<double> velField(velocitySpec,pdGridData.numPoints);
	uOwnedField.setValue(0.0);
	velField.setValue(0.0);
	fNField.setValue(0.0);
	double *u1x = uOwnedField.getArray().get()+3;
	double *v1x = velField.getArray().get()+3;
	double *f1x = fNField.getArray().get()+3;

	/*
	 * Stages
	 * 1) load
	 * 2) unload
	 * 3) reload
	 */
	std::vector<StageFunction> stages(3);
	/*
	 * Loading
	 */
	stages[0] = StageFunction(0,epsYield);
	/*
	 * Unloading
	 */
	stages[1] = stages[0].next(-.001275);

	/*
	 * Re-Unloading
	 */
	stages[2] = stages[1].next(.001275);

	int numStepsPerStage = 50;
	double dt = 1.0/numStepsPerStage;

	/*
	 * Create data file
	 */
	std::fstream out("ep.dat", std::fstream::out);
	/*
	 * Write out initial condition
	 */
	double t=0;
	out << 0 << " " << 0 << " " << 0 << std::endl;
	for(std::vector<StageFunction>::iterator stageIter=stages.begin(); stageIter!=stages.end();stageIter++){

		*v1x = stageIter->slope();

		for(int step=0;step<numStepsPerStage;step++){

			*u1x += *v1x * dt;
			op->computeInternalForce(uOwnedField,fNField,false);

			/*
			 * Get sign of "f" -- this works as long as f does not ever land "exactly" on zero
			 * Put a negative sign in front so that loading is "positive"
			 */
			double signF = -*f1x/abs(*f1x);
			t += dt;
			out << t << " " << *u1x << " " << signF*MAGNITUDE(f1x) << std::endl;
			op->advanceStateVariables();


		}
	}
	out.close();
}


bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPdIti_twoPointPlasticity");
  proc->add(BOOST_TEST_CASE(&runPureShear));
  framework::master_test_suite().add(proc);

  return success;
}

bool init_unit_test()
{
  init_unit_test_suite();
  return true;
}

int main
(int argc, char* argv[])
{
	// Initialize MPI and timer
	PdImpRunTime::PimpMpiFixture pimpMPI = PdImpRunTime::PimpMpiFixture::getPimpMPI(argc,argv);
	const Epetra_Comm& comm = pimpMPI.getEpetra_Comm();

	// These are static (file scope) variables
	myRank = comm.MyPID();
	numProcs = comm.NumProc();

  // Initialize UTF
  return unit_test_main(init_unit_test, argc, argv);
}
