/*
 * utPdImp_eigen1d.cxx
 *
 *  Created on: Aug 23, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "../../PdImpMpiFixture.h"
#include "PdNeighborhood.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdNeighborhood.h"
#include "PdZoltan.h"
#include "Field.h"
#include "../../PdImpMaterials.h"
#include "../../PdImpOperator.h"
#include "../../PdImpOperatorUtilities.h"
#include "../../DirichletBcSpec.h"
#include "../../BodyLoadSpec.h"
#include "../../StageFunction.h"
#include "../../Loader.h"
#include "../../StageComponentDirichletBc.h"
#include "../../ComponentDirichletBcSpec.h"
#include "../../IsotropicElasticConstitutiveModel.h"
#include "../../ConstitutiveModel.h"
#include "../../PdImpMaterials.h"
#include "../../ImplicitLinearDynamicsIntegrator.h"
#include "../../NewmarkBetaIntegrator.h"
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
#include <iostream>
#include <fstream>
using namespace std;



using namespace PdQuickGrid;
using namespace PdImp;
using namespace PdNeighborhood;
using namespace Field_NS;
using PdITI::IsotropicElasticConstitutiveModel;
using PdITI::ConstitutiveModel;
using std::tr1::shared_ptr;
using namespace boost::unit_test;

static int numProcs;
static int myRank;

const double cylinderLength = 15.0;
const int nx = 1;
const int ny = nx;
const int nz = 150;
const double lX = 1.5;
const double lY = lX;
const double lZ = cylinderLength;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const double zStart  =  0.0;
const double zLength =  lZ;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const double horizon = 3.05 * zLength/nz;
const int vectorNDF=3;


/*
 * Young's Modulus (MPa) -- standard for aluminum
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
 * Initial velocity of point 1
 */
static double INITIAL_VELOCITY=1.0;

IsotropicHookeSpec getMaterialSpec() {

	/*
	 * 1D material spec requires and adjustment to properties
	 */
	BulkModulus K = IsotropicHookeSpec::bulkModulus(E/9.0);
	PoissonsRatio poissonsRatio = IsotropicHookeSpec::poissonsRatio(nu);
	return IsotropicHookeSpec(K,poissonsRatio);
}

PdGridData getGridData(){

	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData pdGridData =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	pdGridData = getLoadBalancedDiscretization(pdGridData);
	return pdGridData;

}

shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm) {
	PdImp::PdImpOperator *op = new PdImp::PdImpOperator(comm,decomp);
	return shared_ptr<PdImp::PdImpOperator>(op);
}

void initialConditions
(
		TemporalField<double>& displacement,
		TemporalField<double>& velocity,
		TemporalField<double>& acceleration,
		PdGridData& pdGridData
)
{


	/*
	 * Initialize everything to 0
	 */
	double zero(0.0);
	FieldSpec::FieldStep NP1 = FieldSpec::STEP_NP1;
	FieldSpec::FieldStep N = FieldSpec::STEP_N;
	PdITI::SET(displacement.getField(N).getArray().get(),displacement.getField(N).getArray().end(),zero);
	PdITI::SET(displacement.getField(NP1).getArray().get(),displacement.getField(NP1).getArray().end(),zero);
	PdITI::SET(velocity.getField(N).getArray().get(),velocity.getField(N).getArray().end(),zero);
	PdITI::SET(velocity.getField(NP1).getArray().get(),velocity.getField(NP1).getArray().end(),zero);
	PdITI::SET(acceleration.getField(N).getArray().get(),acceleration.getField(N).getArray().end(),zero);
	PdITI::SET(acceleration.getField(NP1).getArray().get(),acceleration.getField(NP1).getArray().end(),zero);

	/*
	 * Assign initial velocity to first eigen function
	 */
	int numPoints = pdGridData.numPoints;
	double *zPtr = pdGridData.myX.get()+2;
	double *vZPtr = velocity.getField(NP1).getArray().get()+2;
	for(int p=0;p<numPoints;p++,zPtr+=3,vZPtr+=3){
		*vZPtr = INITIAL_VELOCITY * cos(M_PI*(*zPtr)/cylinderLength);
	}

}

shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian){
	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();
	const Epetra_BlockMap& colMap = jacobian->getColMap();

	/*
	 * Epetra Graph
	 */
	Pd_shared_ptr_Array<int> numCols = jacobian->getNumColumnsPerRow();
	Epetra_CrsGraph *graph = new Epetra_CrsGraph(Copy,rowMap,colMap,numCols.get());
	for(int row=0;row<rowMap.NumMyElements();row++){
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.getSize();
		int *cols = rowLIDs.get();
		BOOST_CHECK(0==graph->InsertMyIndices(row,numCol,cols));
	}
	BOOST_CHECK(0==graph->FillComplete());
	return shared_ptr<Epetra_CrsGraph>(graph);
}

shared_ptr<Epetra_RowMatrix>
getOperator
(
		const vector<shared_ptr<StageComponentDirichletBc> >& bcArray,
		shared_ptr<Epetra_CrsGraph>& graphPtr,
		shared_ptr<RowStiffnessOperator>& jacobian,
		ImplicitLinearDynamicsIntegrator implicitIntegrator,
		double dt
)
{
	std::cout << "Begin jacobian calculation\n";
	double beta = implicitIntegrator.getNewmarkIntegrator().getBeta();
	double density = implicitIntegrator.getDensity().getValue();


	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();

	/*
	 * Loop over Bc's and create set of ids that can be searched
	 */

	vector<std::set<int> > bcPointIds(bcArray.size());
	{
		vector<shared_ptr<StageComponentDirichletBc> >::const_iterator bcIter = bcArray.begin();
		for(int i=0;bcIter != bcArray.end(); i++,bcIter++){
			StageComponentDirichletBc* stageComponentPtr = bcIter->get();
			const DirichletBcSpec& spec = stageComponentPtr->getSpec();
			const Pd_shared_ptr_Array<int>& ids = spec.getPointIds();
			bcPointIds[i]= std::set<int>(ids.get(),ids.get()+ids.getSize());
		}
	}

	/*
	 * Epetra Matrix
	 * PERHAPS the 'operator' can keep its own copy of Graph, then this
	 * can be a 'View'
	 */

	Epetra_FEVbrMatrix *m = new Epetra_FEVbrMatrix(Copy,*(graphPtr.get()));

	Epetra_SerialDenseMatrix k;
	k.Shape(vectorNDF,vectorNDF);
	for(int row=0;row<rowMap.NumMyElements();row++){
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.getSize();
		int *cols = rowLIDs.get();
		BOOST_CHECK(0==m->BeginReplaceMyValues(row,numCol,cols));

		/*
		 * loop over columns in row and submit block entry
		 */
		Pd_shared_ptr_Array<double> actualK = jacobian->computeRowStiffness(row, rowLIDs);

		/**
		 * 1) First, loop through computed stiffness and scale for implicit dynamics
		 * 2) Also add in density on diagonal
		 */
		{
			/*
			 * These are pointers to the diagonal
			 */
			double *xPtr = actualK.get();
			double *yPtr = actualK.get() + 4;
			double *zPtr = actualK.get() + 8;
			double *kPtr = actualK.get();
			for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){

				/*
				 * Scale Jacobian by -1
				 */
				PdITI::SCALE_BY_VALUE(kPtr,kPtr+9,-1);

				int col = *colPtr;
				double c = density/(beta*dt*dt);

				/*
				 * If this is true, then we are on the diagonal
				 * Add in density
				 */
				if(row==col) {
					*xPtr += c;
					*yPtr += c;
					*zPtr += c;
				}
				/*
				 * Increment pointers to diagonal of 3x3 matrix and also pointer to start of 3x3 matrix
				 */
				xPtr+=9;
				yPtr+=9;
				zPtr+=9;
				kPtr+=9;

			}
		}


		/*
		 * 1) Zero out row as necessary due to BC's
		 * 2) Assemble into matrix -- zero columns as necessary
		 * 3)
		 */
		vector<std::set<int> >::iterator pointSetIter = bcPointIds.begin();
		vector<shared_ptr<StageComponentDirichletBc> >::const_iterator bcIter = bcArray.begin();
		for(;bcIter != bcArray.end(); bcIter++, pointSetIter++){
			const std::set<int>::const_iterator bcIdsEnd = pointSetIter->end();

			/*
			 * Get components to be applied
			 */
			StageComponentDirichletBc* stageComponentPtr = bcIter->get();
			const DirichletBcSpec& spec = stageComponentPtr->getSpec();
			vector<DirichletBcSpec::ComponentLabel> components = spec.getComponents();

			/*
			 * Search for row in bcIds
			 */
			// if this is true, then this row is a bc row
			if(bcIdsEnd != pointSetIter->find(row)) {

				/*
				 * 1) This row has a bc; need to zero out in stiffness
				 * 2) Watch for diagonal: place "1" on the diagonal
				 */

				/*
				 * Create array of pointers to each row/component that BC is applied
				 */
				vector<double*> rowPtrs(components.size(), NULL);
				vector<double*> diagPtrs(components.size(), NULL);
				for(std::size_t r=0;r<components.size();r++)
					rowPtrs[r]= actualK.get()+components[r];

				for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){
					int col = *colPtr;

					for(std::size_t r=0;r<components.size();r++){
						/*
						 * 0) Save diagonal location
						 * 1) Set row element to zero
						 * 2) Move to next column
						 * 3) Note that there are 3 columns in 3x3 matrix
						 */

						/*
						 * Save location of diagonal for this component
						 */
						diagPtrs[r] = rowPtrs[r] + 3 * components[r];

						/*
						 * Zero out row corresponding with component
						 * Increment row pointer to next column
						 */
						*(rowPtrs[r])=0; rowPtrs[r]+=3;
						*(rowPtrs[r])=0; rowPtrs[r]+=3;
						*(rowPtrs[r])=0; rowPtrs[r]+=3;

					}
					/*
					 * If this is true, then we are on the diagonal
					 * Need to put "1" on the strict diagonal by component
					 */
					if(row==col) {
						for(std::size_t r=0;r<components.size();r++){
							*(diagPtrs[r]) = 1.0;
						}

					}

				}

			} else {
				/*
				 * This is not a bc row but we still have to check for columns that may have a bc assigned to them
				 * In this case, we just have to zero the column
				 */
				double *kPtr=actualK.get();
				for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++,kPtr+=9){
					int col = *colPtr;
					if(bcIdsEnd == pointSetIter->find(col)) continue;
					/*
					 * We have a column that must be zero'd
					 */
					/*
					 * Zero out column for each component that is applied
					 */
					double *colPtr=0;
					for(std::size_t r=0;r<components.size();r++){
						colPtr = kPtr+3*components[r];
						for(int r=0;r<3;r++)
							colPtr[r]=0;
					}

				}

			}

		}

		/*
		 * Now just populate the matrix
		 */
		double *kPtr = actualK.get();
		for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){

			/*
			 * Fill 'k'
			 */
			for(int c=0;c<3;c++){
				double *colK = k[c];
				for(int r=0;r<3;r++,kPtr++)
					colK[r] = *kPtr;
			}

			BOOST_CHECK(0==m->SubmitBlockEntry(k));
		}


		/*
		 * Finalize this row
		 */
		BOOST_CHECK(0==m->EndSubmitEntries());
	}

	BOOST_CHECK(0==m->FillComplete());

	std::cout << "END jacobian calculation\n";
	return shared_ptr<Epetra_RowMatrix>(m);
}


void runWave() {
	IsotropicHookeSpec matSpec = getMaterialSpec();
	PdGridData pdGridData =  getGridData();

	/*
	 * Communicator
	 */
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);

	/*
	 * Create PdImp Operator
	 */
	shared_ptr<PdImp::PdImpOperator> op = getPimpOperator(pdGridData,comm);
	shared_ptr<PdITI::ConstitutiveModel> fIntOperator(shared_ptr<ConstitutiveModel>(new PdITI::IsotropicElasticConstitutiveModel(matSpec)));
	op->addConstitutiveModel(fIntOperator);

	/*
	 * Displacement and Internal Force Vectors
	 */
	const FieldSpec displacementSpec(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR3D, "Displacement");
	const FieldSpec velocitySpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D, "v");
	const FieldSpec accelerationSpec(FieldSpec::ACCELERATION,FieldSpec::VECTOR3D, "a");
	const FieldSpec residualSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D, "residual");
	TemporalField<double> ut(displacementSpec,pdGridData.numPoints);
	TemporalField<double> vt(velocitySpec,pdGridData.numPoints);
	TemporalField<double> at(accelerationSpec,pdGridData.numPoints);
	TemporalField<double> rt(residualSpec,pdGridData.numPoints);

	// Initialize above vectors to ZERO and also set initial velocities on left and right end of cylinder
	initialConditions(ut,vt,at,pdGridData);

	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op->getRowStiffnessOperator(ut.getField(FieldSpec::STEP_N),horizon);

	/*
	 * Create graph
	 */
	shared_ptr<Epetra_CrsGraph> graphPtr = getGraph(jacobian);

	/*
	 * Weighted volume
	 */
	Field_NS::Field<double> mField = op->getWeightedVolume();

	/*
	 * Create grid -- since this is a linear case, current coordinates are always initial coordinates;
	 * Set up grid with coordinates x
	 * Write fields for NP1
	 * Note that each time step NP1 is updated but grid and writer always uses the same pointer; Because
	 * the state of u, and v are swapped at the end of each step, it is necessary to re-assign these pointers
	 * inside the time loop.
	 */


	FieldSpec::FieldStep NP1 = FieldSpec::STEP_NP1;
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(pdGridData.myX,pdGridData.numPoints);
	PdVTK::writeField(grid,ut.getField(NP1));
	PdVTK::writeField(grid,vt.getField(NP1));
	PdVTK::writeField(grid,mField);
	PdVTK::CollectionWriter collectionWriter("utPdImp_eigen1D",comm.NumProc(), comm.MyPID(), PdVTK::vtkASCII);


	/*
	 * Initial conditions are in NP1
	 * Write this state
	 * Advance state
	 * Loop on time {
	 *    Integrate forward
	 *    Write state NP1
	 *    Advance state
	 * }
	 */
	double t = 0.0;
	MassDensity density(rho);
	ImplicitLinearDynamicsIntegrator integrator(density);
	const NewmarkBetaIntegrator& newmark = integrator.getNewmarkIntegrator();
	collectionWriter.writeTimeStep(t,grid);
	newmark.advanceState(ut,vt,at);

	/*
	 * N -- number of steps per period
	 */
	std::size_t numStepsPerPeriod = 48;
	/*
	 * nT -- eigen function number
	 */
	std::size_t numPeriods = 1;
	int periodNumber = 1;
	double omega = (periodNumber*M_PI)*sqrt(E/rho)/(cylinderLength);
	double period = 2.0 * M_PI / omega;
	std::size_t numSteps = numStepsPerPeriod*numPeriods;
	double dt = period/numStepsPerPeriod;

	/**
	 * Create boundary conditions spec -- THERE ARE NONE
	 * Create Epetra_RowMatrix
	 */
	vector<shared_ptr<StageComponentDirichletBc> > bcs;
	shared_ptr<Epetra_RowMatrix> mPtr = getOperator(bcs,graphPtr,jacobian,integrator,dt);

	for(std::size_t step=0;step<numSteps;step++){
		/*
		 * advance time step
		 */
		t += dt;

		/*
		 * Compute residual
		 */
		integrator.computeResidual(ut,vt,at,dt,rt);

		/*
		 * Solve linear problem: computes uNp1
		 */
		Epetra_LinearProblem linProblem;
		linProblem.SetOperator(mPtr.get());
		linProblem.AssertSymmetric();

		const Epetra_BlockMap& rangeMap  = mPtr->OperatorRangeMap();
		Epetra_Vector lhs(View,rangeMap,ut.getField(FieldSpec::STEP_NP1).getArray().get());
		Epetra_Vector rhs(View,rangeMap,rt.getField(FieldSpec::STEP_NP1).getArray().get());
		linProblem.SetLHS(&lhs);
		linProblem.SetRHS(&rhs);

		AztecOO solver(linProblem);
		solver.SetAztecOption(AZ_precond, AZ_Jacobi);
		BOOST_CHECK(0==solver.CheckInput());
		solver.Iterate(500,1e-6);

		/*
		 * Integrate step: compute vNp1, aNp1
		 */
		newmark.integrateStep(ut,vt,at,dt);

		/*
		 * Write uNp1, vNp1
		 * NOTE: We have to re-write these each time since the pointer underlying the data is swapped each time step
		 */
		PdVTK::writeField(grid,ut.getField(NP1));
		PdVTK::writeField(grid,vt.getField(NP1));
		collectionWriter.writeTimeStep(t,grid);

		/*
		 * Now advance step and prepare for next time step
		 */
		newmark.advanceState(ut,vt,at);

	}

	/*
	 * This writes the "pvd" collection file
	 */
	collectionWriter.close();
}


bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPdImp_eigen1D");
  proc->add(BOOST_TEST_CASE(&runWave));
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
