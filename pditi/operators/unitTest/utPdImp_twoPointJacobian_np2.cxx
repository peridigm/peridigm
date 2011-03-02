/*
 * utPimp_twoPointJacobian_np2.cxx
 *
 *  Created on: Jun 9, 2010
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
#include "../PdImpMaterials.h"
#include "../PdImpOperator.h"
#include "../PdImpOperatorUtilities.h"
#include "../IsotropicElasticConstitutiveModel.h"
#include "../StageComponentDirichletBc.h"
#include "../ComponentDirichletBcSpec.h"
#include "../DirichletBcSpec.h"
#include "../StageFunction.h"
#include <set>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Vector.h>
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
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;
using std::endl;
using std::vector;
using std::set;
using namespace PdImp;

static int myRank;
static int numProcs;

const int nx = 2;
const int ny = 1;
const int nz = 1;
const double lX = 1.0;
const double theta = 0.0;
const double xStart  = -lX/4.0;
const double xLength =  lX;
const double xMin = xStart+lX/nx/2.0;
const double xMax = xStart+xLength-lX/nx/2.0;
const double yStart  = -lX/ny/2.0;
const double yLength =  lX;
const double zStart  = -lX/nz/2.0;
const double zLength =  lX;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=1.1*lX;
const double delta = 1.0e-2 * lX;
const double bond = lX/2.0;
const PdImp::BulkModulus _K(130000.0);
const PdImp::PoissonsRatio _MU(0.25);
const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
const int vectorNDF=3;

void setRotatedValue(double* x, double magnitude);
shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian);
shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm);
shared_ptr<Epetra_VbrMatrix>
getOperator
(
		const vector<shared_ptr<StageComponentDirichletBc> >& bcArray,
		shared_ptr<Epetra_CrsGraph>& graphPtr,
		shared_ptr<RowStiffnessOperator>& jacobian
);


PdGridData getGrid() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	decomp = getLoadBalancedDiscretization(decomp);
	BOOST_CHECK(1==decomp.numPoints);
	BOOST_CHECK(2==decomp.globalNumPoints);
	int myPoint = *(decomp.myGlobalIDs.get());

	if(myPoint==1){
		/*
		 * rotate point 1 by theta degrees
		 */
		double *x1 = decomp.myX.get();
		setRotatedValue(x1,lX);
	}

	return decomp;
}

void computeJacobian(){
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	PdGridData decomp = getGrid();
	decomp = getLoadBalancedDiscretization(decomp);

	/*
	 * Create Pimp Operator
	 */
	shared_ptr<PdImp::PdImpOperator> op = getPimpOperator(decomp,comm);

	/*
	 * Get points for bc's
	 */
	PdNeighborhood::CoordinateLabel axis = PdNeighborhood::X;
	Pd_shared_ptr_Array<int> bcIdsFixed = PdNeighborhood::getPointsInNeighborhoodOfAxisAlignedMinimumValue(axis,decomp.myX,decomp.numPoints,1.0e-3*horizon,xMin);
	std::sort(bcIdsFixed.get(),bcIdsFixed.get()+bcIdsFixed.getSize());
	Pd_shared_ptr_Array<int> bcIdsApplied = PdNeighborhood::getPointsInNeighborhoodOfAxisAlignedMaximumValue(axis,decomp.myX,decomp.numPoints,1.0e-3*horizon,xMax);
	std::sort(bcIdsApplied.get(),bcIdsApplied.get()+bcIdsApplied.getSize());

	/**
	 * Create boundary conditions spec
	 */
	vector<shared_ptr<StageComponentDirichletBc> > bcs(3);
	ComponentDirichletBcSpec fixedSpec1 = ComponentDirichletBcSpec::getAllComponents(bcIdsFixed);
	StageFunction constStageFunction(0.0,0.0);
	shared_ptr<StageComponentDirichletBc> bcFixed1(new StageComponentDirichletBc(fixedSpec1,constStageFunction));
	bcs[0] = bcFixed1;
	std::vector< DirichletBcSpec::ComponentLabel > c2(2);
	c2[0] = DirichletBcSpec::Y;
	c2[1] = DirichletBcSpec::Z;
	ComponentDirichletBcSpec fixedSpec2(c2,bcIdsApplied);
	shared_ptr<StageComponentDirichletBc> bcFixed2(new StageComponentDirichletBc(fixedSpec2,constStageFunction));
	bcs[1] = bcFixed2;

	std::vector< DirichletBcSpec::ComponentLabel > c3(1);
	c3[0] = DirichletBcSpec::X;
	ComponentDirichletBcSpec appliedSpec(c3,bcIdsApplied);
	StageFunction dispStageFunction(1.0e-3,1.0e-3);
	shared_ptr<StageComponentDirichletBc> bcApplied(new StageComponentDirichletBc(appliedSpec,dispStageFunction));
	bcs[2] = bcApplied;



	/*
	 * Create Jacobian -- note that SCOPE of jacobian is associated with the PimpOperator "op"
	 */
	Field<double> uOwnedField(DISPL3D,decomp.numPoints);
	uOwnedField.setValue(0.0);
	bcApplied->applyKinematics(1.0,uOwnedField);
	/*
	 * NOTE: once the jacobian has been computed with the given displacement field, it is safe to
	 * change the displacement field -- BUT NOT BEFORE jacobian has been computed otherwise since
	 * the jacobian is a linearization about the displacement field
	 */
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op->getRowStiffnessOperator(uOwnedField,horizon);


	/*
	 * Create graph
	 */
	shared_ptr<Epetra_CrsGraph> graphPtr = getGraph(jacobian);

	/*
	 * Create Epetra_RowMatrix
	 */
	shared_ptr<Epetra_VbrMatrix> mPtr = getOperator(bcs,graphPtr,jacobian);

	/*
	 * Create force field
	 * IN this case,
	 * 1) Compute internal force with displacement vector that has kinematics applied
	 * 2) Negate internal force since we are using it as a residual on the RHS
	 * 3) Apply kinematics to this vector so that solution properly includes the applied kinematics
	 */
	Field_NS::FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fN(fNSpec,decomp.numPoints);
	fN.setValue(0.0);

	op->computeInternalForce(uOwnedField,fN);
	PdITI::SCALE_BY_VALUE(fN.getArray().get(),fN.getArray().end(),-1.0);
	bcApplied->applyKinematics(1.0,fN);
	bcFixed1->applyHomogeneousForm(fN);
	bcFixed2->applyHomogeneousForm(fN);
	uOwnedField.setValue(0.0);
//	bcApplied->applyKinematics(1.0,uOwnedField);


	Epetra_LinearProblem linProblem;
	linProblem.SetOperator(mPtr.get());
	linProblem.AssertSymmetric();

	const Epetra_BlockMap& rangeMap  = mPtr->OperatorRangeMap();
	const Epetra_BlockMap& domainMap  = mPtr->OperatorDomainMap();
	const Epetra_BlockMap& rowMap = mPtr->RowMap();
	const Epetra_BlockMap& colMap = mPtr->ColMap();
//	if(0==myRank){
////		cout << "column Map = " << colMap << endl;
////		cout << "row Map = " << rowMap << endl;
////		mPtr->Print(cout);
//
//		Field<double> xOverlapField(DISPL3D,2);
//		xOverlapField.setValue(2.0);
//		Epetra_Vector y(rowMap);
//		Epetra_Vector x(View,colMap,xOverlapField.getArray().get());
//		BOOST_CHECK(0==mPtr->Multiply1(false,x,y));
//		cout << y << endl;
//	} else {
////		cout << "column Map = " << colMap << endl;
////		cout << "row Map = " << rowMap << endl;
////		mPtr->Print(cout);
//		Field<double> xOverlapField(DISPL3D,2);
//		xOverlapField.setValue(2.0);
//		Epetra_Vector y(rowMap);
//		Epetra_Vector x(View,colMap,xOverlapField.getArray().get());
//		BOOST_CHECK(0==mPtr->Multiply1(false,x,y));
//		cout << y << endl;
//	}



	double *f = fN.getArray().get();
	Epetra_Vector rhs(View,rangeMap,f);
	Epetra_Vector lhs(View,domainMap,uOwnedField.getArray().get());

	linProblem.SetRHS(&rhs);
	linProblem.SetLHS(&lhs);
	BOOST_CHECK(0==linProblem.CheckInput());

	AztecOO solver(linProblem);
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
//	solver.SetAztecOption(AZ_conv,AZ_noscaled);
//	solver.SetAztecOption(AZ_precond, AZ_none);
	BOOST_CHECK(0==solver.CheckInput());
	solver.Iterate(500,1e-6);
	/*
	 * Write problem set up parameters to file
	 */
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,decomp.numPoints);
	PdVTK::writeField(grid,uOwnedField);
	PdVTK::writeField(grid,fN);
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utPimp_twoPointJacobian_np2.pvtu", comm.NumProc(), comm.MyPID());
	PdVTK::write(writer,grid);


}

shared_ptr<Epetra_VbrMatrix>
getOperator
(
		const vector<shared_ptr<StageComponentDirichletBc> >& bcArray,
		shared_ptr<Epetra_CrsGraph>& graphPtr,
		shared_ptr<RowStiffnessOperator>& jacobian
)
{
	std::cout << "Begin jacobian calculation\n";
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
	 * Create a searchable set for bc
	 */

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
					 * Put -1 because later the whole matrix is negated
					 */
					if(row==col) {
						for(std::size_t r=0;r<components.size();r++){
							*(diagPtrs[r]) = -1.0;
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
//		PRINT_3x3MATRIX(kPtr,std::cout);
//		PRINT_3x3MATRIX(kPtr+9,std::cout);
//		cout << endl;
		for(int* colPtr=cols; colPtr!=cols+numCol;colPtr++){
			/*
			 * Fill 'k'
			 */
			for(int c=0;c<3;c++){
				double *colK = k[c];
				for(int r=0;r<3;r++,kPtr++)
					colK[r] = - *kPtr;
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
	return shared_ptr<Epetra_VbrMatrix>(m);
}


shared_ptr<PdImp::PdImpOperator> getPimpOperator(PdGridData& decomp,Epetra_MpiComm& comm) {
	PdImp::PdImpOperator *op = new PdImp::PdImpOperator(comm,decomp);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op->addConstitutiveModel(fIntOperator);
	return shared_ptr<PdImp::PdImpOperator>(op);
}

void setRotatedValue(double* x, double magnitude){
	*(x+0) = magnitude*cos(theta);
	*(x+1) = magnitude*sin(theta);
	// no change to z
}

shared_ptr<Epetra_CrsGraph> getGraph(shared_ptr<RowStiffnessOperator>& jacobian){
	const Epetra_BlockMap& rowMap   = jacobian->getRowMap();
	const Epetra_BlockMap& colMap = jacobian->getColMap();

	/*
	 * Epetra Graph
	 */
	Pd_shared_ptr_Array<int> numCols = jacobian->getNumColumnsPerRow();
	Epetra_CrsGraph *graph = new Epetra_CrsGraph(Copy,rowMap,colMap,numCols.get(),true);
	for(int row=0;row<jacobian->getNumRows();row++){
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.getSize();
		/*
		 * Each row should have 4 entries
		 */
		BOOST_CHECK(2==numCol);
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

	test_suite* proc = BOOST_TEST_SUITE( "utPimp_twoPointJacobian_np2" );
	proc->add(BOOST_TEST_CASE( &computeJacobian ));
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
	 * This test only make sense for numProcs == 2
	 */
	if(2 != numProcs){
		std::cerr << "Unit test runtime ERROR: utPimp_twoPointJacobian_np2 is intended for 2 processors only." << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 2 ./utPimp_twoPointJacobian_np2" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
