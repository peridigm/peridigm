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
#include <tr1/memory>
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdBondFilter.h"
#include "PdVTK.h"
#include "Field.h"
#include "../utPdITI.h"
#include "../../PdImpMpiFixture.h"
#include "../../../pdneigh/NeighborhoodList.h"
#include "../../PdImpMaterials.h"
#include "../../PdITI_Operator.h"
#include "../../PdITI_Utilities.h"
#include "../../IsotropicElasticConstitutiveModel.h"
#include "../../DirichletBcSpec.h"
#include "../../StageFunction.h"
#include "../../StageComponentDirichletBc.h"
#include "../../ComponentDirichletBcSpec.h"
#include <Epetra_LinearProblem.h>
#include <AztecOO.h>
#include "Epetra_RowMatrix.h"
#include <Epetra_FEVbrMatrix.h>

#include <iostream>


using namespace PdQuickGrid;
using namespace PdBondFilter;
using namespace PdVTK;
using namespace Field_NS;
using namespace PdImp;
using namespace PdITI;
using std::tr1::shared_ptr;
const int vectorNDF=3;
//using namespace boost::unit_test;
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
		if(0!=graph->InsertMyIndices(row,numCol,cols)){
			std::string message("graph->InsertMyIndices(row,numCol,cols)\n");
			message += "\t0!=graph->FillComplete()";
			throw std::runtime_error(message);
		}
	}

	if(0!=graph->FillComplete()){
		std::string message("getGraph(shared_ptr<RowStiffnessOperator>& jacobian)\n");
		message += "\t0!=graph->FillComplete()";
		throw std::runtime_error(message);
	}
	return shared_ptr<Epetra_CrsGraph>(graph);
}

shared_ptr<Epetra_RowMatrix>
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

		if(0!=m->BeginReplaceMyValues(row,numCol,cols)){
			std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
			message += "\t0!=m->BeginReplaceMyValues(row,numCol,cols)";
			throw std::runtime_error(message);
		}

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
					double *stiffPtr=0;
					for(std::size_t r=0;r<components.size();r++){
						stiffPtr = kPtr+3*components[r];
						for(int r=0;r<3;r++)
							stiffPtr[r]=0;
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
					colK[r] = - *kPtr;
			}

			if(0!=m->SubmitBlockEntry(k)){
				std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
				message += "\t0!=m->SubmitBlockEntry(k)";
				throw std::runtime_error(message);
			}

		}


		/*
		 * Finalize this row
		 */
		if(0!=m->EndSubmitEntries()){
			std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
			message += "\t0!=m->EndSubmitEntries()";
			throw std::runtime_error(message);
		}

	}

	if(0!=m->FillComplete()){
		std::string message("utPdITI::getOperator(bcArray,graphPtr,jacobian)\n");
		message += "\t0!=m->FillComplete()";
		throw std::runtime_error(message);
	}


	std::cout << "END jacobian calculation\n";
	return shared_ptr<Epetra_RowMatrix>(m);
}



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

	FinitePlane crackPlane=getYZ_CrackPlane();
//	shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new FinitePlaneFilter(crackPlane));
	shared_ptr<BondFilter> filterPtr = shared_ptr<BondFilter>(new BondFilterDefault());
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

	/*
	 * Get points for bc's
	 */
	/*
	 * Note that we are looking for the first three planes of points on ends
	 */
	double searchRadius=1.1*dx;
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
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getJacobian(uOwnedField);

	/*
	 * Create graph
	 */
	shared_ptr<Epetra_CrsGraph> graphPtr = getGraph(jacobian);

	/*
	 * Create Epetra_RowMatrix
	 */
	shared_ptr<Epetra_RowMatrix> mPtr = getOperator(bcs,graphPtr,jacobian);

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
//	BOOST_CHECK(0==linProblem.CheckInput());

	AztecOO solver(linProblem);
	solver.SetAztecOption(AZ_precond, AZ_Jacobi);
//	BOOST_CHECK(0==solver.CheckInput());
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
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("utCrackOpeningDemo_npX.pvtu", comm.NumProc(), comm.MyPID());
	PdVTK::write(writer,grid);

}

//bool init_unit_test_suite()
//{
//	// Add a suite for each processor in the test
//	bool success=true;
//	test_suite* proc = BOOST_TEST_SUITE( "utCrackOpeningDemo_npX" );
//	proc->add(BOOST_TEST_CASE( &crackOpeningDemo ));
//	framework::master_test_suite().add( proc );
//	return success;
//}
//
//
//bool init_unit_test()
//{
//	init_unit_test_suite();
//	return true;
//}

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

	crackOpeningDemo();

	// Initialize UTF
//	return unit_test_main( init_unit_test, argc, argv );
}
