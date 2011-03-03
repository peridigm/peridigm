/*
 * utPimp_twoPointJacobian.cxx
 *
 *  Created on: Apr 6, 2010
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
#include "../PdITI_Utilities.h"
#include "../IsotropicElasticConstitutiveModel.h"
#include <set>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialSymDenseMatrix.h>
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


static int myRank;
static int numProcs;

const int nx = 2;
const int ny = 1;
const int nz = 1;
const double lX = 1.0;
const double theta = M_PI/6.0;
const double xStart  = -lX/4.0;
const double xLength =  lX;
const double yStart  = -0.5;
const double yLength =  1.0;
const double zStart  = -0.5;
const double zLength =  1.0;
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


void assertPoint_0_at_Origin(double* x0);
void assertPoint_1_at_lX_by_2(double* x1);
void setRotatedValue(double* x, double magnitude);
void assertRotatedPoint_1(double *x1);
void assertForce(const Field<double>& fOwnedField);
Pd_shared_ptr_Array<double> computeAnalytical3x3Stiffness(const Field<double>& uOwnedField);

/*
 * Unit Test Description
 * 1) Begins with single bond along x-axis
 * 2) Creates displacement that is composed of:
 * 2a) elongation of bond
 * 2b) rigid body rotation of bond
 * 3) Computes force in deformed configuration
 * 4) Computes jacobian at deformed configuration
 */
void twoPointJacobian() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	decomp = getLoadBalancedDiscretization(decomp);
	int numPoints = decomp.numPoints;
	BOOST_CHECK(2==numPoints);
	Field<double> uOwnedField = Field<double>(DISPL3D,numPoints);
	PdITI::SET(uOwnedField.getArray().get(),uOwnedField.getArray().get()+uOwnedField.getArray().getSize(),0.0);

	{
		/*
		 * Point 0 and Point 1
		 */
		double *x0 = decomp.myX.get();
		double *x1 = x0+3;
		assertPoint_0_at_Origin(x0);
		assertPoint_1_at_lX_by_2(x1);

		/*
		 * Create displacement field at point 1
		 * 1) create final configuration
		 * 2) difference final configuration with starting configuration to get displacement
		 */
		double y[]={0.0,0.0,0.0};
		double yMag = x1[0]+delta;
		setRotatedValue(y,yMag);

		/*
		 * Set displacement on point 1
		 */
		double *u1 = uOwnedField.getArray().get()+3;
		u1[0] = y[0] - x1[0];
		u1[1] = y[1] - x1[1];
		u1[2] = y[2] - x1[2];
	}


	/*
	 * Create parallel communication maps
	 */
	const int vectorNDF=3;
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	Epetra_BlockMap rowMap   = PdQuickGrid::getOwnedMap  (comm, decomp, vectorNDF);
	Epetra_BlockMap colMap = PdQuickGrid::getOverlapMap(comm, decomp, vectorNDF);
	BOOST_CHECK(rowMap.NumMyElements()==numPoints);
	PdImp::PdImpOperator op(comm,decomp);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);


	/*
	 * Compute force
	 */
	Field_NS::TemporalField<double> force = Field_NS::TemporalField<double>(Field_NS::FORCE3D,numPoints);
	op.computeInternalForce(uOwnedField,force.getField(Field_NS::FieldSpec::STEP_NP1));
	assertForce(force.getField(Field_NS::FieldSpec::STEP_NP1));
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getRowStiffnessOperator(uOwnedField,horizon);

	/*
	 * Assert Diagonal Entries
	 */
	{
		int row(0),col(0);
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		BOOST_CHECK(2 == rowLIDs.getSize());
		BOOST_CHECK(0 == *(rowLIDs.get_shared_ptr().get()));
		BOOST_CHECK(1 == *(rowLIDs.get_shared_ptr().get()+1));
		const double *k = jacobian->computeRowStiffness(row, rowLIDs).get()+9*col;
		/*
		 * This is the hand calculated 'diagonal' entry
		 */
		const Pd_shared_ptr_Array<double> kA = computeAnalytical3x3Stiffness(uOwnedField);
		const double *kAns = kA.get();
//		Pimp::PRINT_3x3MATRIX(k,std::cout);
//		Pimp::PRINT_3x3MATRIX(kAns,std::cout);
		 double tolerance = 1.0e-15;
		for(int i=0;i<8;i++,k++,kAns++){
			BOOST_CHECK_SMALL((*k)+(*kAns),tolerance);
		}
		tolerance = 1.0e-10;
		for(int i=8;i<9;i++,k++,kAns++){
			BOOST_CHECK_CLOSE((*k),-(*kAns),tolerance);
		}


//		std::cout << "ROW = " << row <<  std::endl;
//		int *cols = rowLIDs.get();
//		for(std::size_t col=0;col<rowLIDs.getSize();col++){
//			std::cout << "\tCOL = " << cols[col] << std::endl;
//			double *c = stiffness.get()+9*col;
//			Pimp::PRINT_3x3MATRIX(c,std::cout);
//		}
//		std::cout << endl;
	}
	{
		int row(1);
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		BOOST_CHECK(2 == rowLIDs.getSize());
//		cout << "Row = " << row << "; first col id = " << *(rowLIDs.get_shared_ptr().get()) << endl;
//		cout << "Row = " << row << "; 2nd col id = " << *(rowLIDs.get_shared_ptr().get()+1) << endl;
		BOOST_CHECK(0 == *(rowLIDs.get_shared_ptr().get()));
		BOOST_CHECK(1 == *(rowLIDs.get_shared_ptr().get()+1));
		Pd_shared_ptr_Array<double> stiffness = jacobian->computeRowStiffness(row, rowLIDs);

		/*
		 * This analytical jacobian computes the diagonal entry;
		 * Therefore, we have to offset the pointer returned by 9 for comparison purposes.
		 */
		const double *k = jacobian->computeRowStiffness(row, rowLIDs).get()+9;
		const Pd_shared_ptr_Array<double> kA = computeAnalytical3x3Stiffness(uOwnedField);
		const double *kAns = kA.get();

		double tolerance = 1.0e-15;
		for(int i=0;i<8;i++,k++,kAns++){
			BOOST_CHECK_SMALL((*k)+(*kAns),tolerance);
		}
		tolerance = 1.0e-10;
		for(int i=8;i<9;i++,k++,kAns++){
			BOOST_CHECK_CLOSE((*k),-(*kAns),tolerance);
		}

//		std::cout << "Diagonal of Jacobian by hand:" << endl;
//		PdImp::PRINT_3x3MATRIX(kA.get(),std::cout);
		std::cout << "Computed row:" << std::endl;
		std::cout << "ROW = " << row <<  std::endl;
		int *cols = rowLIDs.get();
		for(std::size_t col=0;col<rowLIDs.getSize();col++){
			std::cout << "\tCOL = " << cols[col] << std::endl;
			double *c = stiffness.get()+9*col;
			PdITI::PRINT_3x3MATRIX(c,std::cout);
		}
		std::cout << endl;
	}
	{
//		Pd_shared_ptr_Array<double> kAnswer = computeAnalytical3x3Stiffness(uOwnedField);
//		double *c = kAnswer.get();
//		Pimp::PRINT_3x3MATRIX(c,std::cout);
	}

	{
		/*
		 * Probe for tangent
		 */
		/*
		 * Create displacement field at point 1
		 * 1) create final configuration
		 * 2) difference final configuration with starting configuration to get displacement
		 */
		double y[]={0.0,0.0,0.0};
		double *x1 = decomp.myX.get()+3;
		double yMag = *(x1)+delta;
		setRotatedValue(y,yMag);

		/*
		 * Set displacement on point 1
		 */
		double *u1 = uOwnedField.getArray().get()+3;
		u1[0] = y[0] - x1[0];
		u1[1] = y[1] - x1[1];
		u1[2] = y[2] - x1[2];
		op.computeInternalForce(uOwnedField,force.getField(Field_NS::FieldSpec::STEP_N));
		double *fN = force.getField(Field_NS::FieldSpec::STEP_N).getArray().get();
		/*
		 * Now probe
		 * NOTE: The following two rows of forces should match with the
		 * two matrices printed for ROW = 1;  Note that the matrices
		 * need to be placed side by side but it is easy to see that they
		 * are the same when things are correct.
		 *
		 */
		// X-Dir
//		std::cout << "PROBE " << std::endl;
		u1[0] = y[0] - x1[0] + 1.0e-6*delta;
		op.computeInternalForce(uOwnedField,force.getField(Field_NS::FieldSpec::STEP_NP1));
		double *fNP1X = force.getField(Field_NS::FieldSpec::STEP_NP1).getArray().get();
		std::cout <<  std::scientific << "\t";
		for(int r=0;r<6;r++)
			std::cout <<   (*(fNP1X+r) - *(fN+r))/(1.0e-6*delta) << " ";

		// Y-Dir
		std::cout <<  std::scientific << "\n";
		u1[0] = y[0] - x1[0];
		u1[1] = y[1] - x1[1] + 1.0e-6*delta;
		op.computeInternalForce(uOwnedField,force.getField(Field_NS::FieldSpec::STEP_NP1));
		fNP1X = force.getField(Field_NS::FieldSpec::STEP_NP1).getArray().get();
		std::cout <<  std::scientific << "\t";
		for(int r=0;r<6;r++)
			std::cout <<   (*(fNP1X+r) - *(fN+r))/(1.0e-6*delta) << " ";

		std::cout << std::endl;


		{
			u1[0] = y[0] - x1[0];
			u1[1] = y[1] - x1[1];
			jacobian = op.getRowStiffnessOperator(uOwnedField,horizon);
			int row=0;
			Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
			Pd_shared_ptr_Array<double> stiffness = jacobian->computeRowStiffness(row, rowLIDs);
//			std::cout << "ROW = " << row <<  std::endl;
//			int *cols = rowLIDs.get();
//			for(std::size_t col=0;col<rowLIDs.getSize();col++){
//				std::cout << "\tCOL = " << cols[col] << std::endl;
//				double *c = stiffness.get()+9*col;
//				PdImp::PRINT_3x3MATRIX(c,std::cout);
//			}
//			std::cout << endl;
		}

	}
}

Pd_shared_ptr_Array<double> computeAnalytical3x3Stiffness(const Field<double>& uOwnedField){

	double x1 = bond;
	double y0Mag = bond+delta;
	double y0[]={0.0,0.0,0.0};
	y0[2]=0.0;
	setRotatedValue(y0,y0Mag);
	const double tolerance = 1.0e-15;
	BOOST_CHECK_CLOSE(PdITI::NORMALIZE(y0),y0Mag,tolerance);
	Pd_shared_ptr_Array<double> MY0xMY0(9);
	double *k = MY0xMY0.get();

	PdITI::TENSOR_PRODUCT(y0,y0,k);
	PdITI::SCALE_BY_VALUE(k,k+9,-1.0);
	k[0] += 1.0;
	k[4] += 1.0;
	k[8] += 1.0;
	PdITI::SCALE_BY_VALUE(k,k+9,-x1/y0Mag);
	k[0] += 1.0;
	k[4] += 1.0;
	k[8] += 1.0;

	double c = 18.0*isotropicSpec.bulkModulus()/x1/x1;
	PdITI::SCALE_BY_VALUE(k,k+9,c);

	return MY0xMY0;
}

void assertForce(const Field<double>& fOwnedField){
	double x1 = bond;
	double y = bond+delta;
	double fMag = 18.0*isotropicSpec.bulkModulus()*(y-x1)/x1/x1;
	double fAnswer[3] = {0,0,0};
	setRotatedValue(fAnswer,fMag);
	const double *f = fOwnedField.getArray().get();
	const double tolerance = 1.0e-15;
	BOOST_CHECK_CLOSE(*(f+0),fAnswer[0],tolerance*100);
	BOOST_CHECK_SMALL(*(f+1)-fAnswer[1],tolerance);
	BOOST_CHECK_SMALL(*(f+2)-fAnswer[2],tolerance);
	BOOST_CHECK_CLOSE(*(f+3),-fAnswer[0],tolerance*100);
	BOOST_CHECK_SMALL(*(f+4)+fAnswer[1],tolerance);
	BOOST_CHECK_SMALL(*(f+5)+fAnswer[2],tolerance);
}


void setRotatedValue(double* x, double magnitude){
	*(x+0) = magnitude*cos(theta);
	*(x+1) = magnitude*sin(theta);
	// no change to z
}

void assertRotatedPoint_1(double *x1){
	const double tolerance = 1.0e-15;
	BOOST_CHECK_CLOSE(*(x1+0),(lX/2.0)*cos(theta),tolerance);
	BOOST_CHECK_CLOSE(*(x1+1),(lX/2.0)*sin(theta),tolerance);
	BOOST_CHECK_CLOSE(*(x1+2),0.0,tolerance);
}

void assertPoint_0_at_Origin(double* x0){
	const double tolerance = 1.0e-15;
	BOOST_CHECK_CLOSE(*(x0+0),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x0+1),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x0+2),0.0,tolerance);
}

void assertPoint_1_at_lX_by_2(double* x1){
	const double tolerance = 1.0e-15;
	BOOST_CHECK_CLOSE(*(x1+0),lX/2.0,tolerance);
	BOOST_CHECK_CLOSE(*(x1+1),0.0,tolerance);
	BOOST_CHECK_CLOSE(*(x1+2),0.0,tolerance);
}

void rotatePoint(double* x1){

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPimp_twoPointJacobian" );
	proc->add(BOOST_TEST_CASE( &twoPointJacobian ));
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
		std::cerr << "Unit test runtime ERROR: utPimp_twoPointJacobian is intended for \"serial\" run only and makes sense on 1 processor" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 1 ./utPimp_twoPointJacobian" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
