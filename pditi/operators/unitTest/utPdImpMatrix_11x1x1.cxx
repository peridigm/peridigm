/*
 * utPimpMatrix_11x1x1.cxx
 *
 *  Created on: Mar 17, 2010
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
#include "../../pdneigh/NeighborhoodList.h"
#include "../PdImpMaterials.h"
#include "../PdITI_Operator.h"
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


static int myRank;
static int numProcs;

const int nx = 11;
const int ny = 1;
const int nz = 1;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = 0.0;
const double zLength = 1.0;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=.1;

void axialBarLinearSpacing() {

	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	decomp = getLoadBalancedDiscretization(decomp);
	int numPoints = decomp.numPoints;

	/*
	 * Create parallel communication maps
	 */
	const int vectorNDF=3;
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	Epetra_BlockMap rowMap   = PdQuickGrid::getOwnedMap  (comm, decomp, vectorNDF);
	Epetra_BlockMap colMap = PdQuickGrid::getOverlapMap(comm, decomp, vectorNDF);
	BOOST_CHECK(rowMap.NumMyElements()==numPoints);

	const PdImp::BulkModulus _K(130000.0);
	const PdImp::PoissonsRatio _MU(0.25);
	const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);
	PdITI::PdITI_Operator op(comm,list,decomp.cellVolume);

//	PdImp::PdImpOperator op(comm,decomp);
	op.addConstitutiveModel(fIntOperator);
	Field<double> uOwnedField = Field<double>(DISPL3D,list.get_num_owned_points());
	PdITI::SET(uOwnedField.getArray().get(),uOwnedField.getArray().get()+uOwnedField.getArray().getSize(),0.0);
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getJacobian(uOwnedField);

//	int n0[]  = {0,1,2};
//	int n1[]  = {0,1,2,3};
//	int n2[]  = {0,1,2,3,4};
//	int n3[]  = {1,2,3,4,5};
//	int n4[]  = {2,3,4,5,6};
//	int n5[]  = {3,4,5,6,7};
//	int n6[]  = {4,5,6,7,8};
//	int n7[]  = {5,6,7,8,9};
//	int n8[]  = {6,7,8,9,10};
//	int n9[]  = {7,8,9,10};
//	int n10[] = {8,9,10};
//	int* nPtr[]  = {n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10};
	const std::size_t s[] = {3,4,5,5,5,5,5,5,5,4,3};

	/*
	 * Epetra Graph
	 */
	Pd_shared_ptr_Array<int> numCols = jacobian->getNumColumnsPerRow();
	Epetra_CrsGraph graph(Copy,rowMap,colMap,numCols.get());
	for(int row=0;row<rowMap.NumMyElements();row++){
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.getSize();
		BOOST_CHECK(s[row]==numCol);
		int *cols = rowLIDs.get();
		BOOST_CHECK(0==graph.InsertMyIndices(row,numCol,cols));
	}
	BOOST_CHECK(0==graph.FillComplete());

	/*
	 * Epetra Matrix
	 * PERHAPS the 'operator' can keep its own copy of Graph, then this
	 * can be a 'View'
	 */
	Epetra_FEVbrMatrix m(Copy,graph);
	Epetra_SerialDenseMatrix k;
	k.Shape(vectorNDF,vectorNDF);
	for(int row=0;row<rowMap.NumMyElements();row++){
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numCol = rowLIDs.getSize();
		BOOST_CHECK(s[row]==numCol);
		int *cols = rowLIDs.get();
		BOOST_CHECK(0==m.BeginReplaceMyValues(row,numCol,cols));
		/*
		 * loop over columns in row and submit block entry
		 */
		Pd_shared_ptr_Array<double> actualK = jacobian->computeRowStiffness(row, rowLIDs);
		double *kPtr = actualK.get();
		for(std::size_t col=0;col<numCol;col++){
			/*
			 * Fill 'k'
			 */
			for(int c=0;c<3;c++){
				double *colK = k[c];
				for(int r=0;r<3;r++,kPtr++)
					colK[r] = *kPtr;
			}
			BOOST_CHECK(0==m.SubmitBlockEntry(k));
		}
		/*
		 * Finalize this row
		 */
		BOOST_CHECK(0==m.EndSubmitEntries());
	}

	BOOST_CHECK(0==m.FillComplete());

}

void probe() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	decomp = getLoadBalancedDiscretization(decomp);
	int numPoints = decomp.numPoints;
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);

	/*
	 * Create parallel communication maps
	 */
	const int vectorNDF=3;
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);

	const PdImp::BulkModulus _K(130000.0);
	const PdImp::PoissonsRatio _MU(0.25);
	const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon);
	PdITI::PdITI_Operator op(comm,list,decomp.cellVolume);
	op.addConstitutiveModel(fIntOperator);

	Epetra_BlockMap rowMap   = list.getOwnedMap(comm,vectorNDF);
	Epetra_BlockMap colMap = list.getOverlapMap(comm,vectorNDF);
	BOOST_CHECK(rowMap.NumMyElements()==numPoints);

	Field_NS::Field<double> u(Field_NS::DISPL3D,rowMap.NumMyElements());
	PdITI::SET(u.getArray().get(),u.getArray().get()+u.getArray().getSize(),0.0);

	/*
	 * Get first row of stiffness
	 */
	double delta = .00001;
	*(u.getArray().get())=delta;
	Field_NS::Field<double> force = Field_NS::getFORCE3D(rowMap.NumMyElements());
	double *f = force.getArray().get();
	PdITI::SET(f,f+force.getArray().getSize(),0.0);
	op.computeInternalForce(u,force);
	PdITI::SCALE_BY_VALUE(f,f+force.getArray().getSize(),1/delta);
//	for(;f!=fOwnedPtr.get()+3*4;f+=3)
//		std::cout << std::scientific << *f << " " << *(f+1) << " " << *(f+2) << " ";

}

void rowMap() {
	int NUM_RAND=1000;
	srand ( time(NULL) );

	std::vector< std::pair<int,int> > pairs(NUM_RAND);
	for(std::size_t i=0;i<pairs.size();i++)
		pairs[i] = std::make_pair(rand(),i);

	std::map<int,int> map(pairs.begin(),pairs.end());

	for(int i=0;i<NUM_RAND;i++){
		int randNum = pairs[i].first;
		BOOST_CHECK(i == map[randNum]);
	}

}

void TENSOR_PRODUCT(){
	double ux[] = {1,0,0};
	double m[9];
	PdITI::TENSOR_PRODUCT(ux,ux,m);
	const double tolerance = 1.0e-15;
	/*
	 * Column 1
	 */
	BOOST_CHECK_CLOSE(m[0],1.0,tolerance);
	BOOST_CHECK_CLOSE(m[1],0.0,tolerance);
	BOOST_CHECK_CLOSE(m[2],0.0,tolerance);

	/*
	 * Column 2
	 */
	BOOST_CHECK_CLOSE(m[3],0.0,tolerance);
	BOOST_CHECK_CLOSE(m[4],0.0,tolerance);
	BOOST_CHECK_CLOSE(m[5],0.0,tolerance);

	/*
	 * Column 3
	 */
	BOOST_CHECK_CLOSE(m[6],0.0,tolerance);
	BOOST_CHECK_CLOSE(m[7],0.0,tolerance);
	BOOST_CHECK_CLOSE(m[8],0.0,tolerance);

}

void SET(){
	double m[9];
	const double tolerance = 1.0e-15;
	double one(1.0f);
	PdITI::SET(m,m+9,one);
	for(int i=0;i<9;i++)
		BOOST_CHECK_CLOSE(m[i],one,tolerance);
}

void BOND(){
	double ux1[] = {0,0,0};
	double ux2[] = {0,0,0};
	double bond[3];
	double answers[3];
	for(int i=0;i<3;i++){
		ux1[i]=rand();
		ux2[i]=rand();
		answers[i] = ux2[i]-ux1[i];
	}

	PdITI::BOND(ux1,ux2,bond);
	const double tolerance = 1.0e-15;
	for(int i=0;i<3;i++)
		BOOST_CHECK_CLOSE(answers[i],bond[i],tolerance);

}

void UPDATE_GEOMETRY(){
	double x[] = {0,0,0};
	double u[] = {0,0,0};
	double y[3];
	double answers[3];
	for(int i=0;i<3;i++){
		x[i]=rand();
		u[i]=rand();
		answers[i] = x[i]+u[i];
	}

	PdITI::UPDATE_GEOMETRY(x,u,y);
	const double tolerance = 1.0e-15;
	for(int i=0;i<3;i++)
		BOOST_CHECK_CLOSE(answers[i],y[i],tolerance);

}

void MAGNITUDE(){
	double x[] = {0,0,0};
	double sum(0.0);
	for(int i=0;i<3;i++){
		x[i]=rand()%100;
		sum += x[i]*x[i];
	}
	double answer = sqrt(sum);
	const double tolerance = 1.0e-15;
	double norm = PdITI::MAGNITUDE(x);
	BOOST_CHECK_CLOSE(answer,norm,tolerance);
}

void NORMALIZE(){
	double x[] = {0,0,0};
	double answer[3];
	double sum(0.0);
	for(int i=0;i<3;i++){
		x[i]=rand()%100;
		sum += x[i]*x[i];
	}
	double norm = sqrt(sum);
	for(int i=0;i<3;i++){
		answer[i]=x[i]/norm;
	}
	const double tolerance = 1.0e-15;
	PdITI::NORMALIZE(x);
	for(int i=0;i<3;i++){
		BOOST_CHECK_CLOSE(answer[i],x[i],tolerance);
	}
}

void SUMINTO(){
	double x[] = {0,0,0};
	double intoMe[3] = {0,0,0};
	double answer[3];
	for(int i=0;i<3;i++){
		x[i]=rand()%100;
		intoMe[i]=rand()%100;
		answer[i]=x[i]+intoMe[i];
	}
	const double tolerance = 1.0e-15;
	PdITI::SUMINTO(x,x+3,intoMe);
	for(int i=0;i<3;i++){
		BOOST_CHECK_CLOSE(answer[i],intoMe[i],tolerance);
	}
}

void SUBTRACTINTO(){
	double x[] = {0,0,0};
	double intoMe[3] = {0,0,0};
	double answer[3];
	for(int i=0;i<3;i++){
		x[i]=rand()%100;
		intoMe[i]=rand()%100;
		answer[i]=intoMe[i]-x[i];
	}
	const double tolerance = 1.0e-15;
	PdITI::SUBTRACTINTO(x,x+3,intoMe);
	for(int i=0;i<3;i++){
		BOOST_CHECK_CLOSE(answer[i],intoMe[i],tolerance);
	}
}

void COPY(){
	double x[] = {0,0,0};
	double intoMe[3] = {0,0,0};
	double answer[3];
	for(int i=0;i<3;i++){
		x[i]=rand()%100;
		answer[i]=x[i];
	}
	const double tolerance = 1.0e-15;
	PdITI::COPY(x,x+3,intoMe);
	for(int i=0;i<3;i++){
		BOOST_CHECK_CLOSE(answer[i],intoMe[i],tolerance);
	}
}

void SCALE_BY_VALUE(){
	double x[] = {0,0,0};
	double answer[3];
	double value= rand();
	for(int i=0;i<3;i++){
		x[i]=rand()%100;
		answer[i]=value*x[i];
	}
	const double tolerance = 1.0e-15;
	PdITI::SCALE_BY_VALUE(x,x+3,value);
	for(int i=0;i<3;i++){
		BOOST_CHECK_CLOSE(answer[i],x[i],tolerance);
	}
}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPimpMatrix_11x1x1" );
	proc->add(BOOST_TEST_CASE( &axialBarLinearSpacing ));


	proc->add(BOOST_TEST_CASE( &BOND ));
	proc->add(BOOST_TEST_CASE( &UPDATE_GEOMETRY ));
	proc->add(BOOST_TEST_CASE( &MAGNITUDE ));
	proc->add(BOOST_TEST_CASE( &NORMALIZE ));
	proc->add(BOOST_TEST_CASE( &TENSOR_PRODUCT ));
	proc->add(BOOST_TEST_CASE( &SET ));
	proc->add(BOOST_TEST_CASE( &SUMINTO ));
	proc->add(BOOST_TEST_CASE( &SUBTRACTINTO ));
	proc->add(BOOST_TEST_CASE( &COPY ));
	proc->add(BOOST_TEST_CASE( &SCALE_BY_VALUE ));
	proc->add(BOOST_TEST_CASE( &rowMap ));
	proc->add(BOOST_TEST_CASE( &probe ));
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
		std::cerr << "Unit test runtime ERROR: utPimpMatrix_11x1x1 is intended for \"serial\" run only and makes sense on 1 processor" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 1 ./utPimpMatrix_11x1x1" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
