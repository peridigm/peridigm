/*
 * utPimp_probeJacobian.cxx
 *
 *  Created on: Apr 13, 2010
 *      Author: jamitch
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include "../PdImpMpiFixture.h"
#include "PdNeighborhood.h"
#include "../../pdneigh/NeighborhoodList.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdNeighborhood.h"
#include "PdZoltan.h"
#include "Field.h"
#include "../PdImpMaterials.h"
#include "../PdImpOperator.h"
#include "../PdITI_Operator.h"
#include "../PdImpOperatorUtilities.h"
#include "../IsotropicElasticConstitutiveModel.h"
#include "PdVTK.h"
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

const int nx = 4;
const int ny = 4;
const int nz = 8;
const double lX = 1.0;
const double lY = 1.0;
const double lZ = 1.0;
const double xStart  = -lX/2.0/nx;
const double xLength =  lX;
const double yStart  = -lY/2.0/ny;
const double yLength =  lY;
const double zStart  = -lZ/2.0/nz;
const double zLength =  lZ;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;
const double horizon=1.01*sqrt(pow(lX/nx,2)+pow(lY/ny,2)+pow(lZ/nz,2));
const PdImp::BulkModulus _K(130000.0);
const PdImp::PoissonsRatio _MU(0.0);
const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);

using PdVTK::writeField;

void probe() {
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec,PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	decomp = getLoadBalancedDiscretization(decomp);
	int numPoints = decomp.numPoints;
	PDNEIGH::NeighborhoodList list(decomp.zoltanPtr.get(),numPoints,decomp.myGlobalIDs,decomp.myX,horizon);

	BOOST_CHECK(numCells==numPoints);
	Epetra_MpiComm comm = Epetra_MpiComm(MPI_COMM_WORLD);
	Field<double> uOwnedField = PdITI::getPureShearXY(Field_NS::Field<double>(Field_NS::COORD3D,decomp.myX,numPoints));
	Field<double> probeField = Field_NS::Field<double>(Field_NS::COORD3D,numPoints);

	/*
	 * Compute Tangent
	 * For each row, probe force operator and compare
	 */
	Field_NS::TemporalField<double> force = Field_NS::TemporalField<double>(Field_NS::FORCE3D,numPoints);
	PdImp::PdImpOperator op(comm,decomp);
	PdITI::PdITI_Operator pditiOp(comm,list,decomp.cellVolume);
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op.addConstitutiveModel(fIntOperator);
	pditiOp.addConstitutiveModel(fIntOperator);


	/*
	 * Compute jacobian in deformed state
	 */
	std::tr1::shared_ptr<RowStiffnessOperator> jacobian = op.getRowStiffnessOperator(uOwnedField,horizon);

	for(std::size_t row=0;row<numPoints;row++){
		Pd_shared_ptr_Array<int> rowLIDs = jacobian->getColumnLIDs(row);
		std::size_t numColsRow = rowLIDs.getSize();
		std::vector< std::pair<int,int> > pairs(numColsRow);
		{
			int *cols = rowLIDs.get();
			for(std::size_t i=0;i<pairs.size();i++,cols++){
				pairs[i] = std::make_pair(*cols,i);
			}
		}
		std::map<int,int> map(pairs.begin(),pairs.end());
		const double *k = jacobian->computeRowStiffness(row, rowLIDs).get();

		/*
		 * Initialize probe displacement
		 */
		int loc = 3*row;
		double delta = horizon*1.0e-4;
		double kProbe[9];
		double *u = uOwnedField.getArray().get();
		{
			/*
			 * probe in x-dir
			 */
			PdITI::SET(probeField.getArray().get(),probeField.getArray().get()+probeField.getArray().getSize(),0.0);
			PdITI::SUMINTO(u,uOwnedField.getArray().end(),probeField.getArray().get());
			double *p = probeField.getArray().get()+loc;
			*p += delta;
			pditiOp.computeInternalForce(probeField,force.getField(Field_NS::FieldSpec::STEP_NP1));
			double *f1 = force.getField(Field_NS::FieldSpec::STEP_NP1).getArray().get();

			PdITI::SET(probeField.getArray().get(),probeField.getArray().get()+probeField.getArray().getSize(),0.0);
			PdITI::SUMINTO(u,uOwnedField.getArray().end(),probeField.getArray().get());
			p = probeField.getArray().get()+loc;
			*p -= delta;
			pditiOp.computeInternalForce(probeField,force.getField(Field_NS::FieldSpec::STEP_N));
			double *f0 = force.getField(Field_NS::FieldSpec::STEP_N).getArray().get();

			kProbe[0] = (f1[loc+0]-f0[loc+0])/2.0/delta;
			kProbe[3] = (f1[loc+1]-f0[loc+1])/2.0/delta;
			kProbe[6] = (f1[loc+2]-f0[loc+2])/2.0/delta;

		}

		{
			/*
			 * probe in y-dir
			 */
			PdITI::SET(probeField.getArray().get(),probeField.getArray().get()+probeField.getArray().getSize(),0.0);
			PdITI::SUMINTO(u,uOwnedField.getArray().end(),probeField.getArray().get());
			double *p = probeField.getArray().get()+loc+1;
			*p += delta;
			pditiOp.computeInternalForce(probeField,force.getField(Field_NS::FieldSpec::STEP_NP1));
			double *f1 = force.getField(Field_NS::FieldSpec::STEP_NP1).getArray().get();

			PdITI::SET(probeField.getArray().get(),probeField.getArray().get()+probeField.getArray().getSize(),0.0);
			PdITI::SUMINTO(u,uOwnedField.getArray().end(),probeField.getArray().get());
			p = probeField.getArray().get()+loc+1;
			*p -= delta;
			pditiOp.computeInternalForce(probeField,force.getField(Field_NS::FieldSpec::STEP_N));
			double *f0 = force.getField(Field_NS::FieldSpec::STEP_N).getArray().get();

			kProbe[1] = (f1[loc+0]-f0[loc+0])/2.0/delta;
			kProbe[4] = (f1[loc+1]-f0[loc+1])/2.0/delta;
			kProbe[7] = (f1[loc+2]-f0[loc+2])/2.0/delta;
		}

		{
			/*
			 * probe in z-dir
			 */
			PdITI::SET(probeField.getArray().get(),probeField.getArray().get()+probeField.getArray().getSize(),0.0);
			PdITI::SUMINTO(u,uOwnedField.getArray().end(),probeField.getArray().get());
			double *p = probeField.getArray().get()+loc+2;
			*p += delta;
			pditiOp.computeInternalForce(probeField,force.getField(Field_NS::FieldSpec::STEP_NP1));
			double *f1 = force.getField(Field_NS::FieldSpec::STEP_NP1).getArray().get();

			PdITI::SET(probeField.getArray().get(),probeField.getArray().get()+probeField.getArray().getSize(),0.0);
			PdITI::SUMINTO(u,uOwnedField.getArray().end(),probeField.getArray().get());
			p = probeField.getArray().get()+loc+2;
			*p -= delta;
			pditiOp.computeInternalForce(probeField,force.getField(Field_NS::FieldSpec::STEP_N));
			double *f0 = force.getField(Field_NS::FieldSpec::STEP_N).getArray().get();

			kProbe[2] = (f1[loc+0]-f0[loc+0])/2.0/delta;
			kProbe[5] = (f1[loc+1]-f0[loc+1])/2.0/delta;
			kProbe[8] = (f1[loc+2]-f0[loc+2])/2.0/delta;

		}

//		{
//			std::cout << "ROW # " << (row+1) << std::endl;
//			int *cols = rowLIDs.get();
//			for(;cols!=rowLIDs.end();cols++){
//				int col = map[*cols];
//				double *kIC = k+col*9;
//				std::cout <<"\tCOL = " << (*cols+1) << "\n";
//				Pimp::PRINT_3x3MATRIX(kIC,std::cout);
//			}
//		}

		int diagCol = map[row];
		const double *kII = k+diagCol*9;
//		std::cout << "Analytical KII\n";
//		Pimp::PRINT_3x3MATRIX(kII,std::cout);
//		std::cout << "Probe KII\n";
//		Pimp::PRINT_3x3MATRIX(kProbe,std::cout);

		double tolerance = 1.0e-3;
		for(int i=0;i<9;i++){
			BOOST_CHECK_CLOSE(kII[i],kProbe[i],tolerance);
		}
	}


}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;

	test_suite* proc = BOOST_TEST_SUITE( "utPimp_probeJacobian" );
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
		std::cerr << "Unit test runtime ERROR: utPimp_probeJacobianNu.0.0 is intended for \"serial\" run only and makes sense on 1 processor" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 1 ./utPimp_probeJacobianNu.0.0" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}
