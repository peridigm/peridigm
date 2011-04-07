/*! \file utVTK_ioExample_np4.cxx */

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ************************************************************************
//@HEADER

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <tr1/memory>
#include "PdZoltan.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdutMpiFixture.h"
#include "PdVTK.h"
#include <iostream>
#include "PdMaterialUtilities.h"


#include "Field.h"
using namespace Field_NS;

using namespace PdQuickGrid;
using namespace Pdut;
using std::tr1::shared_ptr;
using namespace boost::unit_test;
using std::cout;

static int myRank;
static int numProcs;
const int nx = 10;
const int ny = 10;
const int nz = 1;
const double xStart = 0.0;
const double xLength = 1.0;
const double yStart = 0.0;
const double yLength = 1.0;
const double zStart = -0.5;
const double zLength = 1.0;
const PdQPointSet1d xSpec(nx,xStart,xLength);
const PdQPointSet1d ySpec(ny,yStart,yLength);
const PdQPointSet1d zSpec(nz,zStart,zLength);
const int numCells = nx*ny*nz;

#ifndef PIMPMPIFIXTURE_H_
#define PIMPMPIFIXTURE_H_
#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

namespace PdImpRunTime {

struct PimpMpiFixture {
private:
	static bool finalized;
#ifdef HAVE_MPI
	Epetra_MpiComm Comm;
#else
	Epetra_SerialComm Comm;
#endif

#ifdef HAVE_MPI
	PimpMpiFixture(int argc, char* argv[] ) : Comm(MPI_COMM_WORLD) {
		Comm = Epetra_MpiComm(MPI_COMM_WORLD);
	}
#else
	PimpMpiFixture() : Comm()  {
		Comm = Epetra_MpiComm(MPI_COMM_WORLD);
	}
#endif
public:

	static PimpMpiFixture getPimpMPI(int argc, char* argv[]){
#ifdef HAVE_MPI
		MPI_Init(&argc, &argv);
		return PimpMpiFixture(argc,argv);
#else
		return PimpMpiFixture();
#endif
	}
	~PimpMpiFixture(){
#ifdef HAVE_MPI
		if(!finalized){
			MPI_Finalize();
			finalized=true;
		}
#endif
	}

	const Epetra_Comm& getEpetra_Comm() { return Comm; }
};

bool PimpMpiFixture::finalized=false;

} // namespace PimpRunTime

#endif // PIMPMPIFIXTURE_H_


Field<double> getPureShearXY(double gamma, const Field<double>& X, Field<double>& U){
	std::size_t numPoints = X.getNumPoints();
	double *u = U.getArray().get();
	const double *x = X.getArray().get();

	for(std::size_t i=0;i<numPoints;i++){
		int p=3*i;
		u[p]=gamma*x[p+1];
		u[p+1]=0;
		u[p+2]=0;
	}
	return U;
}

PdGridData getGrid() {

	double dx = xSpec.getCellSize();
	double dy = ySpec.getCellSize();
	double horizon = sqrt(dx*dx+dy*dy);
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);

	// This reload balances
	decomp = getLoadBalancedDiscretization(decomp);
	return decomp;
}



void utVTK_ioExample()
{
	PdGridData pdGridData = getGrid();
	int numPoints = pdGridData.numPoints;

	/*
	 * NOTE
	 *
	 * FieldSpec is IMMUTABLE
	 *
	 */

	/*
	 * Create Spec(s) From Scratch
	 */
	const FieldSpec myRankSpec(FieldSpec::DEFAULT_FIELDTYPE,FieldSpec::SCALAR,"MyRank");
	const FieldSpec displacementSpec(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR3D, "Displacement");
	const FieldSpec velocitySpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D, "v");
	const FieldSpec accelerationSpec(FieldSpec::ACCELERATION,FieldSpec::VECTOR3D, "a");

	/*
	 * Use existing spec (these are equivalent to the above -- but different 'names' on output
	 */
	const FieldSpec uSpec(DISPL3D);
	const FieldSpec vSpec(VELOC3D);
	const FieldSpec aSpec(ACCEL3D);

	const FieldSpec volSpec(VOLUME);
	const FieldSpec wSpec(WEIGHTED_VOLUME);
	const FieldSpec thetaSpec(DILATATION);

	/*
	 * This is not required in Peridigm
	 */
	Field<double> X(COORD3D,pdGridData.myX,numPoints);
	Field<double> uField(uSpec,numPoints), vField(vSpec,numPoints),aField(aSpec,numPoints);
	Field<double> wField(wSpec,numPoints), thetaField(thetaSpec,numPoints);
	Field<int> rankField(myRankSpec,numPoints);
	uField.setValue(0.0); vField.setValue(0.0); aField.setValue(0.0);
	wField.setValue(0.0); thetaField.setValue(0.0);
	rankField.setValue(myRank);

	/*
	 * RAW POINTERS; GET THESE from Epetra_Vector
	 */
	const double *xPtr = X.getArray().get();
	double *uPtr = uField.getArray().get();
	double *vPtr = vField.getArray().get();
	double *aPtr = aField.getArray().get();
	double *volPtr = pdGridData.cellVolume.get();
	double *wPtr = wField.getArray().get();
	double *thetaPtr = thetaField.getArray().get();
	const int *neighPtr = pdGridData.neighborhood.get();

	/*
	 * Create VTK unstructured grid
	 */
	vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(pdGridData.myX,numPoints);

	/*
	 * Create example 'collection' writers
	 * 1) ascii
	 * 2) binary
	 */
	PdVTK::CollectionWriter asciiWriter("utVTK_ioExample_ascii",numProcs, myRank, PdVTK::vtkASCII);
	PdVTK::CollectionWriter binWriter("utVTK_ioExample_bin",numProcs, myRank, PdVTK::vtkBINARY);

	/*
	 * Write fields
	 * This doesn't actually write until later; Here it just sets the pointers
	 */
	PdVTK::writeField<double>(grid,uSpec,uPtr);
	PdVTK::writeField(grid,vSpec,vPtr);
	PdVTK::writeField(grid,aSpec,aPtr);
	PdVTK::writeField(grid,volSpec,volPtr);
	PdVTK::writeField(grid,wSpec,wPtr);
	PdVTK::writeField(grid,thetaSpec,thetaPtr);
	PdVTK::writeField(grid,rankField);

	/*
	 * Example that loops over time
	 */

	std::size_t numSteps = 10;
	const PdQPointSet1d gammaSpec(numSteps,0,.1);
	double dGamma=gammaSpec.getCellSize();

	/*
	 * Write initial conditions: gamma=0
	 */
	double gamma=0.0;
	asciiWriter.writeTimeStep(gamma,grid);
	binWriter.writeTimeStep(gamma,grid);
	for(std::size_t j=0;j<numSteps;j++){

		/*
		 * Do time integration and physics
		 */
		gamma += dGamma;
		getPureShearXY(gamma, X, uField);


		asciiWriter.writeTimeStep(gamma,grid);
		binWriter.writeTimeStep(gamma,grid);


	}

	/*
	 * This writes the "pvd" collection file
	 */
	asciiWriter.close("mpiexec -np 4 ./utVTK_ioExample_np4\n");
	binWriter.close();



}



bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success=true;
	test_suite* proc = BOOST_TEST_SUITE( "utVTK_ioExample" );
	proc->add(BOOST_TEST_CASE( &utVTK_ioExample ));
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
	 * This test only make sense for numProcs == 4
	 */
	if(4 != numProcs){
		std::cerr << "Unit test runtime ERROR: utVTK_ioExample_np4 only makes sense on 4 processors" << std::endl;
		std::cerr << "\t Re-run unit test $mpiexec -np 4 ./utVTK_ioExample_np4" << std::endl;
		pimpMPI.PimpMpiFixture::~PimpMpiFixture();
		std::exit(-1);
	}

	// Initialize UTF
	return unit_test_main( init_unit_test, argc, argv );
}


