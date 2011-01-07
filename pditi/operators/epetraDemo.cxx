#include "PdImpMpiFixture.h"
#include "PdImpOperator.h"
#include "PdGridData.h"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdVTK.h"
#include "Field.h"
#include "PdZoltan.h"
#include "StageFunction.h"
#include "IsotropicElasticConstitutiveModel.h"
#include "ConstitutiveModel.h"
#include "ImplicitDynamicsIntegrator.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include <tr1/memory>
#include <iostream>
using std::cout;
using std::endl;
using PdVTK::writeField;
using PdImp::IsotropicElasticConstitutiveModel;
using PdImp::ConstitutiveModel;

shared_ptr<double> getAxialExtensionZ(int numOwnedPoints, shared_ptr<double>& xPtr);
PdGridData getCylinderDiscretizaton(int rank, int numProcs);
PdGridData getDemoDiscretizaton(int myRank, int numProcs);
Field_NS::TemporalField<double> getOwnedQ2CylinderInitialConditions(double vr0, double vr1, double vz0, double a, PdGridData& gridData);

int main( int argc, char *argv[]) {
	PdImpRunTime::PimpMpiFixture pimpMPI = PdImpRunTime::PimpMpiFixture::getPimpMPI(argc,argv);
	const Epetra_Comm& comm = pimpMPI.getEpetra_Comm();
	PdGridData decomp = getCylinderDiscretizaton(comm.MyPID(), comm.NumProc());
//	PdGridData decomp = getDemoDiscretizaton(comm.MyPID(), comm.NumProc());
	{
		/*
		 * Load balance and write new decomposition
		 */
		decomp=getLoadBalancedDiscretization(decomp);
		/*
		 * Create an operator
		 */
		const PdImp::BulkModulus _K(130000.0);
		const PdImp::PoissonsRatio _MU(0.25);
		const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
		const PdImp::MassDensity rho(1.0);
		std::tr1::shared_ptr<PdImp::PdImpOperator>  op(new PdImp::PdImpOperator(comm,decomp));
		shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
		op->addConstitutiveModel(fIntOperator);


		/*
		 * Create VTK Grid for Output
		 */
		int numPoints = decomp.numPoints;
		vtkSmartPointer<vtkUnstructuredGrid> grid = PdVTK::getGrid(decomp.myX,numPoints);

		/*
		 * Add processor rank to output
		 */
		Field_NS::Field<int> myRank(Field_NS::PROC_NUM,decomp.numPoints);
		myRank.setValue(comm.MyPID());
		writeField<int>(grid,myRank);

		/*
		 * These fields are independent of the displacement
		 */
		Field_NS::Field<double> vol = Field_NS::getVOLUME(decomp.cellVolume,numPoints);
		Field_NS::Field<double> m   = op->getWeightedVolume();
		writeField<double>(grid,vol);
		writeField<double>(grid,m);

		/*
		 * Set a displacement field using prescribed initial velocity
		 */
		double vr0 = .005;
		double vr1 = .005;
		double vz0 = .5;
		double a = 50.0;
		vr0 = 0;
		vr1 = 0;


//		std::tr1::shared_ptr<double> uPtr = getPureShear();
//		std::tr1::shared_ptr<double> uPtr = getAxialExtensionZ(decomp.numPoints,decomp.myX);


		/*
		 * Create force vector and initialize to zero
		 */
		Field_NS::TemporalField<double> displacement = getOwnedQ2CylinderInitialConditions(vr0,vr1,vz0,a,decomp);
		Field_NS::TemporalField<double> velocity = Field_NS::getVELOC3D(numPoints);
		Field_NS::TemporalField<double> force = Field_NS::TemporalField<double>(Field_NS::FORCE3D,numPoints);

		Field_NS::Field<double> theta = op->computeOwnedDilatation(displacement.getField(Field_NS::FieldSpec::STEP_NP1));
		op->computeInternalForce(displacement.getField(Field_NS::FieldSpec::STEP_NP1),force.getField(Field_NS::FieldSpec::STEP_NP1));
//		residualOp->computeResidual(displacement,velocity,0.0,force);
		displacement.advanceStep();
		writeField<double>(grid, displacement.getField(Field_NS::FieldSpec::STEP_N));
		writeField<double>(grid, theta);
		writeField<double>(grid, force.getField(Field_NS::FieldSpec::STEP_NP1));
		vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer = PdVTK::getWriter("Q2CylinderBalanced.pvtu", comm.NumProc(), comm.MyPID());
		PdVTK::write(writer,grid);
		PdVTK::CollectionWriter collectionWriter("Q2CylinderBalanced",comm.NumProc(), comm.MyPID());
		collectionWriter.writeTimeStep(0.0, grid);
		/*
		 * HERE YOU create a new grid with new values for the coordinates and for field variables
		 */
		collectionWriter.writeTimeStep(0.6667, grid);
		collectionWriter.close();

	}
	return 0;

}

PdGridData getCylinderDiscretizaton(int myRank, int numProcs){
	/*
	 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
	 * 	 double innerRadius = 0.020;
	 *   double outerRadius = 0.025;
	 *   int numRings = 3;
	 */
	double innerRadius = 20.0;
	double outerRadius = 25.0;
	double a = 50.0;
	double cylinderLength = 2.0*a;
	int numRings = 3;

	/*
	 * Create 2d Ring
	 */
	double xC = 0.0;
	double yC = 0.0;
	std::valarray<double> center(0.0,3);
	center[0] = xC;
	center[1] = yC;
	center[2] = 0;
	/*
	 * Note that zStart is used for the 1D spec along cylinder axis
	 */
	double zStart = 0.0;
	PdQuickGrid::PdQRing2d ring2dSpec(center,innerRadius,outerRadius,numRings);

	/*
	 * Create 1d Spec along cylinder axis
     * Compute number of cells along length of cylinder so that aspect ratio
	 * is cells is approximately 1.
	 * Cell sizes along axis are not exactly "cellSize" since last cell
	 * would be a fraction of a cellSize -- so 1 is added to numCellsAlongAxis.
	 * Actual cell sizes are slightly different than "cellSize" because of this.
	 */

	double SCALE=2.51;
	double cellSize = ring2dSpec.getRaySpec().getCellSize();
	double horizon = SCALE*cellSize;
	int numCellsAxis = (int)(cylinderLength/cellSize)+1;
	PdQuickGrid::PdQPointSet1d axisSpec(numCellsAxis,zStart,cylinderLength);

	// Create abstract decomposition iterator
	PdQuickGrid::TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec, PdQuickGrid::SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	return decomp;
}

Field_NS::TemporalField<double> getOwnedQ2CylinderInitialConditions(double vr0, double vr1, double vz0, double a, PdGridData& gridData){
	int numPoints = gridData.numPoints;
	Field_NS::TemporalField<double> v0(Field_NS::getDISPL3D(numPoints));
	double *v = v0.getField(Field_NS::FieldSpec::STEP_NP1).getArray().get();
	double *x = gridData.myX.get();
	for(int i=0;i<numPoints;i++,x+=3){
		double zBya = (*(x+2))/a;
		double vr = vr0 - vr1 * (zBya)*(zBya);
		double vz = vz0 * zBya;
		double xx = *x;
		double yy = *(x+1);
		double r = sqrt(xx*xx+yy*yy);
		double c = xx/r;
		double s = yy/r;
		*v = vr*c; v++;
		*v = vr*s; v++;
		*v = vz; v++;
	}
	return v0;
}

shared_ptr<double> getAxialExtensionZ(int numOwnedPoints, shared_ptr<double>& xPtr){
	shared_ptr<double> uPtr(new double[3*numOwnedPoints],PdQuickGrid::Deleter<double>());
	double *u = uPtr.get();
	double *x = xPtr.get();
	double gamma=1.0;
	for(int i=0;i<numOwnedPoints;i++,x+=3,u+=3){
		u[0]=0;
		u[1]=0;
		u[2]=gamma*x[2];
	}
	return uPtr;
}

void printNeighborList(std::tr1::shared_ptr<int> localNeighborList,int numOwned){
	int *neighPtr = localNeighborList.get();
	for(int n=0;n<numOwned;n++){
		int numNeig = *neighPtr; neighPtr++;
		for(int m=0;m<numNeig;m++,neighPtr++){
			if(0==m%15)
				cout << endl;
			cout << *neighPtr << ", ";
		}
		cout << endl;
	}
}


PdGridData getDemoDiscretizaton(int myRank, int numProcs){
	const int nx = 1;
	const int ny = 1;
	const int nz = 100;
	const double xStart = 0.0;
	const double xLength = 1.0;
	const double yStart = 0.0;
	const double yLength = 1.0;
	const double zStart = 0.0;
	const double zLength = 1.0;
//	const double zLength = nz*xLength/nx; // this makes cell size along z -- same as x, y
	const PdQuickGrid::PdQPointSet1d xSpec(nx,xStart,xLength);
	const PdQuickGrid::PdQPointSet1d ySpec(ny,yStart,yLength);
	const PdQuickGrid::PdQPointSet1d zSpec(nz,zStart,zLength);
	double horizon = .1;
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	return PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
}

