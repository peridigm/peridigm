#include "PdImpMpiFixture.h"
#include "PdITI_Operator.h"
#include "NeighborhoodList.h"
#include "../mesh_input/quick_grid/QuickGrid.h"
#include "../mesh_input/quick_grid/QuickGridData.h"
#include "PdVTK.h"
#include "Field.h"
#include "PdZoltan.h"
#include "Array.h"
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
using PdITI::IsotropicElasticConstitutiveModel;
using PdITI::ConstitutiveModel;


/*
 * Following dataset produces 255 cells (2d) = ring2dSpec.getNumCells()
 * 	 double innerRadius = 0.020;
 *   double outerRadius = 0.025;
 *   int numRings = 3;
 */
const double innerRadius = 20.0;
const double outerRadius = 25.0;
const double a = 50.0;
const double cylinderLength = 2.0*a;
const int numRings = 3;

/*
 * Create 2d Ring
 */
const double xC = 0.0;
const double yC = 0.0;
const double zC = 0.0;
const double centerPoint[] = {xC, yC, zC} ;
const std::valarray<double> center(centerPoint,3);

/*
 * Note that zStart is used for the 1D spec along cylinder axis
 */
const double zStart = 0.0;
const QUICKGRID::SpecRing2D ring2dSpec(center,innerRadius,outerRadius,numRings);

/*
 * Create 1d Spec along cylinder axis
 * Compute number of cells along length of cylinder so that aspect ratio
 * is cells is approximately 1.
 * Cell sizes along axis are not exactly "cellSize" since last cell
 * would be a fraction of a cellSize -- so 1 is added to numCellsAlongAxis.
 * Actual cell sizes are slightly different than "cellSize" because of this.
 */

const double SCALE=2.51;
const double cellSize = ring2dSpec.getRaySpec().getCellSize();
const double horizon = SCALE*cellSize;
const size_t numCellsAxis = (int)(cylinderLength/cellSize)+1;
const QUICKGRID::Spec1D axisSpec(numCellsAxis,zStart,cylinderLength);



shared_ptr<double> getAxialExtensionZ(size_t numOwnedPoints, shared_ptr<double>& xPtr);
QUICKGRID::QuickGridData getCylinderDiscretizaton(int rank, int numProcs);
size_t getDemoDiscretizaton(size_t size_t, int numProcs);
Field_NS::TemporalField<double> getOwnedQ2CylinderInitialConditions(double vr0, double vr1, double vz0, double a, QUICKGRID::QuickGridData& gridData);

int main( int argc, char *argv[]) {
	PdImpRunTime::PimpMpiFixture pimpMPI = PdImpRunTime::PimpMpiFixture::getPimpMPI(argc,argv);
	const Epetra_Comm& comm = pimpMPI.getEpetra_Comm();
	size_t decomp = getCylinderDiscretizaton(comm.MyPID(), comm.NumProc());
//	PdGridData decomp = getDemoDiscretizaton(comm.MyPID(), comm.NumProc());

	/*
	 * Load balance and write new decomposition
	 */
	decomp=QUICKGRID::getLoadBalancedDiscretization(decomp);
	int numPoints = decomp.numPoints;
	PDNEIGH::NeighborhoodList list(comm,decomp.zoltanPtr.get(),numPoints,decomp.myGlobalIDs,decomp.myX,horizon);
	/*
	 * Create an operator
	 */
	const PdImp::BulkModulus _K(130000.0);
	const PdImp::PoissonsRatio _MU(0.25);
	const PdImp::IsotropicHookeSpec isotropicSpec(_K,_MU);
	const PdImp::MassDensity rho(1.0);
	std::tr1::shared_ptr<PdITI::PdITI_Operator> op(new PdITI::PdITI_Operator(comm,list,decomp.cellVolume));
	shared_ptr<ConstitutiveModel> fIntOperator(new IsotropicElasticConstitutiveModel(isotropicSpec));
	op->addConstitutiveModel(fIntOperator);


	/*
	 * Create VTK Grid for Output
	 */

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
	vr0 = 0;
	vr1 = 0;


	//		std::tr1::shared_ptr<double> uPtr = getPureShear();
	//		std::tr1::shared_ptr<double> uPtr = getAxialExtensionZ(decomp.numPoints,decomp.myX);


	/*
	 * Create force vector and initialize to zero
	 */
	Field_NS::TemporalField<double> displacement = getOwnedQ2CylinderInitialConditions(vr0,vr1,vz0,a,decomp);
	Field_NS::TemporalField<double> force = Field_NS::TemporalField<double>(Field_NS::FORCE3D,numPoints);
	Field_NS::Field<double> theta = op->computeOwnedDilatation(displacement.getField(Field_NS::FieldSpec::STEP_NP1));
	op->computeInternalForce(displacement.getField(Field_NS::FieldSpec::STEP_NP1),force.getField(Field_NS::FieldSpec::STEP_NP1));
	/*
	 * this puts the displacement in STEP_N
	 */
	displacement.advanceStep();
	writeField<double>(grid, displacement.getField(Field_NS::FieldSpec::STEP_N));
	writeField<double>(grid, theta);
	writeField<double>(grid, force.getField(Field_NS::FieldSpec::STEP_NP1));
	PdVTK::CollectionWriter collectionWriter("Q2CylinderBalanced",comm.NumProc(), comm.MyPID());
	collectionWriter.writeTimeStep(0.0, grid);
	collectionWriter.close();

	return 0;

}

QUICKGRID::QuickGridData getCylinderDiscretizaton(size_t myRank, size_t numProcs){
	// Create abstract decomposition iterator
	QUICKGRID::TensorProductCylinderMeshGenerator cellPerProcIter(numProcs, horizon,ring2dSpec, axisSpec, PdQuickGrid::SphericalNorm);
	QUICKGRID::QuickGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	return decomp;
}

Field_NS::TemporalField<double> getOwnedQ2CylinderInitialConditions(double vr0, double vr1, double vz0, double a, QUICKGRID::QuickGridData& gridData){
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

UTILITIES::Array<double> getAxialExtensionZ(int numOwnedPoints, shared_ptr<double>& xPtr){
	UTILITIES::Array<double> uPtr(3*numOwnedPoints);
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


QUICKGRID::QuickGridData getDemoDiscretizaton(size_t myRank, size_t numProcs){
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
	const QUICKGRID::Spec1D xSpec(nx,xStart,xLength);
	const QUICKGRID::Spec1D ySpec(ny,yStart,yLength);
	const QUICKGRID::Spec1D zSpec(nz,zStart,zLength);
	double horizon = .1;
	QUICKGRID::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec);
	return QUICKGRID::getDiscretization(myRank, cellPerProcIter);
}

