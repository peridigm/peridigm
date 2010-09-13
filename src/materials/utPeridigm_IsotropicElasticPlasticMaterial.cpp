/*! \file utPeridigm_LinearElasticIsotropicMaterial.cpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include "Peridigm_IsotropicElasticPlasticMaterial.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Epetra_SerialComm.h>
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "Pd_shared_ptr_Array.h"
#include "PdMaterialUtilities.h"
#include <math.h>

using namespace boost::unit_test;
using namespace std;
using namespace PeridigmNS;
using namespace Teuchos;
using namespace PdQuickGrid;
using namespace PdMaterialUtilities;
using std::tr1::shared_ptr;

const double horizon=sqrt(2);

class Control {
private:
	double start, end;

public:
	Control() : start(0), end(0) {}
	Control(double start, double end) : start(start), end(end){}
	double value(double lambda) const {return (1-lambda)*start + lambda*end;}
	double slope() const { return end-start; }
};

inline void SET(double* const start, const double* const end, double value){
	for(double *i = start; i!=end; i++)
			*i = value;
}

inline double MAGNITUDE(const double *x){
	double x1 = * x;
	double x2 = *(x+1);
	double x3 = *(x+2);
	return sqrt(x1*x1+x2*x2+x3*x3);
}

//! Tests state variable count and name accessor functions.
ParameterList getParamList()
{
	/*
	 * Young's Modulus (MPa)
	 */
	double E = 68.9e3;

	/*
	 * Poisson's ratio
	 */
	double nu = .33;

	/*
	 * Yield "Stress" estimate for perfect plasticity (MPa)
	 * 6061-T6 data
	 */
	double Y = 351.79;

	/*
	 * Density of aluminum g/mm^3
	 */
	double rho = 2.7e-3;

	/*
	 * Bulk Modulus
	 */
	double K = E / 3 / (1-2.0 * nu);

	/*
	 * Shear Modulus
	 */
	double mu = E / 2 / (1+nu);
	ParameterList params;
	params.set("Density", rho);
	params.set("Bulk Modulus", K);
	params.set("Shear Modulus", mu);
	params.set("Material Horizon", horizon);
	params.set("Yield Stress",Y);
	IsotropicElasticPlasticMaterial mat(params);
	BOOST_REQUIRE(mat.NumScalarConstitutiveVariables() == 3);
	BOOST_CHECK(mat.ScalarConstitutiveVariableName(0) == "Weighted Volume");
	BOOST_CHECK(mat.ScalarConstitutiveVariableName(1) == "Dilatation");
	BOOST_CHECK(mat.ScalarConstitutiveVariableName(2) == "Damage");
	BOOST_CHECK_THROW(mat.ScalarConstitutiveVariableName(-1), std::range_error);
	BOOST_CHECK_THROW(mat.ScalarConstitutiveVariableName(3), std::range_error);
	BOOST_REQUIRE(mat.NumVectorConstitutiveVariables() == 1);
	BOOST_CHECK(mat.VectorConstitutiveVariableName(0) == "Current Position");
	BOOST_CHECK_THROW(mat.VectorConstitutiveVariableName(-1), std::range_error);
	BOOST_CHECK_THROW(mat.VectorConstitutiveVariableName(1), std::range_error);
	BOOST_CHECK(mat.BondConstitutiveVariableName(0) == "scalarPlasticExtensionState");
	BOOST_CHECK_THROW(mat.BondConstitutiveVariableName(1), std::range_error);
	BOOST_CHECK(mat.NumBondConstitutiveVariables()==1);

	return params;
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

PdGridData getFourPointGridData(){
	const int nx = 2;
	const int ny = 2;
	const int nz = 1;
	const double xStart = -0.5;
	const double xLength = 2.0;
	const double yStart = -0.5;
	const double yLength = 2.0;
	const double zStart = -0.5;
	const double zLength = 1.0;
	const PdQPointSet1d xSpec(nx,xStart,xLength);
	const PdQPointSet1d ySpec(ny,yStart,yLength);
	const PdQPointSet1d zSpec(nz,zStart,zLength);
	const double horizon(1.1*sqrt(1));
	const int myRank = 0;
	const int numProcs = 1;
	PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numProcs,horizon,xSpec,ySpec,zSpec, SphericalNorm);
	PdGridData decomp =  PdQuickGrid::getDiscretization(myRank, cellPerProcIter);
	return decomp;
}


void runPureShear() {
	ParameterList paramList = getParamList();
	IsotropicElasticPlasticMaterial mat(paramList);
	PdGridData pdGridData = getTwoPointGridData();
	int numPoints = pdGridData.numPoints;
	BOOST_CHECK(2 == numPoints);
	BOOST_CHECK(4 == pdGridData.sizeNeighborhoodList);

	/*
	 * Material Props
	 */
	double K = paramList.get<double>("Bulk Modulus");
	double MU = paramList.get<double>("Shear Modulus");
	double DELTA = paramList.get<double>("Material Horizon");
	double Y = paramList.get<double>("Yield Stress");

	/*
	 * yield strain ~.0051 -- 1/2 the engineering strain -- not the same as the
	 * shear strain used here to load
	 */
	double epsYield = .0051;


	/*
	 * Displacement, Velocity, Current Position, Force
	 */
	Pd_shared_ptr_Array<double> uPtr(numPoints*3), vPtr(numPoints*3), yPtr(numPoints*3), fPtr(numPoints*3);
	SET(uPtr.get(),uPtr.end(),0.0);
	SET(vPtr.get(),vPtr.end(),0.0);
	SET(yPtr.get(),yPtr.end(),0.0);
	SET(fPtr.get(),fPtr.end(),0.0);

	/*
	 * Weighted Volume
	 */
	Pd_shared_ptr_Array<double> mPtr(numPoints);
	SET(mPtr.get(),mPtr.end(),0.0);
	computeWeightedVolume(pdGridData.myX.get(),pdGridData.cellVolume.get(),mPtr.get(),numPoints,pdGridData.neighborhood.get());

	/*
	 * Dilatation
	 */
	Pd_shared_ptr_Array<double> thetaPtr(numPoints);
	SET(thetaPtr.get(),thetaPtr.end(),0.0);

	/*
	 * Bond State and deviatoric plastic extension
	 */
	Pd_shared_ptr_Array<double> bondStatePtr(pdGridData.sizeNeighborhoodList-numPoints), edpNPtr(pdGridData.sizeNeighborhoodList-numPoints);
	SET(bondStatePtr.get(),bondStatePtr.end(),0.0);
	SET(edpNPtr.get(),edpNPtr.end(),0.0);

	/*
	 * Track
	 */
	/*
	 * Stages
	 * 1) load
	 * 2) unload
	 * 3) reload
	 */
	std::vector<Control> stages(3);
	/*
	 * Loading
	 */
	stages[0] = Control(0,epsYield);
	/*
	 * Unloading
	 */
	stages[1] = Control(1.25*epsYield,0);

	/*
	 * Re-Unloading
	 */
	stages[2] = Control(0,0.5*epsYield);

	int numStepsPerStage = 50;
	double dt = 1.0/numStepsPerStage;

	/*
	 * Write out initial condition
	 */
	double t=0;
	std::cout << 0 << " " << 0 << " " << 0 << std::endl;
	for(std::vector<Control>::iterator stageIter=stages.begin(); stageIter!=stages.end();stageIter++){
		double vel = stageIter->slope();
		int p1x = 3;
		double *u1x = uPtr.get()+p1x;
		double *v1x = vPtr.get()+p1x; *v1x=vel;
		double *f1x = fPtr.get()+p1x;

		for(int step=0;step<numStepsPerStage;step++){
			double *x = pdGridData.myX.get();
			double *u = uPtr.get();
			double *v = vPtr.get();
			double *y = yPtr.get();
			double *m = mPtr.get();
			double *theta = thetaPtr.get();
			double *bondState = bondStatePtr.get();
			double *vol = pdGridData.cellVolume.get();
			int *neigh = pdGridData.neighborhood.get();
			double *edpN = edpNPtr.get();
			double *f = fPtr.get();

			updateGeometry(x,u,v,y,numPoints*3,dt);

			/*
			 * Do not compute dilatation -- just set it to zero
			 */
			// computeDilatation(x,y,m,vol,bondState,theta,neigh,numPoints);

			SET(fPtr.get(),fPtr.end(),0.0);
//			computeInternalForceIsotropicElasticPlastic(x,y,m,vol,theta,bondState,edpN,f,neigh,numPoints,K,MU,DELTA,Y);

			/*
			 * Next step
			 */
			*u1x += vel*dt;
			t += dt;

			/*
			 * Get sign of "f" -- this works as long as f does not ever land "exactly" on zero
			 * Put a negative sign in front so that loading is "positive"
			 */
			double signF = -*f1x/abs(*f1x);

			/*
			 * Length of bar
			 * NOTE: original coordinates of point 1 were
			 * x=w, y=w
			 * where w=1.0
			 * l: New length of bar
			 * L: Original length of bar
			 * NOTE: yielding occurs wrt to the shear strain -- not the same as the strain along the axis of the
			 * bond.  This distinction is important to remember.  In this case, the shear strain is equivalent
			 * to the displacement along the x-axis
			 */

			double wx=1;
			double wy=1;
			double l = sqrt((wx+*u1x)*(wx+*u1x)+wy*wy);
			double L = sqrt(wx*wx+wy*wy);
			double engineeringStrain=(l-L)/L;
			double stretch=l/L;
			double trueStrain=std::log(stretch);
			std::cout << t << " " << *u1x << " " << signF*MAGNITUDE(f1x) << std::endl;
		}
	}

}


bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPeridigm_IsotropicElasticPlasticMaterial");
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
  // Initialize UTF
  return unit_test_main(init_unit_test, argc, argv);
}
