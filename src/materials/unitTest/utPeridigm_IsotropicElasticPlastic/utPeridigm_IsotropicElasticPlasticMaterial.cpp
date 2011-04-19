/*! \file utPeridigm_LinearElasticIsotropicMaterial.cpp */

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
#include "Peridigm_IsotropicElasticPlasticMaterial.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Epetra_SerialComm.h>
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "Pd_shared_ptr_Array.h"
#include "PdMaterialUtilities.h"
#include "ordinary_elastic_plastic.h"
#include <math.h>
#include "Field.h"
#include <fstream>

using namespace boost::unit_test;
using namespace std;
using namespace PeridigmNS;
using namespace Teuchos;
using namespace PdQuickGrid;
using namespace PdMaterialUtilities;
using namespace MATERIAL_EVALUATION;
using std::tr1::shared_ptr;
using namespace Field_NS;

const double horizon=sqrt(2);

class StageFunction {

private:
	double start, end;

public:
	StageFunction() : start(0), end(0) {}

	StageFunction(double start, double end) : start(start), end(end){}

	/**
	 * Control function which produces step or proportional loading;
	 * Step loading is given with start==end
	 * Proportional loading is given by start != end
	 * @param lambda Load parameter for stage; it is assumed that 0<=lambda <=1.0
	 * @return Function value
	 */
	double value(double lambda) const {
		return (1-lambda)*start + lambda*end;
	}

	/**
	 * @return slope of function
	 */
	double slope() const { return end-start; }


	/**
	 * StageFunction for next stage
	 * @return new StageFunction will hold constant
	 */
	StageFunction next() const {
		double fStartNew = this->end;
		double fEndNew = fStartNew;
		return StageFunction(fStartNew,fEndNew);
	}



	/**
	 * StageFunction for next stage
	 * @param <code>endVal </code>value of loader at end of load step
	 * @return new proportional<code>StageFunction</code> with starting value of this end value and end value <code>endVal</code>
	 */
	StageFunction next(double endVal) const {
		double fStartNew = this->end;
		double fEndNew = endVal;
		return StageFunction(fStartNew,fEndNew);
	}

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

    // \todo check field specs

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
	 * Displacement and Internal Force Vectors
	 */
	FieldSpec uSpec(FieldSpec::DISPLACEMENT,FieldSpec::VECTOR3D,"u");
	Field<double> uOwnedField(uSpec,pdGridData.numPoints);
	FieldSpec fNSpec(FieldSpec::FORCE,FieldSpec::VECTOR3D,"fN");
	Field_NS::Field<double> fNField(fNSpec,pdGridData.numPoints);
	const FieldSpec velocitySpec(FieldSpec::VELOCITY,FieldSpec::VECTOR3D, "velocity");
	Field<double> velField(velocitySpec,pdGridData.numPoints);
	const FieldSpec ySpec(FieldSpec::COORDINATES,FieldSpec::VECTOR3D, "CURRENT_COORDINATES");
	Field<double> yField(ySpec,pdGridData.numPoints);
	const FieldSpec dsfSpec = Field_NS::SHEAR_CORRECTION_FACTOR;
	Field<double> dsfField(dsfSpec,pdGridData.numPoints);
	uOwnedField.setValue(0.0);
	velField.setValue(0.0);
	fNField.setValue(0.0);
	yField.setValue(0.0);
	dsfField.setValue(1.0);
	double *u1x = uOwnedField.getArray().get()+3;
	double *v1x = velField.getArray().get()+3;
	double *f1x = fNField.getArray().get()+3;

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
	Pd_shared_ptr_Array<double> bondStatePtr(pdGridData.sizeNeighborhoodList-numPoints);
	SET(bondStatePtr.get(),bondStatePtr.end(),0.0);
	TemporalField<double> edpTemporalField(DEVIATORIC_PLASTIC_EXTENSION,pdGridData.sizeNeighborhoodList-numPoints);
	Field<double> edpNField = edpTemporalField.getField(FieldSpec::STEP_N);
	Field<double> edpNP1Field = edpTemporalField.getField(FieldSpec::STEP_NP1);
	SET(edpNField.getArray().get(),edpNField.getArray().end(),0.0);
	TemporalField<double> lambdaTemporalField(Field_NS::LAMBDA,pdGridData.numPoints);
	Field<double> lambdaNField = lambdaTemporalField.getField(FieldSpec::STEP_N);
	Field<double> lambdaNP1Field = lambdaTemporalField.getField(FieldSpec::STEP_NP1);
	SET(lambdaNField.getArray().get(),lambdaNField.getArray().end(),0.0);
	/*
	 * Track
	 */
	/*
	 * Stages
	 * 1) load
	 * 2) unload
	 * 3) reload
	 */
	std::vector<StageFunction> stages(3);
	/*
	 * Loading
	 */
	stages[0] = StageFunction(0,epsYield);
	/*
	 * Unloading
	 */
	stages[1] = stages[0].next(-.001275);

	/*
	 * Re-Unloading
	 */
	stages[2] = stages[1].next(.001275);

	int numStepsPerStage = 50;
	double dt = 1.0/numStepsPerStage;

	/*
	 * Pointers to data that don't change in this test
	 */
	double *x = pdGridData.myX.get();
	double *u = uOwnedField.getArray().get();
	double *v = velField.getArray().get();
	double *y = yField.getArray().get();
	double *m = mPtr.get();
	double *theta = thetaPtr.get();
	double *bondState = bondStatePtr.get();
	double* dsfOwned = dsfField.getArray().get();
	double *vol = pdGridData.cellVolume.get();
	int *neigh = pdGridData.neighborhood.get();

	/*
	 * Create data file
	 */
	std::fstream out("ep.dat", std::fstream::out);

	/*
	 * Write out initial condition
	 */
	double t=0;
	out << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
	for(std::vector<StageFunction>::iterator stageIter=stages.begin(); stageIter!=stages.end();stageIter++){

		double vel = stageIter->slope();
		*v1x = vel;


		for(int step=0;step<numStepsPerStage;step++){
			Field<double> edpNField = edpTemporalField.getField(FieldSpec::STEP_N);
			Field<double> edpNP1Field = edpTemporalField.getField(FieldSpec::STEP_NP1);
			double *edpN = edpNField.getArray().get();
			double *edpNP1 = edpNP1Field.getArray().get();
			Field<double> lambdaNField = lambdaTemporalField.getField(FieldSpec::STEP_N);
			Field<double> lambdaNP1Field = lambdaTemporalField.getField(FieldSpec::STEP_NP1);
			double *lambdaN = lambdaNField.getArray().get();
			double *lambdaNP1 = lambdaNP1Field.getArray().get();
			fNField.setValue(0.0);
			double *f = fNField.getArray().get();

			t += dt;

			updateGeometry(x,u,v,y,numPoints*3,dt);

			/*
			 * Do not compute dilatation -- just set it to zero
			 */
			computeInternalForceIsotropicElasticPlastic(x,y,m,vol,theta,bondState,dsfOwned,edpN,edpNP1,lambdaN,lambdaNP1,f,neigh,numPoints,K,MU,DELTA,Y);


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
			/*
			 * Next step; this is a bit squirely -- has to do with updateGeometry
			 * Update geometry takes velocity and existing displacement field (N)
			 * Update geometry y = X + U(N) +Velocity*dt = X + U(NP1)
			 * So, that is why we need this at the bottom of this loop
			 */
			*u1x += vel*dt;

			out << t << " " << *u1x << " " << *lambdaNP1 << " "<< signF*MAGNITUDE(f1x) << std::endl;
			edpTemporalField.advanceStep();
			lambdaTemporalField.advanceStep();
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
