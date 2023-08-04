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

#include "Peridigm_ElasticPlasticMaterial.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include <Epetra_SerialComm.h>
#include "QuickGrid.h"
#include "Field.h"
#include "Array.h"
#include "material_utilities.h"
#include "elastic_plastic.h"
#include <math.h>

#include <fstream>

using namespace MATERIAL_EVALUATION;
using std::shared_ptr;
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

void updateGeometry
(
		const double* xOverlap,
		const double* uOverlap,
		const double* velocityOverlap,
		double* yOverlap,
		int overLapLength,
		double dt
)
{
	const double* x = xOverlap;
	const double* u = uOverlap;
	const double* v = velocityOverlap;
	double*       y = yOverlap;

	int length = overLapLength;
	for(;x!=xOverlap+length;x++,u++,v++,y++)
		*y = *x + *u + *v * dt;
}

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


//

void runPureShear(std::string output_file_name,  TemporalField<double> edpTemporalField, TemporalField<double> lambdaTemporalField, Field_NS::Field<double> fNField, double *x, double *u, double *v, double *y, int numPoints,  double *m, double *theta, double *bondState, double *vol, int *neigh, double K, double MU,  double Y, double *u1x, double *v1x, double *f1x, double DELTA, std::vector<StageFunction> stages, int numStepsPerStage){

     double dt = 1.0/numStepsPerStage;

     /*
      * Create data file
      */


     std::fstream out(output_file_name.c_str(), std::fstream::out);

     /*
	 * Write out initial condition
	 */
	double t=0;
	out << 0 << " " << 0 << " " << 0 << std::endl;
	for(std::vector<StageFunction>::iterator stageIter=stages.begin(); stageIter!=stages.end();stageIter++){

		*v1x = stageIter->slope();
		double vel = *v1x;

		for(int step=0;step<numStepsPerStage;step++){
			Field<double> edpNField = edpTemporalField.getField(Field_ENUM::STEP_N);
			Field<double> edpNP1Field = edpTemporalField.getField(Field_ENUM::STEP_NP1);
			double *edpN = edpNField.get();
			double *edpNP1 = edpNP1Field.get();
			Field<double> lambdaNField = lambdaTemporalField.getField(Field_ENUM::STEP_N);
			Field<double> lambdaNP1Field = lambdaTemporalField.getField(Field_ENUM::STEP_NP1);
			double *lambdaN = lambdaNField.get();
			double *lambdaNP1 = lambdaNP1Field.get();
			fNField.set(0.0);
			double *f = fNField.get();

			t += dt;

			{
				updateGeometry(x,u,v,y,numPoints*3,dt);
			}

			// Dummy arguments for planar case
			bool isPlanarProblem = false;
	        double thickness = 0.0;

			/*
			 * Do not compute dilatation -- just set it to zero
			 */
			computeInternalForceIsotropicElasticPlastic(x,y,m,vol,theta,bondState,edpN,edpNP1,lambdaN,lambdaNP1,f,neigh,numPoints,K,MU,DELTA,Y,isPlanarProblem,thickness);


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

			/*
			 * Next step; this is a bit squirrely -- has to do with updateGeometry
			 * Update geometry takes velocity and existing displacement field (N)
			 * Update geometry y = X + U(N) +Velocity*dt = X + U(NP1)
			 * So, that is why we need this at the bottom of this loop
			 */
			*u1x += vel*dt;

			out << t << " " << *u1x << " " << signF * MAGNITUDE(f1x)
								<< std::endl;
			edpTemporalField.advanceStep();
			lambdaTemporalField.advanceStep();
		}
	}

}

//! Tests state variable count and name accessor functions.
Teuchos::ParameterList getParamList()
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
	double Y = 351.79 / sqrt(3);

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
	Teuchos::ParameterList params;
	params.set("Density", rho);
	params.set("Bulk Modulus", K);
	params.set("Shear Modulus", mu);
	params.set("Horizon", horizon);
	params.set("Yield Stress",Y);
	PeridigmNS::ElasticPlasticMaterial mat(params);

    // \todo check field specs

	return params;
}

QUICKGRID::QuickGridData getTwoPointGridData(){
	int numCells = 2;
	int dimension = 3;
	QUICKGRID::QuickGridData pdGridData = QUICKGRID::allocatePdGridData(numCells,dimension);
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
	 UTILITIES::Array<int> neighborhood(pdGridData.sizeNeighborhoodList);
//	 shared_ptr<int> neighborhood = shared_ptr<int>(new int[pdGridData.sizeNeighborhoodList],PdQuickGrid::Deleter<int>());
	 pdGridData.neighborhood = neighborhood.get_shared_ptr();
	 int *neigh = neighborhood.get();
	 *(neigh+0)=1;
	 *(neigh+1)=1;
	 *(neigh+2)=1;
	 *(neigh+3)=0;

	 UTILITIES::Array<int> neighborhoodPtr(pdGridData.numPoints);
//	 shared_ptr<int> neighborhoodPtr = shared_ptr<int>(new int[pdGridData.numPoints],PdQuickGrid::Deleter<int>());
	 pdGridData.neighborhoodPtr=neighborhoodPtr.get_shared_ptr();
	 int *nPtr = neighborhoodPtr.get();
	 *(nPtr+0)=0;
	 *(nPtr+1)=2;

	 return pdGridData;
}

TEUCHOS_UNIT_TEST(ElasticPlasticMaterial,RunPureShear) {

	Teuchos::ParameterList paramList = getParamList();
	PeridigmNS::ElasticPlasticMaterial mat(paramList);
	QUICKGRID::QuickGridData pdGridData = getTwoPointGridData();
	int numPoints = pdGridData.numPoints;
	TEST_ASSERT(2 == numPoints);
	TEST_ASSERT(4 == pdGridData.sizeNeighborhoodList);

	/*
	 * Material Props
	 */
	double K = paramList.get<double>("Bulk Modulus");
	double MU = paramList.get<double>("Shear Modulus");
	double DELTA = paramList.get<double>("Horizon");
	double Y = paramList.get<double>("Yield Stress");

	/*
	 * yield strain ~.0051 -- 1/2 the engineering strain -- not the same as the
	 * shear strain used here to load
	 */
	double epsYield = .0051;


	/*
	 * Displacement and Internal Force Vectors
	 */
	FieldSpec uSpec = DISPL3D;
	Field<double> uOwnedField(uSpec,pdGridData.numPoints);
	FieldSpec fNSpec = FORCE_DENSITY3D;
	Field_NS::Field<double> fNField(fNSpec,pdGridData.numPoints);
	FieldSpec velocitySpec = VELOC3D;
	Field<double> velField(velocitySpec,pdGridData.numPoints);
	FieldSpec ySpec = CURCOORD3D;
	Field<double> yField(ySpec,pdGridData.numPoints);
	uOwnedField.set(0.0);
	velField.set(0.0);
	fNField.set(0.0);
	yField.set(0.0);
	double *u1x = uOwnedField.get()+3;
	double *v1x = velField.get()+3;
	double *f1x = fNField.get()+3;

	/*
	 * Weighted Volume
	 */
	UTILITIES::Array<double> mPtr(numPoints);
	mPtr.set(0.0);
	MATERIAL_EVALUATION::computeWeightedVolume(pdGridData.myX.get(),pdGridData.cellVolume.get(),mPtr.get(),numPoints,pdGridData.neighborhood.get(),horizon);

	/*
	 * Dilatation
	 */
	UTILITIES::Array<double> thetaPtr(numPoints);
	thetaPtr.set(0.0);

	/*
	 * Bond State and deviatoric plastic extension
	 */
	UTILITIES::Array<double> bondStatePtr(pdGridData.sizeNeighborhoodList-numPoints);
	bondStatePtr.set(0.0);
	TemporalField<double> edpTemporalField(DEVIATORIC_PLASTIC_EXTENSION,pdGridData.sizeNeighborhoodList-numPoints);
	Field<double> edpNField = edpTemporalField.getField(Field_ENUM::STEP_N);
	Field<double> edpNP1Field = edpTemporalField.getField(Field_ENUM::STEP_NP1);
	edpNField.set(0.0);
	TemporalField<double> lambdaTemporalField(Field_NS::LAMBDA,pdGridData.numPoints);
	Field<double> lambdaNField = lambdaTemporalField.getField(Field_ENUM::STEP_N);
	Field<double> lambdaNP1Field = lambdaTemporalField.getField(Field_ENUM::STEP_NP1);
	lambdaNField.set(0.0);
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
	stages[1] = stages[0].next(-.0005);

	/*
	 * Re-Unloading
	 */
	stages[2] = stages[1].next(.001275);

	int numStepsPerStage = 50;
	//double dt = 1.0/numStepsPerStage;

	/*
	 * Pointers to data that don't change in this test
	 */
	double *x = pdGridData.myX.get();
	double *u = uOwnedField.get();
	double *v = velField.get();
	double *y = yField.get();
	double *m = mPtr.get();
	double *theta = thetaPtr.get();
	double *bondState = bondStatePtr.get();
	double *vol = pdGridData.cellVolume.get();
	int *neigh = pdGridData.neighborhood.get();

       runPureShear( "ep.dat", edpTemporalField, lambdaTemporalField, fNField, x, u, v, y,numPoints, m, theta, bondState, vol, neigh, K, MU, Y, u1x, v1x, f1x, DELTA, stages, numStepsPerStage);


}


int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
