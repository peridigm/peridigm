/*! \file twoPoint_SLS_Relaxation.cpp */

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

#include "Peridigm_ViscoelasticMaterial.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include <Epetra_SerialComm.h>
#include "QuickGrid.h"
#include "Field.h"
#include "Array.h"
#include "material_utilities.h"
#include "elastic.h"
#include "viscoelastic.h"
#include <math.h>

#include <fstream>

using namespace MATERIAL_EVALUATION;
using std::shared_ptr;
using namespace Field_NS;

/*
 * Young's Modulus (MPa)
 */
const double E = 68.9e3;

/*
 * Poisson's ratio
 */
const double nu = .33;

/*
 * Applied shear and Initial Condition
 */
const double my_gamma = 1.0e-6;


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

Teuchos::ParameterList getParamList(double lambda)
{
	/*
	 * Horizon for this problem
	 */

	double horizon=1.01 * sqrt(2);

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

	/*
	 * Material time constant: Relaxation time
	 */
	double tau_b = 2.0;

	/*
	 * Relaxation Modulus
	 */
	double lambda_i = lambda;

	Teuchos::ParameterList params;

	params.set("Density", rho);
	params.set("Bulk Modulus", K);
	params.set("Shear Modulus", mu);
	params.set("Horizon", horizon);
	params.set("lambda_i",lambda_i);
	params.set("tau b",tau_b);
	PeridigmNS::ViscoelasticMaterial mat(params);

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
	 Array<int> neighborhood(pdGridData.sizeNeighborhoodList);
	 pdGridData.neighborhood = neighborhood.get_shared_ptr();
	 int *neigh = neighborhood.get();
	 *(neigh+0)=1;
	 *(neigh+1)=1;
	 *(neigh+2)=1;
	 *(neigh+3)=0;

	 Array<int> neighborhoodPtr(pdGridData.numPoints);
	 pdGridData.neighborhoodPtr=neighborhoodPtr.get_shared_ptr();
	 int *nPtr = neighborhoodPtr.get();
	 *(nPtr+0)=0;
	 *(nPtr+1)=2;

	 return pdGridData;
}

void updateGeometry
(
		const double* xOverlap,
		const double* uOverlap,
		double* yOverlap,
		int overLapLength
)
{
	const double* x = xOverlap;
	const double* u = uOverlap;
	double*       y = yOverlap;

	int length = overLapLength;
	for(;x!=xOverlap+length;x++,u++,y++)
		*y = *x + *u;
}

// Set material properties

void setMaterialProperties( Teuchos::ParameterList& paramList, double &K, double &MU, double &horizon, double &lambda_i, double &tau_b){

         K = paramList.get<double>("Bulk Modulus");
	 MU = paramList.get<double>("Shear Modulus");
	 horizon = paramList.get<double>("Horizon");
	 lambda_i = paramList.get<double>("lambda_i");
	 tau_b = paramList.get<double>("tau b");
}


// Data initialization
     

void dataInitialization( QUICKGRID::QuickGridData& pdGridData, int numOwnedPoints, double horizon, Field<double> &uOwnedField, Field_NS::Field<double> &fNField, UTILITIES::Array<double> &mPtr, TemporalField<double> &dilatationTemporalField, UTILITIES::Array<double> &bondDamagePtr, TemporalField<double> &edbTemporalField){

        /*
	 * Displacement and Internal Force Vectors
	 */

        uOwnedField.set(0.0);
	fNField.set(0.0);

        /*
	 * Weighted Volume initialization
	 */

        mPtr.set(0.0);
	MATERIAL_EVALUATION::computeWeightedVolume(pdGridData.myX.get(),pdGridData.cellVolume.get(),mPtr.get(),numOwnedPoints,pdGridData.neighborhood.get(),horizon);

        /*
	 * Dilatation: intialize to zero
	 */

        {
		Field<double> N = dilatationTemporalField.getField(Field_ENUM::STEP_N);
		N.set(0.0);
	}

        /*
	 * Bond State and deviatoric back extension state
	 */
   
        bondDamagePtr.set(0.0);

        {
		Field<double> N = edbTemporalField.getField(Field_ENUM::STEP_N);
		N.set(0.0);
	}
	

}

/**
 * Function returns value of force at last time step
 */

double runPureShear(std::string output_file_name, double K, double MU, double lambda_i, double tau_b, int numOwnedPoints, int numSteps_stage_1, double t_start, double dt, TemporalField<double> yTemporalField, Field_NS::Field<double> fNField, TemporalField<double> dilatationTemporalField, TemporalField<double> edbTemporalField, double *u1x, double *f1x, double *xOverlapPtr, double *yN_OverlapPtr, double *yNp1_OverlapPtr, double *mOwned, double *volumeOverlapPtr, double *bondDamage, int *localNeighborList, double *fInternalOverlapPtr){

        double t=t_start;
	double signF = -*f1x/fabs(*f1x);

        /*
	 * Create data file
	 */


        std::fstream out(output_file_name.c_str(), std::fstream::out);
	out << std::scientific;
	out.precision(2);
	out << t << " ";
	out.precision(0);
	out << *u1x << " ";
	out.precision(15);
	out << signF*MAGNITUDE(f1x) << std::endl;

	

	for(size_t n=0;n<static_cast<unsigned int>(numSteps_stage_1);n++){

		t += dt;

		/*
		 * Do not compute dilatation -- just set it to zero
		 */
		Field<double> thetaNField = dilatationTemporalField.getField(Field_ENUM::STEP_N);
		Field<double> thetaNP1Field = dilatationTemporalField.getField(Field_ENUM::STEP_NP1);
		thetaNP1Field.set(0.0);
		const double* dilatationN_Owned   = thetaNField.get();
		double* dilatationNp1_Owned       = thetaNP1Field.get();

		Field<double> edbNField = edbTemporalField.getField(Field_ENUM::STEP_N);
		Field<double> edbNP1Field = edbTemporalField.getField(Field_ENUM::STEP_NP1);
		const double* deviatoricBackExtensionState_N   = edbNField.get();
		double* deviatoricBackExtensionState_Np1       = edbNP1Field.get();

		/*
		 * initialize force to zero
		 */
		fNField.set(0.0);

		MATERIAL_EVALUATION::computeInternalForceViscoelasticStandardLinearSolid
		(
				dt,
				xOverlapPtr,
				yN_OverlapPtr,
				yNp1_OverlapPtr,
				mOwned,
				volumeOverlapPtr,
				dilatationN_Owned,
				dilatationNp1_Owned,
				bondDamage,
				deviatoricBackExtensionState_N,
				deviatoricBackExtensionState_Np1,
				fInternalOverlapPtr,
				localNeighborList,
				numOwnedPoints,
				K,
				MU,
				lambda_i,
				tau_b
		);

                 
		/*
		 * Get sign of "f" -- this works as long as f does not ever land "exactly" on zero
		 * Put a negative sign in front so that loading is "positive"
		 */
		double signF = -*f1x/abs(*f1x);
		out.precision(2);
		out << t << " ";
		out.precision(0);
		out << *u1x << " ";
		out.precision(15);
		out << signF*MAGNITUDE(f1x) << std::endl;

		edbTemporalField.advanceStep();
		dilatationTemporalField.advanceStep();
		yTemporalField.advanceStep();

	}

	out.close();

	/**
	 * this is value of force for last time step
	 */
      
       return MAGNITUDE(f1x);

}

TEUCHOS_UNIT_TEST(TwoPoint_SLS_Relaxation,Case_1) {

	/*
	 * This produces the 'elastic' response
	 */

	double scale=0.01;
	Teuchos::ParameterList paramList = getParamList(scale);
	PeridigmNS::ViscoelasticMaterial mat(paramList);

        QUICKGRID::QuickGridData pdGridData = getTwoPointGridData();
	int numOwnedPoints = pdGridData.numPoints;
	TEST_ASSERT(2 == numOwnedPoints);
	TEST_ASSERT(4 == pdGridData.sizeNeighborhoodList);

        /*
	 * Material Props
	 */

        double K, MU, horizon, lambda_i, tau_b;

        setMaterialProperties( paramList, K, MU, horizon, lambda_i, tau_b);

        /*
	 * Time stepping data
	 */
	double t_start=0.0;
	double t_end = 3.0 * tau_b;
	size_t numSteps_stage_1(100);

	double dt = (t_end - t_start) / numSteps_stage_1;

        /*
	 * Displacement and Internal Force Vectors
	 */

	FieldSpec uSpec = DISPL3D;
	Field<double> uOwnedField(uSpec,pdGridData.numPoints);
	FieldSpec fNSpec = FORCE_DENSITY3D;
	Field_NS::Field<double> fNField(fNSpec,pdGridData.numPoints);

        /*
	 * Weighted Volume
	 */
	UTILITIES::Array<double> mPtr(numOwnedPoints);
	
	
	/*
	 * Dilatation: intialize to zero
	 */
	TemporalField<double> dilatationTemporalField(DILATATION,numOwnedPoints);
	
	/*
	 * Bond State and deviatoric back extension state
	 */
	UTILITIES::Array<double> bondDamagePtr(pdGridData.sizeNeighborhoodList-numOwnedPoints);
	
	TemporalField<double> edbTemporalField(DEVIATORIC_BACK_EXTENSION,pdGridData.sizeNeighborhoodList-numOwnedPoints);

        // Data initialization

        dataInitialization( pdGridData, numOwnedPoints, horizon, uOwnedField, fNField, mPtr,  dilatationTemporalField, bondDamagePtr, edbTemporalField);

        double *u1x = uOwnedField.get()+ 3;
	double *f1x = fNField.get()+ 3;

	/*
	 * INITIAL CONDITION
	 */
	*u1x = my_gamma;

        /*
	 * Pointers to data that don't change in this test
	 */

        double *xOverlapPtr = pdGridData.myX.get();
	double *mOwned = mPtr.get();
	double *bondDamage = bondDamagePtr.get();
	double *volumeOverlapPtr = pdGridData.cellVolume.get();
	int *localNeighborList = pdGridData.neighborhood.get();
	double *fInternalOverlapPtr = fNField.get();

        /*
	 * Current coordinates are fixed; Assign here
	 */

	TemporalField<double> yTemporalField(CURCOORD3D,numOwnedPoints);
	updateGeometry(xOverlapPtr,uOwnedField.get(),yTemporalField.getField(Field_ENUM::STEP_N).get(),numOwnedPoints*3);
	updateGeometry(xOverlapPtr,uOwnedField.get(),yTemporalField.getField(Field_ENUM::STEP_NP1).get(),numOwnedPoints*3);
	double *yN_OverlapPtr = yTemporalField.getField(Field_ENUM::STEP_N).get();
	double *yNp1_OverlapPtr = yTemporalField.getField(Field_ENUM::STEP_NP1).get();
        
	double* nullPtr = NULL;

        /*
	 * Compute initial force with elastic constitutive model
	 */
	MATERIAL_EVALUATION::computeInternalForceLinearElastic
	(
			xOverlapPtr,
			yN_OverlapPtr,
			mOwned,
			volumeOverlapPtr,
			dilatationTemporalField.getField(Field_ENUM::STEP_N).get(),
			bondDamage,
			fInternalOverlapPtr,
			nullPtr,
			localNeighborList,
			numOwnedPoints,
			K,
			MU,
            horizon
	);

        /*
	 * this is the correct analytical value for the initial condition
	 * ??? Why can't we get this closer than 1.0e-6 ???
	 */

	double f0 = 2.0 * 15.0 * E * my_gamma / 4.0 / (1+nu) / std::sqrt(2.0);
	double tolerance=1.0e-6;
	double rel_diff = std::fabs(f0-MAGNITUDE(f1x))/f0;
	
        TEST_COMPARE( rel_diff, <=, tolerance);

        

	double f = runPureShear("twoPoint_SLS_Elastic.dat", K, MU, lambda_i, tau_b, numOwnedPoints, numSteps_stage_1, t_start, dt, yTemporalField, fNField, dilatationTemporalField,  edbTemporalField, u1x, f1x,xOverlapPtr, yN_OverlapPtr, yNp1_OverlapPtr, mOwned, volumeOverlapPtr, bondDamage, localNeighborList, fInternalOverlapPtr);

	
	/*
	 * Last value computed: tests time integrator against exact value
	 */


	/*
	 * step input deviatoric extension state
	 */

	double ed0=1.0e-6/sqrt(2);
	double m = 2.0;
	double alpha = 15.0 * E / (1+nu) / 2.0 / m;
	double tEnd=6.0;
	double fEnd = 2.0 * ed0*((1-lambda_i) * alpha +lambda_i * alpha * std::exp(-tEnd/tau_b));
	rel_diff = std::fabs(fEnd-f)/fEnd;
	
        TEST_COMPARE( rel_diff, <=, tolerance);

}

TEUCHOS_UNIT_TEST(TwoPoint_SLS_Relaxation,Case_2) {

	/*
	 * This produces the 'elastic' response
	 */

	double scale=0.5;
	Teuchos::ParameterList paramList = getParamList(scale);
	PeridigmNS::ViscoelasticMaterial mat(paramList);

        QUICKGRID::QuickGridData pdGridData = getTwoPointGridData();
	int numOwnedPoints = pdGridData.numPoints;
	TEST_ASSERT(2 == numOwnedPoints);
	TEST_ASSERT(4 == pdGridData.sizeNeighborhoodList);

        /*
	 * Material Props
	 */

        double K, MU, horizon, lambda_i, tau_b;

        setMaterialProperties( paramList, K, MU, horizon, lambda_i, tau_b);

        /*
	 * Time stepping data
	 */
	double t_start=0.0;
	double t_end = 3.0 * tau_b;
	size_t numSteps_stage_1(100);

	double dt = (t_end - t_start) / numSteps_stage_1;

        /*
	 * Displacement and Internal Force Vectors
	 */

	FieldSpec uSpec = DISPL3D;
	Field<double> uOwnedField(uSpec,pdGridData.numPoints);
	FieldSpec fNSpec = FORCE_DENSITY3D;
	Field_NS::Field<double> fNField(fNSpec,pdGridData.numPoints);

        /*
	 * Weighted Volume
	 */
	UTILITIES::Array<double> mPtr(numOwnedPoints);
	
	
	/*
	 * Dilatation: intialize to zero
	 */
	TemporalField<double> dilatationTemporalField(DILATATION,numOwnedPoints);
	
	/*
	 * Bond State and deviatoric back extension state
	 */
	UTILITIES::Array<double> bondDamagePtr(pdGridData.sizeNeighborhoodList-numOwnedPoints);
	
	TemporalField<double> edbTemporalField(DEVIATORIC_BACK_EXTENSION,pdGridData.sizeNeighborhoodList-numOwnedPoints);

        // Data initialization

        dataInitialization( pdGridData, numOwnedPoints, horizon, uOwnedField, fNField, mPtr, dilatationTemporalField, bondDamagePtr, edbTemporalField);

        double *u1x = uOwnedField.get()+ 3;
	double *f1x = fNField.get()+ 3;

	/*
	 * INITIAL CONDITION
	 */
	*u1x = my_gamma;

        /*
	 * Pointers to data that don't change in this test
	 */

        double *xOverlapPtr = pdGridData.myX.get();
	double *mOwned = mPtr.get();
	double *bondDamage = bondDamagePtr.get();
	double *volumeOverlapPtr = pdGridData.cellVolume.get();
	int *localNeighborList = pdGridData.neighborhood.get();
	double *fInternalOverlapPtr = fNField.get();

        /*
	 * Current coordinates are fixed; Assign here
	 */

	TemporalField<double> yTemporalField(CURCOORD3D,numOwnedPoints);
	updateGeometry(xOverlapPtr,uOwnedField.get(),yTemporalField.getField(Field_ENUM::STEP_N).get(),numOwnedPoints*3);
	updateGeometry(xOverlapPtr,uOwnedField.get(),yTemporalField.getField(Field_ENUM::STEP_NP1).get(),numOwnedPoints*3);
	double *yN_OverlapPtr = yTemporalField.getField(Field_ENUM::STEP_N).get();
	double *yNp1_OverlapPtr = yTemporalField.getField(Field_ENUM::STEP_NP1).get();
        
	double* nullPtr = NULL;

        /*
	 * Compute initial force with elastic constitutive model
	 */
	MATERIAL_EVALUATION::computeInternalForceLinearElastic
	(
			xOverlapPtr,
			yN_OverlapPtr,
			mOwned,
			volumeOverlapPtr,
			dilatationTemporalField.getField(Field_ENUM::STEP_N).get(),
			bondDamage,
			fInternalOverlapPtr,
			nullPtr,
			localNeighborList,
			numOwnedPoints,
			K,
			MU,
            horizon
	);

        /*
	 * this is the correct analytical value for the initial condition
	 * ??? Why can't we get this closer than 1.0e-6 ???
	 */


	double f0 = 2.0 * 15.0 * E * my_gamma / 4.0 / (1+nu) / std::sqrt(2.0);
	double tolerance=1.0e-6;
	double rel_diff = std::fabs(f0-MAGNITUDE(f1x))/f0;

	
        TEST_COMPARE( rel_diff, <=, tolerance);


	double f = runPureShear("twoPoint_SLS_Relaxation.dat", K, MU, lambda_i, tau_b, numOwnedPoints, numSteps_stage_1, t_start, dt, yTemporalField, fNField, dilatationTemporalField,  edbTemporalField, u1x, f1x,xOverlapPtr, yN_OverlapPtr, yNp1_OverlapPtr, mOwned, volumeOverlapPtr, bondDamage, localNeighborList, fInternalOverlapPtr);

	
	/*
	 * Last value computed: tests time integrator against exact value
	 */


	/*
	 * step input deviatoric extension state
	 */

	double ed0=1.0e-6/sqrt(2);
	double m = 2.0;
	double alpha = 15.0 * E / (1+nu) / 2.0 / m;
	double tEnd=6.0;
	double fEnd = 2.0 * ed0*((1-lambda_i) * alpha +lambda_i * alpha * std::exp(-tEnd/tau_b));
	rel_diff = std::fabs(fEnd-f)/fEnd;
	
        TEST_COMPARE( rel_diff, <=, tolerance);

}

TEUCHOS_UNIT_TEST(TwoPoint_SLS_Relaxation,Case_3) {

	/*
	 * This produces the 'elastic' response
	 */

	double scale=0.99;
	Teuchos::ParameterList paramList = getParamList(scale);
	PeridigmNS::ViscoelasticMaterial mat(paramList);

        QUICKGRID::QuickGridData pdGridData = getTwoPointGridData();
	int numOwnedPoints = pdGridData.numPoints;
	TEST_ASSERT(2 == numOwnedPoints);
	TEST_ASSERT(4 == pdGridData.sizeNeighborhoodList);

        /*
	 * Material Props
	 */

        double K, MU, horizon, lambda_i, tau_b;

        setMaterialProperties( paramList, K, MU, horizon, lambda_i, tau_b);

        /*
	 * Time stepping data
	 */
	double t_start=0.0;
	double t_end = 3.0 * tau_b;
	size_t numSteps_stage_1(100);

	double dt = (t_end - t_start) / numSteps_stage_1;

        /*
	 * Displacement and Internal Force Vectors
	 */

	FieldSpec uSpec = DISPL3D;
	Field<double> uOwnedField(uSpec,pdGridData.numPoints);
	FieldSpec fNSpec = FORCE_DENSITY3D;
	Field_NS::Field<double> fNField(fNSpec,pdGridData.numPoints);

        /*
	 * Weighted Volume
	 */
	UTILITIES::Array<double> mPtr(numOwnedPoints);
	
	/*
	 * Dilatation: intialize to zero
	 */
	TemporalField<double> dilatationTemporalField(DILATATION,numOwnedPoints);
	
	/*
	 * Bond State and deviatoric back extension state
	 */
	UTILITIES::Array<double> bondDamagePtr(pdGridData.sizeNeighborhoodList-numOwnedPoints);
	
	TemporalField<double> edbTemporalField(DEVIATORIC_BACK_EXTENSION,pdGridData.sizeNeighborhoodList-numOwnedPoints);

        // Data initialization

        dataInitialization( pdGridData, numOwnedPoints, horizon, uOwnedField, fNField, mPtr, dilatationTemporalField, bondDamagePtr, edbTemporalField);

        double *u1x = uOwnedField.get()+ 3;
	double *f1x = fNField.get()+ 3;

	/*
	 * INITIAL CONDITION
	 */
	*u1x = my_gamma;

        /*
	 * Pointers to data that don't change in this test
	 */

        double *xOverlapPtr = pdGridData.myX.get();
	double *mOwned = mPtr.get();
	double *bondDamage = bondDamagePtr.get();
	double *volumeOverlapPtr = pdGridData.cellVolume.get();
	int *localNeighborList = pdGridData.neighborhood.get();
	double *fInternalOverlapPtr = fNField.get();

        /*
	 * Current coordinates are fixed; Assign here
	 */

	TemporalField<double> yTemporalField(CURCOORD3D,numOwnedPoints);
	updateGeometry(xOverlapPtr,uOwnedField.get(),yTemporalField.getField(Field_ENUM::STEP_N).get(),numOwnedPoints*3);
	updateGeometry(xOverlapPtr,uOwnedField.get(),yTemporalField.getField(Field_ENUM::STEP_NP1).get(),numOwnedPoints*3);
	double *yN_OverlapPtr = yTemporalField.getField(Field_ENUM::STEP_N).get();
	double *yNp1_OverlapPtr = yTemporalField.getField(Field_ENUM::STEP_NP1).get();
        
	double* nullPtr = NULL;

        /*
	 * Compute initial force with elastic constitutive model
	 */
	MATERIAL_EVALUATION::computeInternalForceLinearElastic
	(
			xOverlapPtr,
			yN_OverlapPtr,
			mOwned,
			volumeOverlapPtr,
			dilatationTemporalField.getField(Field_ENUM::STEP_N).get(),
			bondDamage,
			fInternalOverlapPtr,
			nullPtr,
			localNeighborList,
			numOwnedPoints,
			K,
			MU,
            horizon
	);

        /*
	 * this is the correct analytical value for the initial condition
	 * ??? Why can't we get this closer than 1.0e-6 ???
	 */
 

	double f0 = 2.0 * 15.0 * E * my_gamma / 4.0 / (1+nu) / std::sqrt(2.0);
	double tolerance=1.0e-6;
	double rel_diff = std::fabs(f0-MAGNITUDE(f1x))/f0;
	
        TEST_COMPARE( rel_diff, <=, tolerance);


	double f = runPureShear("twoPoint_Maxwell_Relaxation.dat", K, MU, lambda_i, tau_b, numOwnedPoints, numSteps_stage_1, t_start, dt, yTemporalField, fNField, dilatationTemporalField,  edbTemporalField, u1x, f1x,xOverlapPtr, yN_OverlapPtr, yNp1_OverlapPtr, mOwned, volumeOverlapPtr, bondDamage, localNeighborList, fInternalOverlapPtr);

	
	/*
	 * Last value computed: tests time integrator against exact value
	 */


	/*
	 * step input deviatoric extension state
	 */

	double ed0=1.0e-6/sqrt(2);
	double m = 2.0;
	double alpha = 15.0 * E / (1+nu) / 2.0 / m;
	double tEnd=6.0;
	double fEnd = 2.0 * ed0*((1-lambda_i) * alpha +lambda_i * alpha * std::exp(-tEnd/tau_b));
	rel_diff = std::fabs(fEnd-f)/fEnd;
	
        TEST_COMPARE( rel_diff, <=, tolerance);

}



int main
(int argc, char* argv[])
{
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
