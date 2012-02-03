/*! \file ordinary_iso_hard_plastic.cxx */

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
#include <cmath>
#include <Sacado.hpp>
#include "ordinary_elastic_plastic.h"
#include "ordinary_iso_hard_plastic.h"

namespace MATERIAL_EVALUATION {

/*
 * Simple sign() function
 */
int sign(double v) {
    return v > 0.0 ? 1 : (v < 0.0 ? -1 : 0);
}

/*
 * function for finding the RHS of the root solve equation
 */
double lambdaRHS
(
		double tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS,
        double lambdaNP1
)
{
    //return tdNorm*tdNorm/2.0/(1.0+alpha*(lambdaNP1 - lambdaN))/(1.0+alpha*(lambdaNP1 - lambdaN))
                    //- pointWiseYieldValue;
    return tdNorm*tdNorm/2.0/(1.0+alpha*(lambdaNP1 - lambdaN))/(1.0+alpha*(lambdaNP1 - lambdaN))
                    - HARD_MODULUS*lambdaNP1 - pointWiseYieldValue;
}

/*
 * Finds the value of lambda at (n+1) using the solution of a backward-Euler
 * finite difference discretization of delta lambda
 */
double updateLambdaNP1
(
		double tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS
)
{
    // We'll start with a bisection root finder.
    double a = 0.0;
    double b = 1.0e5;

    double TOL1 = 1.0e-35;

    double FA = lambdaRHS ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, a);
    double p = 0.0;
                
    //Converge to the first tolerance using bisection method
    for( int i=0; i < 1000000; i++){
        
        // Calculate the midpoint of the interval.
        p = a + (b - a)/2;
        // The function evaluated at the midpoint.
        double FP = lambdaRHS ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, p);
        
        //Stop conditions
        if( FP == 0 || (b-a)/2.0 < TOL1) 
            break;

        //Checks to see if the function changes sign over the interval [a,p]
        if( sign(FA)*sign(FP) > 0){
           a = p;
           FA = FP;
        } else {
           b = p; 
        }

        if( i == 999999) 
            std::cout << "failed to converge" << std::endl;
    };

    //Now we switch to the secant method to polish off the root.
    //double TOL2 = 1.0e-6;
    //double p=0.0;
    //double p0=a;
    //double p1=b;
    //double q0=lambdaRHS ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, p0);
    //double q1=lambdaRHS ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, p1);

    //for( int i=0; i < 10000; i++){
    
        ////Find the root of the secant line between points p1 and p2
        //p = p1 - q1*(p1-p0)/(q1-q0);
        
        ////Stop condition
        //if (fabs(p) < 1e-10 || fabs(p-p1) < TOL2) break;
        
        //p0=p1;
        //q0=q1;
        //p1=p;
        //q1=lambdaRHS ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, p);

        //if(i == 9999){
            //p = 0.0;
        //}
    //};
    
    return p;
}
    

void computeInternalForceIsotropicHardeningPlastic
(
		const double* xOverlap,
		const double* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		const double* deviatoricPlasticExtensionStateN,
		double* deviatoricPlasticExtensionStateNp1,
		const double* lambdaN,
		double* lambdaNP1,
		double* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
		double HORIZON,
		double yieldStress,
		double HARD_MODULUS
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;
	double OMEGA=1.0;
	double DELTA=HORIZON;
	double H=HARD_MODULUS;
	/*
	 * 3d variety of yield value
	 */
	double yieldValue = 75.0 * yieldStress * yieldStress / 8 / M_PI / pow(DELTA,5);
	/*
	 * Planar variety of yield value
	 */
//		double THICKNESS=1.0;
//		double yieldValue = 225.0 * yieldStress * yieldStress / 8 / M_PI / THICKNESS / pow(DELTA,4);
//		double yieldValue = 0.5 * pow(15*yieldStress/weightedVol,2) * M_PI * THICKNESS * pow(DELTA,4) / 16.0;



	const double *xOwned = xOverlap;
	const double *yOwned = yNP1Overlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const double *theta = dilatationOwned;
	double *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, dx, dy, dz, zeta, dY, t, ti, td, ed, edpN, tdTrial;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, m++, theta++, lambdaN++, lambdaNP1++, dsfOwned++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *Y = yOwned;
		double weightedVol = *m;
		alpha = *dsfOwned * 15.0*MU/weightedVol;
		double selfCellVolume = v[p];
		double c = 3 * K * (*theta) * OMEGA / weightedVol;
		double deltaLambda=0.0;

		/*
		 * Compute norm of trial stress
		 */
		double tdNorm = 0.0;
		tdNorm = computeDeviatoricForceStateNorm(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,X,Y,xOverlap,yNP1Overlap,v,alpha,OMEGA);

		double pointWiseYieldValue = *dsfOwned * (*dsfOwned) * yieldValue;
        /*
         * Compute lambdaNP1 using a backward Euler implicit scheme
        */
        if (tdNorm * tdNorm / 2 - pointWiseYieldValue > 0){
            *lambdaNP1 =  updateLambdaNP1(tdNorm, *lambdaN, pointWiseYieldValue, alpha, H);
        } else {
            *lambdaNP1 = 0.0;
        }
        /*
		 * Evaluate yield function
		 */
        double f = tdNorm * tdNorm / 2 - pointWiseYieldValue - HARD_MODULUS*(*lambdaNP1);
        //double f = tdNorm * tdNorm / 2 - pointWiseYieldValue;
		bool elastic = true;

//		std::cout << "Point id = " << p << std::endl;
//		std::cout << "\tyieldStress/m^(4/5) = " << yieldStress/pow(weightedVol,4/5) << std::endl;
//		std::cout << "\tYield Value = " << yieldValue << "; tdNorm * tdNorm / 2 = " << tdNorm * tdNorm / 2 << std::endl;
		if(f>0){
			/*
			 * This step is incrementally plastic
			 */
			//			std::cout << "\t PLASTIC" << std::endl;
			elastic = false;
            deltaLambda= (*lambdaNP1 - *lambdaN);
            //if(deltaLambda > 0){
                //std::cout << "DeltaLambda is " << deltaLambda << std::endl;
            //};
            //deltaLambda=( tdNorm / sqrt(2.0*pointWiseYieldValue) - 1.0 ) / alpha;
            *lambdaNP1 = *lambdaN + deltaLambda;
		} else {
//			std::cout << "\t ELASTIC" << std::endl;
			*lambdaNP1 = *lambdaN;
		}

		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++, deviatoricPlasticExtensionStateN++, deviatoricPlasticExtensionStateNp1++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const double *YP = &yNP1Overlap[3*localId];
			dx = XP[0]-X[0];
			dy = XP[1]-X[1];
			dz = XP[2]-X[2];
			zeta = sqrt(dx*dx+dy*dy+dz*dz);
			dx = YP[0]-Y[0];
			dy = YP[1]-Y[1];
			dz = YP[2]-Y[2];
			dY = sqrt(dx*dx+dy*dy+dz*dz);
			/*
			 * Deviatoric extension state
			 */
			ed = dY-zeta-*theta*zeta/3;

			/*
			 * Deviatoric plastic extension state from last step
			 */
			edpN = *deviatoricPlasticExtensionStateN;

			/*
			 * Compute trial stress
			 */
			tdTrial = alpha * OMEGA * (ed - edpN);

			/*
			 * Evaluate yield function
			 */
			if(elastic){
				/*
				 * Elastic case
				 */
				td = tdTrial;

				/*
				 * Therefore edpNp1 = edpN
				 */
				*deviatoricPlasticExtensionStateNp1 = *deviatoricPlasticExtensionStateN;

			} else {
				/*
				 * Compute deviatoric force state
				 */
                //td = sqrt(2.0*pointWiseYieldValue) * tdTrial / tdNorm;
                td = tdTrial / (1+alpha*deltaLambda);

				/*
				 * Update deviatoric plastic deformation state
				 */
				*deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;

//				std::cout << "Neighbor Id = " << localId << "; Updating deviatoricPlasticExtensionState = " << *deviatoricPlasticExtensionState << std::endl;
			}
//			std::cout << "\tNeighbor Id = " << localId << "\n\ttd = " << td;
			/*
			 * Compute isotropic part of force state
			 */
			ti = c * zeta;

			/*
			 * Force state (with damage)
			 */
			double d=(1.0-*bondDamage);
			t = d*(ti + d*td);

			/*
			 * Assemble pair wise force function
			 */
			double fx = t * dx / dY;
			double fy = t * dy / dY;
			double fz = t * dz / dY;

			*(fOwned+0) += fx*cellVolume;
			*(fOwned+1) += fy*cellVolume;
			*(fOwned+2) += fz*cellVolume;
			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
		}

	}
}

template<typename ScalarT>
void computeInternalForceIsotropicHardeningPlasticAD
(
		const double* xOverlap,
		const ScalarT* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const ScalarT* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		const double* deviatoricPlasticExtensionStateN,
		ScalarT* deviatoricPlasticExtensionStateNp1,
		const double* lambdaN,
		ScalarT* lambdaNP1,
		ScalarT* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
		double HORIZON,
		double yieldStress,
		double HARD_MODULUS
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;
	double OMEGA=1.0;
	double DELTA=HORIZON;
	double H=HARD_MODULUS;
	/*
	 * 3d variety of yield value
	 */
	double yieldValue = 75.0 * yieldStress * yieldStress / 8 / M_PI / pow(DELTA,5);
	/*
	 * Planar variety of yield value
	 */
//		double THICKNESS=1.0;
//		double yieldValue = 225.0 * yieldStress * yieldStress / 8 / M_PI / THICKNESS / pow(DELTA,4);
//		double yieldValue = 0.5 * pow(15*yieldStress/weightedVol,2) * M_PI * THICKNESS * pow(DELTA,4) / 16.0;


	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yNP1Overlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, dx_X, dy_X, dz_X, zeta, edpN;
    ScalarT dx_Y, dy_Y, dz_Y, dY, ed, tdTrial, t, ti, td;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, m++, theta++, lambdaN++, lambdaNP1++, dsfOwned++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		double weightedVol = *m;
		alpha = *dsfOwned * 15.0*MU/weightedVol;
		double selfCellVolume = v[p];
		ScalarT c = 3 * K * (*theta) * OMEGA / weightedVol;
		ScalarT deltaLambda=0.0;

		/*
		 * Compute norm of trial stress
		 */
		ScalarT tdNorm = 0.0;
		tdNorm = computeDeviatoricForceStateNormAD(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,X,Y,xOverlap,yNP1Overlap,v,alpha,OMEGA);

		/*
		 * Evaluate yield function
		 */
		double pointWiseYieldValue = *dsfOwned * (*dsfOwned) * yieldValue;

        /*
         * Compute lambdaNP1 using a backward Euler implicit scheme
        */
        if (tdNorm * tdNorm / 2 - pointWiseYieldValue > 0){
            *lambdaNP1 =  updateLambdaNP1AD(tdNorm, *lambdaN, pointWiseYieldValue, alpha, H);
        } else {
            *lambdaNP1 = 0.0;
        }


        ScalarT f = tdNorm * tdNorm / 2 - pointWiseYieldValue - HARD_MODULUS*(*lambdaNP1);
		//ScalarT f = tdNorm * tdNorm / 2 - pointWiseYieldValue;
		bool elastic = true;

//		std::cout << "Point id = " << p << std::endl;
//		std::cout << "\tyieldStress/m^(4/5) = " << yieldStress/pow(weightedVol,4/5) << std::endl;
//		std::cout << "\tYield Value = " << yieldValue << "; tdNorm * tdNorm / 2 = " << tdNorm * tdNorm / 2 << std::endl;
		if(f>0){
			/*
			 * This step is incrementally plastic
			 */
			//			std::cout << "\t PLASTIC" << std::endl;
			elastic = false;
            deltaLambda= (*lambdaNP1 - *lambdaN);
			//deltaLambda=( tdNorm / sqrt(2.0*pointWiseYieldValue) - 1.0 ) / alpha;
			*lambdaNP1 = *lambdaN + deltaLambda;
		} else {
// 			std::cout << "\t ELASTIC" << std::endl;
			*lambdaNP1 = *lambdaN;
		}

		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++, deviatoricPlasticExtensionStateN++, deviatoricPlasticExtensionStateNp1++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yNP1Overlap[3*localId];
			dx_X = XP[0]-X[0];
			dy_X = XP[1]-X[1];
			dz_X = XP[2]-X[2];
			zeta = sqrt(dx_X*dx_X+dy_X*dy_X+dz_X*dz_X);
			dx_Y = YP[0]-Y[0];
			dy_Y = YP[1]-Y[1];
			dz_Y = YP[2]-Y[2];
			dY = sqrt(dx_Y*dx_Y+dy_Y*dy_Y+dz_Y*dz_Y);
			/*
			 * Deviatoric extension state
			 */
			ed = dY-zeta-*theta*zeta/3;

			/*
			 * Deviatoric plastic extension state from last step
			 */
			edpN = *deviatoricPlasticExtensionStateN;

			/*
			 * Compute trial stress
			 */
			tdTrial = alpha * OMEGA * (ed - edpN);

			/*
			 * Evaluate yield function
			 */
			if(elastic){
				/*
				 * Elastic case
				 */
				td = tdTrial;

				/*
				 * Therefore edpNp1 = edpN
				 */
				*deviatoricPlasticExtensionStateNp1 = *deviatoricPlasticExtensionStateN;

			} else {
				/*
				 * Compute deviatoric force state
				 */
				//td = sqrt(2.0*pointWiseYieldValue) * tdTrial / tdNorm;
                td = tdTrial / (1+alpha*deltaLambda);

				/*
				 * Update deviatoric plastic deformation state
				 */
				*deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;

//				std::cout << "Neighbor Id = " << localId << "; Updating deviatoricPlasticExtensionState = " << *deviatoricPlasticExtensionState << std::endl;
			}
//			std::cout << "\tNeighbor Id = " << localId << "\n\ttd = " << td;
			/*
			 * Compute isotropic part of force state
			 */
			ti = c * zeta;

			/*
			 * Force state (with damage)
			 */
			double d=(1.0-*bondDamage);
			t = d*(ti + d*td);

			/*
			 * Assemble pair wise force function
			 */
			ScalarT fx = t * dx_Y / dY;
			ScalarT fy = t * dy_Y / dY;
			ScalarT fz = t * dz_Y / dY;

			*(fOwned+0) += fx*cellVolume;
			*(fOwned+1) += fy*cellVolume;
			*(fOwned+2) += fz*cellVolume;
			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
		}

	}
}

/*
 * Simple signAD() function
 */
template<typename ScalarT>
ScalarT signAD(ScalarT v) {
    return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0);
}

template<typename ScalarT>
ScalarT lambdaRHSAD
(
		ScalarT tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS,
        ScalarT lambdaNP1
)
{
    //return tdNorm*tdNorm/2.0/(1.0+alpha*(lambdaNP1 - lambdaN))/(1.0+alpha*(lambdaNP1 - lambdaN))
                    //- pointWiseYieldValue;
    return tdNorm*tdNorm/2.0/(1.0+alpha*(lambdaNP1 - lambdaN))/(1.0+alpha*(lambdaNP1 - lambdaN))
                    - HARD_MODULUS*lambdaNP1 - pointWiseYieldValue;
}

/*
 * Finds the value of lambda at (n+1) using the solution of a backward-Euler
 * finite difference discretization of delta lambda
 */
template<typename ScalarT>
ScalarT updateLambdaNP1AD
(
		ScalarT tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS
)
{
    // We'll start with a bisection root finder.
    ScalarT a = 0.0;
    ScalarT b = 1.0e5;

    double TOL1 = 1.0e-35;

    ScalarT FA = lambdaRHSAD ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, a);
    ScalarT p = 0.0;
                
    //Converge to the first tolerance using bisection method
    for( int i=0; i < 100000; i++){
        
        // Calculate the midpoint of the interval.
        p = a + (b - a)/2;
        // The function evaluated at the midpoint.
        ScalarT FP = lambdaRHSAD ( tdNorm, lambdaN, pointWiseYieldValue, alpha, HARD_MODULUS, p);
        
        //Stop conditions
        if( FP == 0 || (b-a)/2.0 < TOL1) 
            break;

        //Checks to see if the function changes sign over the interval [a,p]
        if( signAD(FA)*signAD(FP) > 0){
           a = p;
           FA = FP;
        } else {
           b = p; 
        }

        if( i == 99999) 
            std::cout << "failed to converge" << std::endl;
    };

}

/** Explicit template instantiation for double. */
template void computeInternalForceIsotropicHardeningPlasticAD<double>
(
		const double* xOverlap,
		const double* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		const double* deviatoricPlasticExtensionStateN,
		double* deviatoricPlasticExtensionStateNp1,
		const double* lambdaN,
		double* lambdaNP1,
		double* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
		double HORIZON,
		double yieldStress,
		double HARD_MODULUS
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceIsotropicHardeningPlasticAD<Sacado::Fad::DFad<double> >
(
		const double* xOverlap,
		const Sacado::Fad::DFad<double>* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const Sacado::Fad::DFad<double>* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		const double* deviatoricPlasticExtensionStateN,
		Sacado::Fad::DFad<double>* deviatoricPlasticExtensionStateNp1,
		const double* lambdaN,
		Sacado::Fad::DFad<double>* lambdaNP1,
		Sacado::Fad::DFad<double>* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
		double HORIZON,
		double yieldStress,
		double HARD_MODULUS
);

/** Explicit template instantiation for int. */
template double signAD<double> 
(
        double v
);

/** Explicit template instantiation for Sacado::Fad::DFad<int>. */
template Sacado::Fad::DFad<double> signAD<Sacado::Fad::DFad<double> >
(
        Sacado::Fad::DFad<double> v
);

/** Explicit template instantiation for double. */
template double updateLambdaNP1AD<double>
(
		double tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS
);


/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template Sacado::Fad::DFad<double> updateLambdaNP1AD<Sacado::Fad::DFad<double> >
(
		Sacado::Fad::DFad<double> tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS
);


/** Explicit template instantiation for double. */
template double lambdaRHSAD<double>
(
		double tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS,
        double lambdaNP1
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template Sacado::Fad::DFad<double> lambdaRHSAD<Sacado::Fad::DFad<double> >
(
		Sacado::Fad::DFad<double> tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS,
        Sacado::Fad::DFad<double> lambdaNP1
);

}
