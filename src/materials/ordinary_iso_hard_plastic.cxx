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
    /*Sorry about the run-on solution, copied from Mathematica*/
    return (4*alpha*(2*HARD_MODULUS*(-1 + alpha*lambdaN) + alpha*pointWiseYieldValue) - 
     (4*pow(2,2.0/3.0)*pow(alpha,2)*
        pow(HARD_MODULUS*(-1 + alpha*lambdaN) - alpha*pointWiseYieldValue,2))/
      pow(-4*pow(alpha,3)*pow(HARD_MODULUS,3) + 
        4*pow(alpha,6)*pow(HARD_MODULUS*lambdaN - pointWiseYieldValue,3) - 
        12*pow(alpha,5)*HARD_MODULUS*
         pow(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue,2) + 
        3*pow(alpha,4)*pow(HARD_MODULUS,2)*
         (4*HARD_MODULUS*lambdaN - 4*pointWiseYieldValue + 9*pow(tdNorm,2)) + 
        3*sqrt(3)*sqrt(pow(alpha,7)*pow(HARD_MODULUS,2)*pow(tdNorm,2)*
           (8*pow(HARD_MODULUS,3)*pow(-1 + alpha*lambdaN,3) + 
             24*pow(alpha,2)*HARD_MODULUS*(-1 + alpha*lambdaN)*
              pow(pointWiseYieldValue,2) - 
             8*pow(alpha,3)*pow(pointWiseYieldValue,3) - 
             3*alpha*pow(HARD_MODULUS,2)*
              (8*pow(-1 + alpha*lambdaN,2)*pointWiseYieldValue - 9*pow(tdNorm,2))))
        ,1.0/3.0) - 2*pow(2,1.0/3.0)*
      pow(-4*pow(alpha,3)*pow(HARD_MODULUS,3) + 
        4*pow(alpha,6)*pow(HARD_MODULUS*lambdaN - pointWiseYieldValue,3) - 
        12*pow(alpha,5)*HARD_MODULUS*
         pow(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue,2) + 
        3*pow(alpha,4)*pow(HARD_MODULUS,2)*
         (4*HARD_MODULUS*lambdaN - 4*pointWiseYieldValue + 9*pow(tdNorm,2)) + 
        3*sqrt(3)*sqrt(pow(alpha,7)*pow(HARD_MODULUS,2)*pow(tdNorm,2)*
           (8*pow(HARD_MODULUS,3)*pow(-1 + alpha*lambdaN,3) + 
             24*pow(alpha,2)*HARD_MODULUS*(-1 + alpha*lambdaN)*
              pow(pointWiseYieldValue,2) - 
             8*pow(alpha,3)*pow(pointWiseYieldValue,3) - 
             3*alpha*pow(HARD_MODULUS,2)*
              (8*pow(-1 + alpha*lambdaN,2)*pointWiseYieldValue - 9*pow(tdNorm,2))))
        ,1.0/3.0))/(12.*pow(alpha,2)*HARD_MODULUS);

}

//template<typename ScalarT>
//void updateLambdaNP1AD
//(
		//double tdNorm,
		//const double lambdaN,
        //double pointWiseYieldValue,
		//double alpha,
        //double HARD_MODULUS,
        //double dt
//)
//{
    //[>Sorry about the run-on solution, copied from Mathematica<]
    //return (4*alpha*(-2*dt*HARD_MODULUS + alpha*(2*HARD_MODULUS*lambdaN + pointWiseYieldValue)) - (4*pow(2,2.0/3.0)*pow(alpha,2)*pow(dt*HARD_MODULUS + alpha*(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue),2)) /pow(-4*pow(alpha,3)*pow(dt,3)*pow(HARD_MODULUS,3) + 4*pow(alpha,6)*pow(HARD_MODULUS*lambdaN - pointWiseYieldValue,3) - 12*pow(alpha,5)*dt*HARD_MODULUS* pow(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue,2) + 3*pow(alpha,4)*pow(dt,2)*pow(HARD_MODULUS,2)* (4*HARD_MODULUS*lambdaN - 4*pointWiseYieldValue + 9*pow(tdNorm,2)) + 3*sqrt(3)*sqrt(pow(alpha,7)*pow(dt,2)*pow(HARD_MODULUS,2)*pow(tdNorm,2)* (-8*pow(dt,3)*pow(HARD_MODULUS,3) + 8*pow(alpha,3)*pow(HARD_MODULUS*lambdaN - pointWiseYieldValue,3) - 24*pow(alpha,2)*dt*HARD_MODULUS* pow(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue,2) + 3*alpha*pow(dt,2)*pow(HARD_MODULUS,2)* (8*HARD_MODULUS*lambdaN - 8*pointWiseYieldValue + 9*pow(tdNorm,2)))), 1.0/3.0) - 2*pow(2,1.0/3.0)* pow(-4*pow(alpha,3)*pow(dt,3)*pow(HARD_MODULUS,3) + 4*pow(alpha,6)*pow(HARD_MODULUS*lambdaN - pointWiseYieldValue,3) - 12*pow(alpha,5)*dt*HARD_MODULUS* pow(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue,2) + 3*pow(alpha,4)*pow(dt,2)*pow(HARD_MODULUS,2)* (4*HARD_MODULUS*lambdaN - 4*pointWiseYieldValue + 9*pow(tdNorm,2)) + 3*sqrt(3)*sqrt(pow(alpha,7)*pow(dt,2)*pow(HARD_MODULUS,2)*pow(tdNorm,2)* (-8*pow(dt,3)*pow(HARD_MODULUS,3) + 8*pow(alpha,3)*pow(HARD_MODULUS*lambdaN - pointWiseYieldValue,3) - 24*pow(alpha,2)*dt*HARD_MODULUS* pow(-(HARD_MODULUS*lambdaN) + pointWiseYieldValue,2) + 3*alpha*pow(dt,2)*pow(HARD_MODULUS,2)* (8*HARD_MODULUS*lambdaN - 8*pointWiseYieldValue + 9*pow(tdNorm,2)))), 1.0/3.0))/(12.*pow(alpha,2)*HARD_MODULUS);

//};

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
        *lambdaNP1 =  updateLambdaNP1(tdNorm, *lambdaN, pointWiseYieldValue, alpha, H);
		
        /*
		 * Evaluate yield function
		 */
		double f = tdNorm * tdNorm / 2 - pointWiseYieldValue - HARD_MODULUS*(*lambdaNP1);
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
			//*lambdaNP1 = *lambdaN + deltaLambda;
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
};

//template<typename ScalarT>
//void computeInternalForceIsotropicHardeningPlasticAD
//(
		//const double* xOverlap,
		//const ScalarT* yNP1Overlap,
		//const double* mOwned,
		//const double* volumeOverlap,
		//const ScalarT* dilatationOwned,
		//const double* bondDamage,
		//const double* dsfOwned,
		//const double* deviatoricPlasticExtensionStateN,
		//ScalarT* deviatoricPlasticExtensionStateNp1,
		//const double* lambdaN,
		//ScalarT* lambdaNP1,
		//ScalarT* fInternalOverlap,
		//const int*  localNeighborList,
		//int numOwnedPoints,
		//double BULK_MODULUS,
		//double SHEAR_MODULUS,
		//double HORIZON,
		//double yieldStress,
		//double HARD_MODULUS,
		//double dt
//)
//{

	/*
	 * Compute processor local contribution to internal force
	 */
	//double K = BULK_MODULUS;
	//double MU = SHEAR_MODULUS;
	//double OMEGA=1.0;
	//double DELTA=HORIZON;
	//double H=HARD_MODULUS;
	//double deltaT=dt;
	/*
	 * 3d variety of yield value
	 */
	//double yieldValue = 75.0 * yieldStress * yieldStress / 8 / M_PI / pow(DELTA,5);
	/*
	 * Planar variety of yield value
	 */
////		double THICKNESS=1.0;
////		double yieldValue = 225.0 * yieldStress * yieldStress / 8 / M_PI / THICKNESS / pow(DELTA,4);
////		double yieldValue = 0.5 * pow(15*yieldStress/weightedVol,2) * M_PI * THICKNESS * pow(DELTA,4) / 16.0;



	//const double *xOwned = xOverlap;
	//const ScalarT *yOwned = yNP1Overlap;
	//const double *m = mOwned;
	//const double *v = volumeOverlap;
	//const ScalarT *theta = dilatationOwned;
	//ScalarT *fOwned = fInternalOverlap;

	//const int *neighPtr = localNeighborList;
	//double cellVolume, alpha, dx_X, dy_X, dz_X, zeta, edpN;
    //ScalarT dx_Y, dy_Y, dz_Y, dY, ed, tdTrial, t, ti, td;
	//for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, m++, theta++, lambdaN++, lambdaNP1++, dsfOwned++){

		//int numNeigh = *neighPtr; neighPtr++;
		//const double *X = xOwned;
		//const ScalarT *Y = yOwned;
		//double weightedVol = *m;
		//alpha = *dsfOwned * 15.0*MU/weightedVol;
		//double selfCellVolume = v[p];
		//ScalarT c = 3 * K * (*theta) * OMEGA / weightedVol;
		//ScalarT deltaLambda=0.0;

		/*
		 * Compute norm of trial stress
		 */
		//ScalarT tdNorm = 0.0;
		//tdNorm = computeDeviatoricForceStateNormAD(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,X,Y,xOverlap,yNP1Overlap,v,alpha,OMEGA);

		//double pointWiseYieldValue = *dsfOwned * (*dsfOwned) * yieldValue;
        /*
         * Compute lambdaNP1 using a backward Euler implicit scheme
        */
        //*lambdaNP1 =  updateLambdaNP1(tdNorm, *lambdaN, pointWiseYieldValue, alpha, H, deltaT);
		
        /*
		 * Evaluate yield function
		 */
		//SaclarT f = tdNorm * tdNorm / 2 - pointWiseYieldValue - HARD_MODULUS*(*lambdaNP1);
		//bool elastic = true;

////		std::cout << "Point id = " << p << std::endl;
////		std::cout << "\tyieldStress/m^(4/5) = " << yieldStress/pow(weightedVol,4/5) << std::endl;
////		std::cout << "\tYield Value = " << yieldValue << "; tdNorm * tdNorm / 2 = " << tdNorm * tdNorm / 2 << std::endl;
		//if(f>0){
			/*
			 * This step is incrementally plastic
			 */
			////			std::cout << "\t PLASTIC" << std::endl;
			//elastic = false;
            //deltaLambda= (*lambdaNP1 - *lambdaN) / dt;
			////deltaLambda=( tdNorm / sqrt(2.0*pointWiseYieldValue) - 1.0 ) / alpha;
			///[>lambdaNP1 = *lambdaN + deltaLambda;
		//} else {
//// 			std::cout << "\t ELASTIC" << std::endl;
			//*lambdaNP1 = *lambdaN;
		//}

		//for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++, deviatoricPlasticExtensionStateN++, deviatoricPlasticExtensionStateNp1++){
			//int localId = *neighPtr;
			//cellVolume = v[localId];
			//const double *XP = &xOverlap[3*localId];
			//const ScalarT *YP = &yNP1Overlap[3*localId];
			//dx_X = XP[0]-X[0];
			//dy_X = XP[1]-X[1];
			//dz_X = XP[2]-X[2];
			//zeta = sqrt(dx_X*dx_X+dy_X*dy_X+dz_X*dz_X);
			//dx_Y = YP[0]-Y[0];
			//dy_Y = YP[1]-Y[1];
			//dz_Y = YP[2]-Y[2];
			//dY = sqrt(dx_Y*dx_Y+dy_Y*dy_Y+dz_Y*dz_Y);
			/*
			 * Deviatoric extension state
			 */
			//ed = dY-zeta-*theta*zeta/3;

			/*
			 * Deviatoric plastic extension state from last step
			 */
			//edpN = *deviatoricPlasticExtensionStateN;

			/*
			 * Compute trial stress
			 */
			//tdTrial = alpha * OMEGA * (ed - edpN);

			/*
			 * Evaluate yield function
			 */
			//if(elastic){
				/*
				 * Elastic case
				 */
				//td = tdTrial;

				/*
				 * Therefore edpNp1 = edpN
				 */
				//*deviatoricPlasticExtensionStateNp1 = *deviatoricPlasticExtensionStateN;

			//} else {
				/*
				 * Compute deviatoric force state
				 */
				////td = sqrt(2.0*pointWiseYieldValue) * tdTrial / tdNorm;
				//td = tdTrial / (1+alpha*deltaLambda);

				/*
				 * Update deviatoric plastic deformation state
				 */
				//*deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;

////				std::cout << "Neighbor Id = " << localId << "; Updating deviatoricPlasticExtensionState = " << *deviatoricPlasticExtensionState << std::endl;
			//}
////			std::cout << "\tNeighbor Id = " << localId << "\n\ttd = " << td;
			/*
			 * Compute isotropic part of force state
			 */
			//ti = c * zeta;

			/*
			 * Force state (with damage)
			 */
			//double d=(1.0-*bondDamage);
			//t = d*(ti + d*td);

			/*
			 * Assemble pair wise force function
			 */
			//ScalarT fx = t * dx_Y / dY;
			//ScalarT fy = t * dy_Y / dY;
			//ScalarT fz = t * dz_Y / dY;

			//*(fOwned+0) += fx*cellVolume;
			//*(fOwned+1) += fy*cellVolume;
			//*(fOwned+2) += fz*cellVolume;
			//fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			//fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			//fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
		//}

	//}
//}



//template void updateLambdaNP1AD<double>
//(
		//double tdNorm,
		//const double lambdaN,
        //double pointWiseYieldValue,
		//double alpha,
        //double HARD_MODULUS,
        //double dt
//);


//template void updateLambdaNP1AD<Sacado::Fad::DFad<double> >
//(
		//double tdNorm,
		//const double lambdaN,
        //double pointWiseYieldValue,
		//double alpha,
        //double HARD_MODULUS,
        //double dt
//);

//[>* Explicit template instantiation for double. <]
//template void computeInternalForceIsotropicHardeningPlasticAD<double>
//(
		//const double* xOverlap,
		//const double* yNP1Overlap,
		//const double* mOwned,
		//const double* volumeOverlap,
		//const double* dilatationOwned,
		//const double* bondDamage,
		//const double* dsfOwned,
		//const double* deviatoricPlasticExtensionStateN,
		//double* deviatoricPlasticExtensionStateNp1,
		//const double* lambdaN,
		//double* lambdaNP1,
		//double* fInternalOverlap,
		//const int*  localNeighborList,
		//int numOwnedPoints,
		//double BULK_MODULUS,
		//double SHEAR_MODULUS,
		//double HORIZON,
		//double yieldStress,
		//double HARD_MODULUS,
		//double dt
//);


//[>* Explicit template instantiation for Sacado::Fad::DFad<double>. <]
//template void computeInternalForceIsotropicHardeningPlasticAD<Sacado::Fad::DFad<double> >
//(
		//const double* xOverlap,
		//const Sacado::Fad::DFad<double>* yNP1Overlap,
		//const double* mOwned,
		//const double* volumeOverlap,
		//const Sacado::Fad::DFad<double>* dilatationOwned,
		//const double* bondDamage,
		//const double* dsfOwned,
		//const double* deviatoricPlasticExtensionStateN,
		//Sacado::Fad::DFad<double>* deviatoricPlasticExtensionStateNp1,
		//const double* lambdaN,
		//Sacado::Fad::DFad<double>* lambdaNP1,
		//Sacado::Fad::DFad<double>* fInternalOverlap,
		//const int*  localNeighborList,
		//int numOwnedPoints,
		//double BULK_MODULUS,
		//double SHEAR_MODULUS,
		//double HORIZON,
		//double yieldStress,
		//double HARD_MODULUS,
		//double dt
//);

}

