/*! \file PdMaterialUtilities.cxx */

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

#include "PdMaterialUtilities.h"
#include "ordinary_elastic_plastic.h"
#include "ordinary_utilities.h"
#include <math.h>
#include <iostream>
#include <string>

#include <functional>
using std::binary_function;

//
//namespace PdMaterialUtilities {
//
//
//const std::string NAMESPACE="PdMaterialUtilities::";
//
//void updateGeometry
//(
//		const double* xOverlap,
//		const double* uOverlap,
//		const double* velocityOverlap,
//		double* yOverlap,
//		int overLapLength,
//		double dt
//)
//{
//	const double* x = xOverlap;
//	const double* u = uOverlap;
//	const double* v = velocityOverlap;
//	double*       y = yOverlap;
//
//	int length = overLapLength;
//	for(;x!=xOverlap+length;x++,u++,v++,y++)
//		*y = *x + *u + *v * dt;
//}
//
//
////! Version of computeDilataion for arbitrary list of ownedIDs
//void computeDilatation
//(
//		const double* xOverlap,
//		const double* yOverlap,
//		const double *mOwned,
//		const double* volumeOverlap,
//		const double* bondDamage,
//		double* dilatationOwned,
//        const int* ownedIDs,
//		const int* localNeighborList,
//		int numOwnedPoints
//)
//{
//	double OMEGA=1.0;
//	const double *xOwned = xOverlap;
//	const double *yOwned = yOverlap;
//	const double *m = mOwned;
//	const double *v = volumeOverlap;
//	double *theta = dilatationOwned;
//	double cellVolume;
//	const int *neighPtr = localNeighborList;
//	for(int p=0; p<numOwnedPoints;p++){
//
//        int ID = ownedIDs[p];
//        xOwned = &xOverlap[3*ID];
//        yOwned = &yOverlap[3*ID];
//        m = &mOwned[ID];
//        theta = &dilatationOwned[ID];
//
//		int numNeigh = *neighPtr; neighPtr++;
//		const double *X = xOwned;
//		const double *Y = yOwned;
//		*theta = 0;
//		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
//			int localId = *neighPtr;
//			cellVolume = v[localId];
//			const double *XP = &xOverlap[3*localId];
//			const double *YP = &yOverlap[3*localId];
//			double dx = XP[0]-X[0];
//			double dy = XP[1]-X[1];
//			double dz = XP[2]-X[2];
//			double zetaSqared = dx*dx+dy*dy+dz*dz;
//			dx = YP[0]-Y[0];
//			dy = YP[1]-Y[1];
//			dz = YP[2]-Y[2];
//			double dY = dx*dx+dy*dy+dz*dz;
//			double d = sqrt(zetaSqared);
//			double e = sqrt(dY)-d;
//			*theta += 3.0*OMEGA*(1.0-*bondDamage)*d*e*cellVolume/(*m);
//		}
//
//	}
//}
//
//
//void computeInternalForceLinearElastic
//(
//		const double* xOverlap,
//		const double* yOverlap,
//		const double* mOwned,
//		const double* volumeOverlap,
//		const double* dilatationOwned,
//		const double* bondDamage,
//		double* fInternalOverlap,
//		const int*  ownedIDs,
//		const int*  localNeighborList,
//		int numOwnedPoints,
//		double BULK_MODULUS,
//		double SHEAR_MODULUS
//)
//{
//
//	/*
//	 * Compute processor local contribution to internal force
//	 */
//	double K = BULK_MODULUS;
//	double MU = SHEAR_MODULUS;
//	double OMEGA=1.0;
//
//	const double *xOwned = xOverlap;
//	const double *yOwned = yOverlap;
//	const double *m = mOwned;
//	const double *v = volumeOverlap;
//	const double *theta = dilatationOwned;
//	double *fOwned = fInternalOverlap;
//
//	const int *neighPtr = localNeighborList;
//	double cellVolume, alpha, dx, dy, dz, zeta, dY, t;
//	for(int p=0;p<numOwnedPoints;p++){
//
//        int ID = ownedIDs[p];
//        xOwned = &xOverlap[3*ID];
//        yOwned = &yOverlap[3*ID];
//        fOwned = &fInternalOverlap[3*ID];
//        m = &mOwned[ID];
//        theta = &dilatationOwned[ID];
//
//		int numNeigh = *neighPtr; neighPtr++;
//		const double *X = xOwned;
//		const double *Y = yOwned;
//		alpha = 15.0*MU/(*m);
//		double selfCellVolume = v[p];
//		double c1 = OMEGA*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
//		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
//			int localId = *neighPtr;
//			cellVolume = v[localId];
//			const double *XP = &xOverlap[3*localId];
//			const double *YP = &yOverlap[3*localId];
//			dx = XP[0]-X[0];
//			dy = XP[1]-X[1];
//			dz = XP[2]-X[2];
//			zeta = sqrt(dx*dx+dy*dy+dz*dz);
//			dx = YP[0]-Y[0];
//			dy = YP[1]-Y[1];
//			dz = YP[2]-Y[2];
//			dY = sqrt(dx*dx+dy*dy+dz*dz);
//			t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * OMEGA * alpha * (dY - zeta));
//			double fx = t * dx / dY;
//			double fy = t * dy / dY;
//			double fz = t * dz / dY;
//
//			*(fOwned+0) += fx*cellVolume;
//			*(fOwned+1) += fy*cellVolume;
//			*(fOwned+2) += fz*cellVolume;
//			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
//			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
//			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
//		}
//
//	}
//}
//
//
//
//
//void computeInternalForceIsotropicElasticPlastic
//(
//		const double* xOverlap,
//		const double* yOverlap,
//		const double* mOwned,
//		const double* volumeOverlap,
//		const double* dilatationOwned,
//		const double* bondDamage,
//		const double* dsfOwned_,
//		const double* deviatoricPlasticExtensionStateN,
//		double* deviatoricPlasticExtensionStateNp1,
//		const double* lambdaN_,
//		double* lambdaNP1_,
//		double* fInternalOverlap,
//		const int*  ownedIDs,
//		const int*  localNeighborList,
//		int numOwnedPoints,
//		double BULK_MODULUS,
//		double SHEAR_MODULUS,
//		double HORIZON,
//		double yieldStress
//)
//{
//
//	/*
//	 * Compute processor local contribution to internal force
//	 */
//	double K = BULK_MODULUS;
//	double MU = SHEAR_MODULUS;
//	double OMEGA=1.0;
//	double DELTA=HORIZON;
//	/*
//	 * 3d variety of yield value
//	 */
//	double yieldValue = 75.0 * yieldStress * yieldStress / 8 / M_PI / pow(DELTA,5);
//	/*
//	 * Planar variety of yield value
//	 */
////		double THICKNESS=1.0;
////		double yieldValue = 225.0 * yieldStress * yieldStress / 8 / M_PI / THICKNESS / pow(DELTA,4);
////		double yieldValue = 0.5 * pow(15*yieldStress/weightedVol,2) * M_PI * THICKNESS * pow(DELTA,4) / 16.0;
//
//
//
//	const double *xOwned = xOverlap;
//	const double *yOwned = yOverlap;
//	const double *m = mOwned;
//	const double *v = volumeOverlap;
//	const double *theta = dilatationOwned;
//	double *fOwned = fInternalOverlap;
//
//	const int *neighPtr = localNeighborList;
//	double cellVolume, alpha, dx, dy, dz, zeta, dY, t, ti, td, ed, edpN, tdTrial;
//	for(int p=0 ; p<numOwnedPoints ; p++){
//
//        int ID = ownedIDs[p];
//        xOwned = &xOverlap[3*ID];
//        yOwned = &yOverlap[3*ID];
//        fOwned = &fInternalOverlap[3*ID];
//        m = &mOwned[ID];
//        theta = &dilatationOwned[ID];
//        const double* lambdaN = &lambdaN_[ID];
//        double* lambdaNP1 = &lambdaNP1_[ID];
//        const double* dsfOwned = &dsfOwned_[ID];
//
//		int numNeigh = *neighPtr; neighPtr++;
//		const double *X = xOwned;
//		const double *Y = yOwned;
//		double weightedVol = *m;
//		alpha = *dsfOwned * 15.0*MU/weightedVol;
//		double selfCellVolume = v[p];
//		double c = 3 * K * (*theta) * OMEGA / weightedVol;
//		double deltaLambda=0.0;
//
//		/*
//		 * Compute norm of trial stress
//		 */
//		double tdNorm = 0.0;
//		tdNorm = MATERIAL_EVALUATION::computeDeviatoricForceStateNorm(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,X,Y,xOverlap,yOverlap,v,alpha,OMEGA);
//
//		/*
//		 * Evaluate yield function
//		 */
//		double pointWiseYieldValue = *dsfOwned * (*dsfOwned) * yieldValue;
//		double f = tdNorm * tdNorm / 2 - pointWiseYieldValue;
//		bool elastic = true;
//
////		std::cout << "Point id = " << p << std::endl;
////		std::cout << "\tyieldStress/m^(4/5) = " << yieldStress/pow(weightedVol,4/5) << std::endl;
////		std::cout << "\tYield Value = " << yieldValue << "; tdNorm * tdNorm / 2 = " << tdNorm * tdNorm / 2 << std::endl;
//		if(f>0){
//			/*
//			 * This step is incrementally plastic
//			 */
//			//			std::cout << "\t PLASTIC" << std::endl;
//			elastic = false;
//			deltaLambda=( tdNorm / sqrt(2.0*pointWiseYieldValue) - 1.0 ) / alpha;
//			*lambdaNP1 = *lambdaN + deltaLambda;
//		} else {
////			std::cout << "\t ELASTIC" << std::endl;
//			*lambdaNP1 = *lambdaN;
//		}
//
//		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++, deviatoricPlasticExtensionStateN++, deviatoricPlasticExtensionStateNp1++){
//			int localId = *neighPtr;
//			cellVolume = v[localId];
//			const double *XP = &xOverlap[3*localId];
//			const double *YP = &yOverlap[3*localId];
//			dx = XP[0]-X[0];
//			dy = XP[1]-X[1];
//			dz = XP[2]-X[2];
//			zeta = sqrt(dx*dx+dy*dy+dz*dz);
//			dx = YP[0]-Y[0];
//			dy = YP[1]-Y[1];
//			dz = YP[2]-Y[2];
//			dY = sqrt(dx*dx+dy*dy+dz*dz);
//			/*
//			 * Deviatoric extension state
//			 */
//			ed = dY-zeta-*theta*zeta/3;
//
//			/*
//			 * Deviatoric plastic extension state from last step
//			 */
//			edpN = *deviatoricPlasticExtensionStateN;
//
//			/*
//			 * Compute trial stress
//			 */
//			tdTrial = alpha * OMEGA * (ed - edpN);
//
//			/*
//			 * Evaluate yield function
//			 */
//			if(elastic){
//				/*
//				 * Elastic case
//				 */
//				td = tdTrial;
//
//				/*
//				 * Therefore edpNp1 = edpN
//				 */
//				*deviatoricPlasticExtensionStateNp1 = *deviatoricPlasticExtensionStateN;
//
//			} else {
//				/*
//				 * Compute deviatoric force state
//				 */
//				td = sqrt(2.0*pointWiseYieldValue) * tdTrial / tdNorm;
//
//				/*
//				 * Update deviatoric plastic deformation state
//				 */
//				*deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;
//
////				std::cout << "Neighbor Id = " << localId << "; Updating deviatoricPlasticExtensionState = " << *deviatoricPlasticExtensionState << std::endl;
//			}
////			std::cout << "\tNeighbor Id = " << localId << "\n\ttd = " << td;
//			/*
//			 * Compute isotropic part of force state
//			 */
//			ti = c * zeta;
//
//			/*
//			 * Force state (with damage)
//			 */
//			double d=(1.0-*bondDamage);
//			t = d*(ti + d*td);
//
//			/*
//			 * Assemble pair wise force function
//			 */
//			double fx = t * dx / dY;
//			double fy = t * dy / dY;
//			double fz = t * dz / dY;
//
//			*(fOwned+0) += fx*cellVolume;
//			*(fOwned+1) += fy*cellVolume;
//			*(fOwned+2) += fz*cellVolume;
//			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
//			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
//			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
//		}
//
//	}
//}


//}
