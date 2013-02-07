//! \file elastic_plastic.cxx

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
#include "elastic_plastic.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
ScalarT computeDeviatoricForceStateNorm
(
		int numNeigh,
		ScalarT theta,
		const int *neighPtr,
		const double *bondDamage,
		const double *deviatoricPlasticExtensionState,
		const double *X,
		const ScalarT *Y,
		const double *xOverlap,
		const ScalarT *yOverlap,
		const double *volumeOverlap,
		double alpha,
		double OMEGA
)
{
	ScalarT norm=0.0;
	const double *v = volumeOverlap;
	double cellVolume, dx_X, dy_X, dz_X, zeta, edpN;
    ScalarT dx_Y, dy_Y, dz_Y, dY, ed, tdTrial;

	for(int n=0;n<numNeigh;n++, neighPtr++, bondDamage++, deviatoricPlasticExtensionState++){
		int localId = *neighPtr;
		cellVolume = v[localId];
		const double *XP = &xOverlap[3*localId];
		const ScalarT *YP = &yOverlap[3*localId];
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
		ed = dY-zeta-theta*zeta/3;

		/*
		 * Deviatoric plastic extension state from last step
		 */
		edpN = *deviatoricPlasticExtensionState;

		/*
		 * Compute trial stress
		 * NOTE: include damage
		 */
		double d=(1.0-*bondDamage);
		tdTrial = d * alpha * OMEGA * (ed - edpN);

		/*
		 * Accumulate norm
		 */
		norm += tdTrial * tdTrial * cellVolume;
	}

	return sqrt(norm);
}

template<typename ScalarT>
void computeInternalForceIsotropicElasticPlastic
(
		const double* xOverlap,
		const ScalarT* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const ScalarT* dilatationOwned,
		const double* bondDamage,
		const double* scfOwned,
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
		bool isPlanarProblem,
		double thickness
)
{
	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;
	double OMEGA=1.0;
	double DELTA=HORIZON;
	double THICKNESS=thickness;
	/*
	 * 2d or 3d variety of yield value (uniaxial stress)
	 */
    double yieldValue = 25.0 * yieldStress * yieldStress / 8 / M_PI / pow(DELTA,5);
	if(isPlanarProblem)
    	double yieldValue = 225.0 / 3. * yieldStress * yieldStress / 8 / M_PI / THICKNESS / pow(DELTA,4);


	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yNP1Overlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, dx_X, dy_X, dz_X, zeta, edpN;
    ScalarT dx_Y, dy_Y, dz_Y, dY, ed, tdTrial, t, ti, td;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, m++, theta++, lambdaN++, lambdaNP1++, scfOwned++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		double weightedVol = *m;
		alpha = *scfOwned * 15.0*MU/weightedVol;
		double selfCellVolume = v[p];
		ScalarT c = 3 * K * (*theta) * OMEGA / weightedVol;
		ScalarT deltaLambda=0.0;

		/*
		 * Compute norm of trial stress
		 */
		ScalarT tdNorm = 0.0;
		tdNorm = computeDeviatoricForceStateNorm(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,X,Y,xOverlap,yNP1Overlap,v,alpha,OMEGA);

		/*
		 * Evaluate yield function
		 */
		double pointWiseYieldValue = *scfOwned * (*scfOwned) * yieldValue;
		ScalarT f = tdNorm * tdNorm / 2 - pointWiseYieldValue;
		bool elastic = true;

		if(f>0){
			/*
			 * This step is incrementally plastic
			 */
			//			std::cout << "\t PLASTIC" << std::endl;
			elastic = false;
			deltaLambda=( tdNorm / sqrt(2.0*pointWiseYieldValue) - 1.0 ) / alpha;
			*lambdaNP1 = *lambdaN + deltaLambda;
		} else {
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
				td = sqrt(2.0*pointWiseYieldValue) * tdTrial / tdNorm;

				/*
				 * Update deviatoric plastic deformation state
				 */
				*deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;
			}
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

/** Explicit template instantiation for double. */
template double computeDeviatoricForceStateNorm<double>
(
		int numNeigh,
		double theta,
		const int *neighPtr,
		const double *bondDamage,
		const double *deviatoricPlasticExtensionState,
		const double *X,
		const double *Y,
		const double *xOverlap,
		const double *yOverlap,
		const double *volumeOverlap,
		double alpha,
		double OMEGA
);

/** Explicit template instantiation for double. */
template void computeInternalForceIsotropicElasticPlastic<double>
(
		const double* xOverlap,
		const double* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* scfOwned,
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
		bool isPlanarProblem,
		double thickness
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template Sacado::Fad::DFad<double> computeDeviatoricForceStateNorm<Sacado::Fad::DFad<double> >
(
		int numNeigh,
		Sacado::Fad::DFad<double> theta,
		const int *neighPtr,
		const double *bondDamage,
		const double *deviatoricPlasticExtensionState,
		const double *X,
		const Sacado::Fad::DFad<double> *Y,
		const double *xOverlap,
		const Sacado::Fad::DFad<double> *yOverlap,
		const double *volumeOverlap,
		double alpha,
		double OMEGA
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceIsotropicElasticPlastic<Sacado::Fad::DFad<double> >
(
		const double* xOverlap,
		const Sacado::Fad::DFad<double>* yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const Sacado::Fad::DFad<double>* dilatationOwned,
		const double* bondDamage,
		const double* scfOwned,
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
		bool isPlanarProblem,
		double thickness
);

}

