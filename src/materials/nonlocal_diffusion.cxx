//! \file nonlocal_diffusion.cxx

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
#include "nonlocal_diffusion.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeInternalFluidFlow
(
		const double*  xOverlap,
 		const ScalarT* yOverlap,
		const ScalarT* fluidPressureYOverlap,
		const double* volumeOverlap,
		const double* bondDamage,
		ScalarT* flowInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double isotropicPermeabilityModulus,
		double isotropicPermeabilityPoissons,
		double fluidDensity,
		double baseDynamicViscosity,
		double permeabilityInflectionDamage,
		double permeabilityAlpha,
		double maxPermeability,
    double horizon,
    double ReynoldsThermalViscosityCoefficient,
    const double* deltaTemperature
)
{

	/*
	 * Compute processor local contribution to internal fluid flow
	 */
	double RHO = fluidDensity; //temperature may affect this as well
	// Reynolds exponential model of fluid viscosity dependence on temperature
	// by default the baseDynamicViscosity defined in the input deck is used
	// when no temperature effect is enabled.
	double MU;
	MU = baseDynamicViscosity;

	const double MAX_PERMEABILITY = maxPermeability; 
	const double CRITICAL_DAMAGE = permeabilityInflectionDamage;
	// This number needs to be corrected for the 3D case, it is only valid for 2D
	const double NONLOCAL_PERMEABILITY_CORRESPONDENCE_COEF = .25;
	const double DAMAGE_ALPHA = permeabilityAlpha;

	//const double* xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
	const ScalarT *fluidPressureYOwned = fluidPressureYOverlap;
	const double *v = volumeOverlap;
	double nonlocalPermeability[9];
  const double *deltaT = deltaTemperature;
	ScalarT bondComponents[3];
	double ownedPermeabilityTrace = 0.0; 
	ScalarT effectivePermeability;

	ScalarT *flowOwned = flowInternalOverlap;
	const int *neighPtr = localNeighborList;
	double cellVolume;
	ScalarT dPressure, Y_dx, Y_dy, Y_dz, dY, q;

// Compute the trace of owned permeability
	// tr(K) = K_11 + K_22 + K_33
	// TODO: this is a placeholder. This statement is not physically correct
	ownedPermeabilityTrace = 3*isotropicPermeabilityModulus;	

	for(int p=0;p<numOwnedPoints;p++, fluidPressureYOwned++, yOwned +=3, flowOwned++, deltaT++){
		int numNeigh = *neighPtr; neighPtr++;
		double selfCellVolume = v[p];
		//const double* X = xOwned;
		const ScalarT *Y = yOwned;
		const ScalarT *pressY = fluidPressureYOwned;
		if (deltaTemperature != NULL)
			MU = baseDynamicViscosity*exp(ReynoldsThermalViscosityCoefficient*(*deltaT));
		else
			MU = baseDynamicViscosity;
	
	for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const ScalarT *pressYP = &fluidPressureYOverlap[localId];
			const ScalarT *YP = &yOverlap[3*localId];
			//const double* XP = &xOverlap[3*localId];
			//X_dx = XP[0]-X[0];
			//X_dy = XP[1]-X[1];
			//X_dz = XP[2]-X[2];
			//zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);

			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);

			// We need the actual components of dY organized for indexing by number.
			bondComponents[0] = Y_dx;
			bondComponents[1] = Y_dy;
			bondComponents[2] = Y_dz;
			//omega = scalarInfluenceFunction(zeta,horizon);

			// Pressure potential
			dPressure = pressYP[0]-pressY[0];
			
			// Determine nonlocal permeability
			// This needs to be corrected for 3D diffusion. The NONLOCAL... coefficient was determined
			// with a 2D problem in mind and may not be correct.
			// K_bar = K - coef*tr(K)I
			// TODO: this statement is a placeholder. It is not physically correct
			for(int col=0; col<3 ; ++col){
				for(int row=0; row<3 ; ++row){
					if(row == col)
						nonlocalPermeability[3*col + row] = isotropicPermeabilityModulus - NONLOCAL_PERMEABILITY_CORRESPONDENCE_COEF*ownedPermeabilityTrace;
					else
						nonlocalPermeability[3*col + row] = 0.0;
				}
			}

			//determine the effective permeability in the bond direction
			// xi *dot* K_bar *dot* xi
			effectivePermeability = 0.0;
			for(int row=0; row<3 ; ++row){
				for(int col=0; col<3 ; ++col){
					effectivePermeability += bondComponents[row]*nonlocalPermeability[3*col+row]*bondComponents[col];
				}
			}

			//NOTE: for fluids, bond damage increases flux
			// This non-scientifically represents the idea that permeability is increased by bond damage but there is an
			// upper limit.
			effectivePermeability = effectivePermeability + (MAX_PERMEABILITY - effectivePermeability)*(1.0/(1.0 + exp(DAMAGE_ALPHA*(*bondDamage- CRITICAL_DAMAGE))));

			//TODO: make correct for 3D problem
			//This is not correct for 3D diffusion, the factor of four was derived for 2D and may be different
			q = (RHO / MU) * (4.0 / (PeridigmNS::value_of_pi()*horizon*horizon)) * (effectivePermeability / pow(dY, 4.0)) * dPressure;

			*(flowOwned) += q*cellVolume;
			flowInternalOverlap[localId] -= q*selfCellVolume;
		}
	}

}

/** Explicit template instantiation for double. */
template void computeInternalFluidFlow<double>
(
		const double*  xOverlap,
 		const double* yOverlap,
		const double* fluidPressureYOverlap,
		const double* volumeOverlap,
		const double* bondDamage,
		double* flowInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double isotropicPermeabilityModulus,
		double isotropicPermeabilityPoissons,
		double fluidDensity,
		double baseDynamicViscosity,
		double permeabilityInflectionDamage,
		double permeabilityAlpha,
		double maxPermeability,
    double horizon,
    double ReynoldsThermalViscosityCoefficient,
    const double* deltaTemperature
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalFluidFlow<Sacado::Fad::DFad<double> >
(
		const double*  xOverlap,
 		const Sacado::Fad::DFad<double>* yOverlap,
		const Sacado::Fad::DFad<double>* fluidPressureYOverlap,
		const double* volumeOverlap,
		const double* bondDamage,
		Sacado::Fad::DFad<double>* flowInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double isotropicPermeabilityModulus,
		double isotropicPermeabilityPoissons,
		double fluidDensity,
		double baseDynamicViscosity,
		double permeabilityInflectionDamage,
		double permeabilityAlpha,
		double maxPermeability,
    double horizon,
    double ReynoldsThermalViscosityCoefficient,
    const double* deltaTemperature
);

//! Compute the pressure driven flow.
//! This simple version of the method ignores the lack of pore damage near the node
//! so that a static equilibrium in an isotropic medium can be achieved for diagnosing
//! problems with the proper diffusion model.
template<typename ScalarT>
void computeInternalFluidFlowDeadSimple
(
		const double*  xOverlap,
 		const ScalarT* yOverlap,
		const ScalarT* fluidPressureYOverlap,
		const double* volumeOverlap,
		const double* bondDamage,
		ScalarT* flowInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double isotropicPermeabilityModulus,
		double isotropicPermeabilityPoissons,
		double fluidDensity,
		double baseDynamicViscosity,
		double permeabilityInflectionDamage,
		double permeabilityAlpha,
		double maxPermeability,
    double horizon,
    double ReynoldsThermalViscosityCoefficient,
    const double* deltaTemperature
)
{

	/*
	 * Compute processor local contribution to internal fluid flow
	 */
//	double RHO = fluidDensity; //temperature may affect this as well
	// Reynolds exponential model of fluid viscosity dependence on temperature
	// by default the baseDynamicViscosity defined in the input deck is used
	// when no temperature effect is enabled.
//	double MU;
// if(deltaTemperature != NULL)
//		MU = baseDynamicViscosity*exp(ReynoldsThermalViscosityCoefficient*deltaTemperature[0]);
//	else
//	    MU = baseDynamicViscosity;

//	const double MAX_PERMEABILITY = maxPermeability; 
//	const double CRITICAL_DAMAGE = permeabilityInflectionDamage;
	// This number needs to be corrected for the 3D case, it is only valid for 2D
//	const double NONLOCAL_PERMEABILITY_CORRESPONDENCE_COEF = .25;
//	const double DAMAGE_ALPHA = permeabilityAlpha;

	//const double* xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
	const ScalarT *fluidPressureYOwned = fluidPressureYOverlap;
	const double *v = volumeOverlap;
// double nonlocalPermeability[9];
//  ScalarT bondComponents[3];
//	double ownedPermeabilityTrace = 0.0; 
//	ScalarT effectivePermeability;

	ScalarT *flowOwned = flowInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume;
	ScalarT dPressure, Y_dx, Y_dy, Y_dz, dY, q;

	// Compute the trace of owned permeability
	// tr(K) = K_11 + K_22 + K_33
	// TODO: this is a placeholder. This statement is not physically correct
//	ownedPermeabilityTrace = 3*isotropicPermeabilityModulus;	

	for(int p=0;p<numOwnedPoints;p++, fluidPressureYOwned++, yOwned +=3, flowOwned++){
		int numNeigh = *neighPtr; neighPtr++;
		double selfCellVolume = v[p];
		//const double* X = xOwned;
		const ScalarT *Y = yOwned;
		const ScalarT *pressY = fluidPressureYOwned;

	for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const ScalarT *pressYP = &fluidPressureYOverlap[localId];
			const ScalarT *YP = &yOverlap[3*localId];
			//const double* XP = &xOverlap[3*localId];
			//X_dx = XP[0]-X[0];
			//X_dy = XP[1]-X[1];
			//X_dz = XP[2]-X[2];
			//zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);

			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);

			// We need the actual components of dY organized for indexing by number.
			// bondComponents[0] = Y_dx;
			// bondComponents[1] = Y_dy;
			// bondComponents[2] = Y_dz;
			//omega = scalarInfluenceFunction(zeta,horizon);

			// Pressure potential
			dPressure = pressYP[0]-pressY[0];
			
			// Determine nonlocal permeability
			// This needs to be corrected for 3D diffusion. The NONLOCAL... coefficient was determined
			// with a 2D problem in mind and may not be correct.
			// K_bar = K - coef*tr(K)I
			// TODO: this statement is a placeholder. It is not physically correct
			/*
			for(int col=0; col<3 ; ++col){
				for(int row=0; row<3 ; ++row){
					if(row == col)
						nonlocalPermeability[3*col + row] = isotropicPermeabilityModulus - NONLOCAL_PERMEABILITY_CORRESPONDENCE_COEF*ownedPermeabilityTrace;
					else
						nonlocalPermeability[3*col + row] = 0.0;
				}
			}
			*/

			//determine the effective permeability in the bond direction
			// xi *dot* K_bar *dot* xi
			/*
			effectivePermeability = 0.0;
			for(int row=0; row<3 ; ++row){
				for(int col=0; col<3 ; ++col){
					effectivePermeability += bondComponents[row]*nonlocalPermeability[3*col+row]*bondComponents[col];
				}
			}
			*/

			//NOTE: Here, permeability not at all affected by bond damage
			// This non-scientifically represents the idea that permeability is increased by bond damage but there is an
			// upper limit.
			//effectivePermeability = effectivePermeability + (MAX_PERMEABILITY - effectivePermeability)*(1.0/(1.0 + exp(DAMAGE_ALPHA*(*bondDamage- CRITICAL_DAMAGE))));

			//TODO: make correct for 3D problem
			//This is not correct for 3D diffusion, the factor of four was derived for 2D and may be different
			//q = (RHO / MU) * (4.0 / (PeridigmNS::value_of_pi()*horizon*horizon)) * (effectivePermeability / pow(dY, 4.0)) * dPressure;
		  //q = (RHO / MU) * (4.0 / (PeridigmNS::value_of_pi()*horizon*horizon)) * (1.0 / pow(dY, 4.0)) * dPressure;
		  q = (horizon*horizon) * (1.0 / pow(dY, 4.0)) * dPressure;

			//NOTE: when pressure outside of a node is greater than pressure inside the node,
			// flow is positive, that is, fluid flows into the node and builds pressure.
			*(flowOwned) += q*cellVolume;
			flowInternalOverlap[localId] -= q*selfCellVolume;
		}
	}
}

/** Explicit template instantiation for double. */
template void computeInternalFluidFlowDeadSimple<double>
(
		const double*  xOverlap,
 		const double* yOverlap,
		const double* fluidPressureYOverlap,
		const double* volumeOverlap,
		const double* bondDamage,
		double* flowInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double isotropicPermeabilityModulus,
		double isotropicPermeabilityPoissons,
		double fluidDensity,
		double baseDynamicViscosity,
		double permeabilityInflectionDamage,
		double permeabilityAlpha,
		double maxPermeability,
    double horizon,
    double ReynoldsThermalViscosityCoefficient,
    const double* deltaTemperature
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalFluidFlowDeadSimple<Sacado::Fad::DFad<double> >
(
		const double*  xOverlap,
 		const Sacado::Fad::DFad<double>* yOverlap,
		const Sacado::Fad::DFad<double>* fluidPressureYOverlap,
		const double* volumeOverlap,
		const double* bondDamage,
		Sacado::Fad::DFad<double>* flowInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double isotropicPermeabilityModulus,
		double isotropicPermeabilityPoissons,
		double fluidDensity,
		double baseDynamicViscosity,
		double permeabilityInflectionDamage,
		double permeabilityAlpha,
		double maxPermeability,
    double horizon,
    double ReynoldsThermalViscosityCoefficient,
    const double* deltaTemperature
);

//! Computes contributions to the internal force resulting from owned points.
template<typename ScalarT>
void computeInternalForceLinearElasticCoupled
(
		const double* xOverlap,
		const ScalarT* yOverlap,
		const ScalarT* fluidPressureYOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const ScalarT* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		ScalarT* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
    double horizon,
    double thermalExpansionCoefficient,
    const double* deltaTemperature
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;

	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
	const ScalarT *fluidPressureYOwned = fluidPressureYOverlap;
  const double *deltaT = deltaTemperature;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const double *dsf = dsfOwned;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, X_dx, X_dy, X_dz, zeta, omega;
	ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e, c1;
	for(int p=0;p<numOwnedPoints;p++, fluidPressureYOwned++, xOwned +=3, yOwned +=3, fOwned+=3, deltaT++, m++, theta++, dsf++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		alpha = 15.0*MU/(*m);
		alpha *= (*dsf);
		double selfCellVolume = v[p];

		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
      e = dY - zeta;
      if(deltaTemperature)
      	e -= thermalExpansionCoefficient*(*deltaT)*zeta;

			omega = scalarInfluenceFunction(zeta,horizon);
			// c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
			// NOTE: set pressure effect to maximum, this is not to be a permanent change
			c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0) -3.0*omega/(*m)*(1.0 /**bondDamage*/)*(*fluidPressureYOwned);
			t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
			fx = t * Y_dx / dY;
			fy = t * Y_dy / dY;
			fz = t * Y_dz / dY;

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
template void computeInternalForceLinearElasticCoupled<double>
(
		const double* xOverlap,
		const double* yOverlap,
		const double* fluidPressureYOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		double* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
    double horizon,
    double thermalExpansionCoefficient,
    const double* deltaTemperature
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceLinearElasticCoupled<Sacado::Fad::DFad<double> >
(
		const double* xOverlap,
		const Sacado::Fad::DFad<double>* yOverlap,
		const Sacado::Fad::DFad<double>* fluidPressureYOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const Sacado::Fad::DFad<double>* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		Sacado::Fad::DFad<double>* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
    double horizon,
    double thermalExpansionCoefficient,
    const double* deltaTemperature
);

//! Computes contributions to the internal force resulting from owned points.
// In this simple version of the method, fluid pressure at a node always affects the 
// dilatation at a node regardless of the lack of bond damage near the node.
template<typename ScalarT>
void computeInternalForceLinearElasticCoupledDeadSimple
(
		const double* xOverlap,
		const ScalarT* yOverlap,
		const ScalarT* fluidPressureYOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const ScalarT* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		ScalarT* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
    double horizon,
    double thermalExpansionCoefficient,
    const double* deltaTemperature
)
{

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;

	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
	const ScalarT *fluidPressureYOwned = fluidPressureYOverlap;
  const double *deltaT = deltaTemperature;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const double *dsf = dsfOwned;
	const ScalarT *theta = dilatationOwned;
	ScalarT *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, alpha, X_dx, X_dy, X_dz, zeta, omega;
	ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e, c1;
	for(int p=0;p<numOwnedPoints;p++, fluidPressureYOwned++, xOwned +=3, yOwned +=3, fOwned+=3, deltaT++, m++, theta++, dsf++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		alpha = 15.0*MU/(*m);
		alpha *= (*dsf);
		double selfCellVolume = v[p];
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
			X_dx = XP[0]-X[0];
			X_dy = XP[1]-X[1];
			X_dz = XP[2]-X[2];
			zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
			Y_dx = YP[0]-Y[0];
			Y_dy = YP[1]-Y[1];
			Y_dz = YP[2]-Y[2];
			dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
            e = dY - zeta;
            if(deltaTemperature)
              e -= thermalExpansionCoefficient*(*deltaT)*zeta;
			omega = scalarInfluenceFunction(zeta,horizon);
			// c1 = omega*(*theta)*(9.0*K-15.0*MU)/(3.0*(*m));
			//c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0) -3.0*omega/(*m)*(*bondDamage)*(*fluidPressureYOwned);
			//NOTE: Notice how regardless of bond damage fluid pressure has an effect on dilatation.
			//c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0) -3.0*omega/(*m)*(*fluidPressureYOwned);
			c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
			t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
			fx = t * Y_dx / dY;
			fy = t * Y_dy / dY;
			fz = t * Y_dz / dY;

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
template void computeInternalForceLinearElasticCoupledDeadSimple<double>
(
		const double* xOverlap,
		const double* yOverlap,
		const double* fluidPressureYOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		double* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
    double horizon,
    double thermalExpansionCoefficient,
    const double* deltaTemperature
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceLinearElasticCoupledDeadSimple<Sacado::Fad::DFad<double> >
(
		const double* xOverlap,
		const Sacado::Fad::DFad<double>* yOverlap,
		const Sacado::Fad::DFad<double>* fluidPressureYOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const Sacado::Fad::DFad<double>* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned,
		Sacado::Fad::DFad<double>* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
    double horizon,
    double thermalExpansionCoefficient,
    const double* deltaTemperature
);


}
