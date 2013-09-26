//! \file material_utilities.cxx

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

#include "material_utilities.h"
#include <cmath>
#include <vector>
#include <Sacado.hpp>

namespace MATERIAL_EVALUATION {

/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 * NOTE: this function will overwrite entries in 'yOverlap'
 * for all of 'X' neighbors
 * OUTPUT: yOverlap such that there is a state of pure
 * shear at 'X'
 */

void set_pure_shear
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		double *yOverlap,
		PURE_SHEAR mode,
		double gamma
)
{

	/*
	 * Pure shear centered at X
	 * X has no displacement
	 */
	int numNeigh=*neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++, neighPtr++){
		int localId = *neighPtr;
		const double *XP = &xOverlap[3*localId];

		double dx = XP[0]-X[0];
		double dy = XP[1]-X[1];
		double dz = XP[2]-X[2];

		/*
		 * Pure shear
		 */
		double zx(0.0), xy(0.0), yz(0.0);
		double xz(0.0), yx(0.0), zy(0.0);
		switch(mode){
		case ZX:
			zx = gamma * dx;
			xz = gamma * dz;
			break;
		case XY:
			xy = gamma * dy;
			yx = gamma * dx;
			break;
		case YZ:
			yz = gamma * dz;
			zy = gamma * dy;
			break;
		}

		double *YP = &yOverlap[3*localId];
		YP[0] = XP[0] + xy + xz;
		YP[1] = XP[1] + yz + yx;
		YP[2] = XP[2] + zx + zy;

	}

}

double scalarInfluenceFunction
(
        double zeta,
        double horizon
)
{
  static PeridigmNS::InfluenceFunction::functionPointer influenceFunction = NULL;
  if(!influenceFunction)
    influenceFunction = PeridigmNS::InfluenceFunction::self().getInfluenceFunction();

  return influenceFunction(zeta, horizon);
}

double computeWeightedVolume
(
		const double *X,
		const double *xOverlap,
		const double* volumeOverlap,
		const int* localNeighborList,
        double horizon,
		const FunctionPointer omega
){

	double m=0.0;
	double cellVolume;
	const int *neighPtr = localNeighborList;
	int numNeigh = *neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++,neighPtr++){
		int localId = *neighPtr;
		cellVolume = volumeOverlap[localId];
		const double *XP = &xOverlap[3*localId];
		double dx = XP[0]-X[0];
		double dy = XP[1]-X[1];
		double dz = XP[2]-X[2];
		double zetaSquared = dx*dx+dy*dy+dz*dz;
		double d = sqrt(zetaSquared);
		m+=omega(d,horizon)*(zetaSquared)*cellVolume;
	}

	return m;
}

void computeDeviatoricDilatation
(
		const double* xOverlap,
		const double* yOverlap,
		const double *mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		const double* epd,
		double* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints,
		double horizon,
		const FunctionPointer OMEGA
)
{
	const double *xOwned = xOverlap;
	const double *yOwned = yOverlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	double *theta = dilatationOwned;
	double cellVolume;
	const int *neighPtr = localNeighborList;
	for(int p=0; p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, m++, theta++){
		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		//const double *Y = yOwned;
		*theta = 0;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,epd++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			//const double *YP = &yOverlap[3*localId];
			double dx = XP[0]-X[0];
			double dy = XP[1]-X[1];
			double dz = XP[2]-X[2];
			double zetaSquared = dx*dx+dy*dy+dz*dz;
			double d = sqrt(zetaSquared);
            double omega = OMEGA(d,horizon);
			double e = (*epd);
			*theta += 3.0*omega*(1.0-*bondDamage)*d*e*cellVolume/(*m);
		}

	}
}

template<typename ScalarT>
void computeDilatation
(
		const double* xOverlap,
		const ScalarT* yOverlap,
		const double *mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		ScalarT* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints,
        double horizon,
		const FunctionPointer OMEGA,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
)
{
	const double *xOwned = xOverlap;
	const ScalarT *yOwned = yOverlap;
    const double *deltaT = deltaTemperature;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	ScalarT *theta = dilatationOwned;
	double cellVolume;
	const int *neighPtr = localNeighborList;
	for(int p=0; p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, deltaT++, m++, theta++){
		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const ScalarT *Y = yOwned;
		*theta = ScalarT(0.0);
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const ScalarT *YP = &yOverlap[3*localId];
			double X_dx = XP[0]-X[0];
			double X_dy = XP[1]-X[1];
			double X_dz = XP[2]-X[2];
			double zetaSquared = X_dx*X_dx+X_dy*X_dy+X_dz*X_dz;
			ScalarT Y_dx = YP[0]-Y[0];
			ScalarT Y_dy = YP[1]-Y[1];
			ScalarT Y_dz = YP[2]-Y[2];
			ScalarT dY = Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz;
			double d = sqrt(zetaSquared);
			ScalarT e = sqrt(dY)-d;
            if(deltaTemperature)
              e -= thermalExpansionCoefficient*(*deltaT)*d;
            double omega = OMEGA(d,horizon);
			*theta += 3.0*omega*(1.0-*bondDamage)*d*e*cellVolume/(*m);
		}

	}
}

/** Explicit template instantiation for double. */
template
void computeDilatation<double>
(
		const double* xOverlap,
		const double* yOverlap,
		const double *mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		double* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints,
        double horizon,
		const FunctionPointer OMEGA,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
 );


/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template
void computeDilatation<Sacado::Fad::DFad<double> >
(
		const double* xOverlap,
		const Sacado::Fad::DFad<double>* yOverlap,
		const double *mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		Sacado::Fad::DFad<double>* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints,
        double horizon,
		const FunctionPointer OMEGA,
        double thermalExpansionCoefficient,
        const double* deltaTemperature
 );

/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 */
double computeDilatation
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		const double *Y,
		const double *yOverlap,
		const double *volumeOverlap,
		double weightedVolume,
		double horizon,
		const FunctionPointer OMEGA
)
{
	double bondDamage=0.0;
	const double *v = volumeOverlap;
	double m = weightedVolume;
	double theta = 0.0;
	int numNeigh=*neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++,neighPtr++){

		int localId = *neighPtr;
		double cellVolume = v[localId];

		const double *XP = &xOverlap[3*localId];
		double dx = XP[0]-X[0];
		double dy = XP[1]-X[1];
		double dz = XP[2]-X[2];
		double zetaSquared = dx*dx+dy*dy+dz*dz;

		const double *YP = &yOverlap[3*localId];
		dx = YP[0]-Y[0];
		dy = YP[1]-Y[1];
		dz = YP[2]-Y[2];
		double dY = dx*dx+dy*dy+dz*dz;
		double d = sqrt(zetaSquared);
        double omega = OMEGA(d,horizon);
		double e = sqrt(dY)-d;
		theta += 3.0*omega*(1.0-bondDamage)*d*e*cellVolume/m;

	}
	return theta;
}


double compute_norm_2_deviatoric_extension
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		const double *Y,
		const double *yOverlap,
		const double *volumeOverlap,
		double weighted_volume,
		double horizon,
		const FunctionPointer OMEGA
)
{

	const double *v = volumeOverlap;
	double cellVolume, dx, dy, dz, zeta, dY, ed;

	/*
	 * Compute weighted volume
	 */
	double m = weighted_volume;
//	double m = computeWeightedVolume(X,xOverlap,volumeOverlap,neighPtr);
//	std::cout << NAMESPACE << "probeShearModulusScaleFactor weighted volume = " << m << std::endl;

	/*
	 * Compute dilatation
	 */
	double theta = computeDilatation(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m,horizon);
//	std::cout << NAMESPACE << "probeShearModulusScaleFactor theta = " << theta << std::endl;

	/*
	 * Pure shear centered at X
	 * X has no displacement
	 */
	double ed_squared=0.0;
	int numNeigh=*neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++, neighPtr++){
		int localId = *neighPtr;
		cellVolume = v[localId];
		const double *XP = &xOverlap[3*localId];
		dx = XP[0]-X[0];
		dy = XP[1]-X[1];
		dz = XP[2]-X[2];
		zeta = sqrt(dx*dx+dy*dy+dz*dz);

		/*
		 * Deformation State
		 */
		const double *YP = &yOverlap[3*localId];
		dx = YP[0]-Y[0];
		dy = YP[1]-Y[1];
		dz = YP[2]-Y[2];
		dY = sqrt(dx*dx+dy*dy+dz*dz);

		/*
		 * Deviatoric extension state
		 */
		ed = dY-zeta-theta*zeta/3;

		/*
		 * Accumulate norm
		 */
		double omega=OMEGA(zeta,horizon);
		ed_squared += ed * omega * ed * cellVolume;

	}

	return ed_squared;
}


void computeShearCorrectionFactor
(
        int numOwnedPoints,
        int lengthYOverlap,
		const double *xOverlap,
		double *yOverlap,
		const double *volumeOverlap,
		const double *owned_weighted_volume,
		const int*  localNeighborList,
		double horizon,
		double *shearCorrectionFactorOwned
){
	double gamma=1.0e-6;
	// currently un-used but may be helpful in further studies
//  double reference = 4.0 * M_PI * gamma * gamma * pow(horizon,5) / 75.0;
	const int *neighPtr = localNeighborList;
	const double *xOwned = xOverlap;
	double *yOwned = yOverlap;
	double *scaleFactor = shearCorrectionFactorOwned;
	PURE_SHEAR mode;
	for(int p=0;p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, scaleFactor++, owned_weighted_volume++){
		int numNeigh = *neighPtr;
		const double *X = xOwned;
		double *Y = yOwned;
		double m = *owned_weighted_volume;
		double scf, max_dsf;

		mode = XY;
		Y[0] = X[0]; Y[1] = X[1]; Y[2] = X[2];
		set_pure_shear(neighPtr,xOwned,xOverlap,yOverlap,mode,gamma);
		scf=compute_norm_2_deviatoric_extension(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m,horizon);
		max_dsf=scf;

		mode = ZX;
		Y[0] = X[0]; Y[1] = X[1]; Y[2] = X[2];
		set_pure_shear(neighPtr,xOwned,xOverlap,yOverlap,mode,gamma);
		scf=compute_norm_2_deviatoric_extension(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m,horizon);
		if(scf>max_dsf) max_dsf=scf;

		mode = YZ;
		Y[0] = X[0]; Y[1] = X[1]; Y[2] = X[2];
		set_pure_shear(neighPtr,xOwned,xOverlap,yOverlap,mode,gamma);
		scf=compute_norm_2_deviatoric_extension(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m,horizon);
		if(scf>max_dsf) max_dsf=scf;

		// We need to reset the location before
		//   iterating to the next point; OTHERWISE its a BUG!
		Y[0] = X[0]; Y[1] = X[1]; Y[2] = X[2];
		scf=max_dsf;
		*scaleFactor = 4.0 * gamma * gamma  * m / scf /15.0;
		neighPtr+=(numNeigh+1);
	}
}

void computeWeightedVolume
(
		const double* xOverlap,
		const double* volumeOverlap,
		double *mOwned,
		int myNumPoints,
		const int* localNeighborList,
        double horizon,
        const FunctionPointer OMEGA
){
	double *m = mOwned;
	const double *xOwned = xOverlap;
	const int *neighPtr = localNeighborList;
	for(int p=0;p<myNumPoints;p++, xOwned+=3, m++){
		int numNeigh = *neighPtr;
		const double *X = xOwned;
		*m=MATERIAL_EVALUATION::computeWeightedVolume(X,xOverlap,volumeOverlap,neighPtr,horizon);
		neighPtr+=(numNeigh+1);
	}
}


namespace WITH_BOND_VOLUME {

/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 */

double computeWeightedVolume
(
		const double *X,
		const double *xOverlap,
		const double* bondVolume,
		const int* localNeighborList,
		double horizon,
		const FunctionPointer OMEGA
){

	double m=0.0;
	const int *neighPtr = localNeighborList;
	const double *bond_volume = bondVolume;
	int numNeigh = *neighPtr; neighPtr++;
//	std::cout << NAMESPACE <<"computeWeightedVolume\n";
//	std::cout << "\tnumber of neighbors = " << numNeigh << std::endl;
	for(int n=0;n<numNeigh;n++,neighPtr++,bond_volume++){
		int localId = *neighPtr;
		const double *XP = &xOverlap[3*localId];
		double dx = XP[0]-X[0];
		double dy = XP[1]-X[1];
		double dz = XP[2]-X[2];
        double zetaSquared = dx*dx+dy*dy+dz*dz;
        double d = sqrt(zetaSquared);
        double omega = OMEGA(d,horizon);
		m+=omega*(zetaSquared)*(*bond_volume);
	}

	return m;
}


/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 * NOTE: bondVolume is layed out like the neighborhood list; length
 * of bondVolume is numNeigh
 */

double computeDilatation
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		const double *Y,
		const double *yOverlap,
		const double *bondVolume,
		double weightedVolume,
		double horizon,
		const FunctionPointer OMEGA
)
{
	double bondDamage=0.0;
	double m = weightedVolume;
	double theta = 0.0;
	int numNeigh=*neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++,neighPtr++, bondVolume++){

		int localId = *neighPtr;

		const double *XP = &xOverlap[3*localId];
		double dx = XP[0]-X[0];
		double dy = XP[1]-X[1];
		double dz = XP[2]-X[2];
		double zetaSquared = dx*dx+dy*dy+dz*dz;

		const double *YP = &yOverlap[3*localId];
		dx = YP[0]-Y[0];
		dy = YP[1]-Y[1];
		dz = YP[2]-Y[2];
		double dY = dx*dx+dy*dy+dz*dz;
		double d = sqrt(zetaSquared);
        double omega = OMEGA(d,horizon);
		double e = sqrt(dY)-d;
		theta += 3.0*omega*(1.0-bondDamage)*d*e*(*bondVolume)/m;

	}
	return theta;
}

double compute_norm_2_deviatoric_extension
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		const double *Y,
		const double *yOverlap,
		const double *bondVolume,
		double weighted_volume,
		double horizon,
		const FunctionPointer OMEGA
)
{

	double dx, dy, dz, zeta, dY, ed;

	/*
	 * Compute weighted volume
	 */
	double m = weighted_volume;
//	double m = computeWeightedVolume(X,xOverlap,volumeOverlap,neighPtr);
//	std::cout << NAMESPACE << "probeShearModulusScaleFactor weighted volume = " << m << std::endl;

	/*
	 * Compute dilatation
	 */
	double theta = computeDilatation(neighPtr,X,xOverlap,Y,yOverlap,bondVolume,m,horizon);
//	std::cout << NAMESPACE << "probeShearModulusScaleFactor theta = " << theta << std::endl;

	/*
	 * Pure shear centered at X
	 * X has no displacement
	 */
	double ed_squared=0.0;
	int numNeigh=*neighPtr; neighPtr++;
	for(int n=0;n<numNeigh;n++, neighPtr++, bondVolume++){
		int localId = *neighPtr;

		const double *XP = &xOverlap[3*localId];
		dx = XP[0]-X[0];
		dy = XP[1]-X[1];
		dz = XP[2]-X[2];
		zeta = sqrt(dx*dx+dy*dy+dz*dz);

		/*
		 * Deformation State
		 */
		const double *YP = &yOverlap[3*localId];
		dx = YP[0]-Y[0];
		dy = YP[1]-Y[1];
		dz = YP[2]-Y[2];
		dY = sqrt(dx*dx+dy*dy+dz*dz);

		/*
		 * Deviatoric extension state
		 */
		ed = dY-zeta-theta*zeta/3;

		/*
		 * Accumulate norm
		 */
		double omega=OMEGA(zeta,horizon);
		ed_squared += ed * omega * ed * (*bondVolume);

	}

	return ed_squared;
}

}


}
