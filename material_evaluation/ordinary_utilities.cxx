
/*! \file ordinary_utilities.cxx */

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

#include "ordinary_utilities.h"
#include <cmath>

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
		switch(mode){
		case ZX:
			zx = gamma * dx;
			break;
		case XY:
			xy = gamma * dy;
			break;
		case YZ:
			yz = gamma * dz;
			break;
		}

		double *YP = &yOverlap[3*localId];
		YP[0] = XP[0] + xy;
		YP[1] = XP[1] + yz;
		YP[2] = XP[2] + zx;

	}

}

double computeWeightedVolume
(
		const double *X,
		const double *xOverlap,
		const double* volumeOverlap,
		const int* localNeighborList
){

	double m=0.0;
	double cellVolume;
	const int *neighPtr = localNeighborList;
	int numNeigh = *neighPtr; neighPtr++;
//	std::cout << NAMESPACE <<"computeWeightedVolume\n";
//	std::cout << "\tnumber of neighbors = " << numNeigh << std::endl;
	for(int n=0;n<numNeigh;n++,neighPtr++){
		int localId = *neighPtr;
		cellVolume = volumeOverlap[localId];
		const double *XP = &xOverlap[3*localId];
		double dx = XP[0]-X[0];
		double dy = XP[1]-X[1];
		double dz = XP[2]-X[2];
		m+=(dx*dx+dy*dy+dz*dz)*cellVolume;
	}

	return m;
}

void computeDilatation
(
		const double* xOverlap,
		const double* yOverlap,
		const double *mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		double* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints
)
{
	double OMEGA=1.0;
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
		const double *Y = yOwned;
		*theta = 0;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP = &xOverlap[3*localId];
			const double *YP = &yOverlap[3*localId];
			double dx = XP[0]-X[0];
			double dy = XP[1]-X[1];
			double dz = XP[2]-X[2];
			double zetaSqared = dx*dx+dy*dy+dz*dz;
			dx = YP[0]-Y[0];
			dy = YP[1]-Y[1];
			dz = YP[2]-Y[2];
			double dY = dx*dx+dy*dy+dz*dz;
			double d = sqrt(zetaSqared);
			double e = sqrt(dY)-d;
			*theta += 3.0*OMEGA*(1.0-*bondDamage)*d*e*cellVolume/(*m);
		}

	}
}

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
		double weightedVolume
)
{
	double OMEGA=1.0;
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
		double zetaSqared = dx*dx+dy*dy+dz*dz;

		const double *YP = &yOverlap[3*localId];
		dx = YP[0]-Y[0];
		dy = YP[1]-Y[1];
		dz = YP[2]-Y[2];
		double dY = dx*dx+dy*dy+dz*dz;
		double d = sqrt(zetaSqared);
		double e = sqrt(dY)-d;
		theta += 3.0*OMEGA*(1.0-bondDamage)*d*e*cellVolume/m;

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
		double weighted_volume
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
	double theta = computeDilatation(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m);
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
		ed_squared += ed * ed * cellVolume;

	}

	return ed_squared;
}


void computeShearCorrectionFactor
(
		int numOwnedPoints,
		const double *xOverlap,
		double *yOverlap_scratch_required_work_space,
		const double *volumeOverlap,
		const double *owned_weighted_volume,
		const int*  localNeighborList,
		double horizon,
		double *shearCorrectionFactorOwned
){
	double gamma=1.0e-6;
	double reference = 4.0 * M_PI * gamma * gamma * pow(horizon,5) / 75.0;
	const int *neighPtr = localNeighborList;
	const double *xOwned = xOverlap;
	const double *yOwned = yOverlap_scratch_required_work_space;
	double *yOverlap = yOverlap_scratch_required_work_space;
	double *scaleFactor = shearCorrectionFactorOwned;
	PURE_SHEAR mode;
	for(int p=0;p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, scaleFactor++, owned_weighted_volume++){
		int numNeigh = *neighPtr;
		const double *X = xOwned;
		const double *Y = yOwned;
		double m = *owned_weighted_volume;
		double dsf, max_dsf;

		mode = XY;
		set_pure_shear(neighPtr,xOwned,xOverlap,yOverlap,mode,gamma);
		dsf=compute_norm_2_deviatoric_extension(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m);
		max_dsf=dsf;

		mode = ZX;
		set_pure_shear(neighPtr,xOwned,xOverlap,yOverlap,mode,gamma);
		dsf=compute_norm_2_deviatoric_extension(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m);
		if(dsf>max_dsf) max_dsf = dsf;

		mode = YZ;
		set_pure_shear(neighPtr,xOwned,xOverlap,yOverlap,mode,gamma);
		dsf=compute_norm_2_deviatoric_extension(neighPtr,X,xOverlap,Y,yOverlap,volumeOverlap,m);
		if(dsf>max_dsf) max_dsf = dsf;

		/*
		 * Better guard against division by zero
		 *
		 */
		double tolerance(1.0e-15);
		if(max_dsf/reference < tolerance)
			max_dsf=1.0;
		else
			max_dsf = reference/max_dsf;

		*scaleFactor = max_dsf;
		neighPtr+=(numNeigh+1);
	}

}

void computeWeightedVolume
(
		const double* xOverlap,
		const double* volumeOverlap,
		double *mOwned,
		int myNumPoints,
		const int* localNeighborList
){
	double *m = mOwned;
	const double *xOwned = xOverlap;
	const int *neighPtr = localNeighborList;
	for(int p=0;p<myNumPoints;p++, xOwned+=3, m++){
		int numNeigh = *neighPtr;
		const double *X = xOwned;
		*m=MATERIAL_EVALUATION::computeWeightedVolume(X,xOverlap,volumeOverlap,neighPtr);
		neighPtr+=(numNeigh+1);
	}
}


}
