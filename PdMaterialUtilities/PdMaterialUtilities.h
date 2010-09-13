// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions?
// David J. Littlewood   djlittl@sandia.gov
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#ifndef PDMATERIALUTILITIES_H_
#define PDMATERIALUTILITIES_H_
namespace PdMaterialUtilities {

//void computeDilatation
//(
//		std::tr1::shared_ptr<double> xOverlapPtr,
//		std::tr1::shared_ptr<double> yOverlapPtr,
//		double *mOwned,
//		std::tr1::shared_ptr<double> volumeOverlapPtr,
//		double* dilatationOwned,
//		std::tr1::shared_ptr<int> localNeighborList,
//		int numOwnedPoints
//);

//void computeInternalForce
//(
//		std::tr1::shared_ptr<double> xOverlapPtr,
//		std::tr1::shared_ptr<double> yOverlapPtr,
//		double *mOwned,
//		std::tr1::shared_ptr<double> volumeOverlapPtr,
//		double *dilatationOwned,
//		std::tr1::shared_ptr<double> fInternalOverlapPtr,
//		std::tr1::shared_ptr<int> localNeighborList,
//		int numOwnedPoints,
//		Pimp::PimpIsotropicHookeSpec matSpec
//);

void computeWeightedVolume
(
		const double* xOverlap,
		const double* volumeOverlap,
		double* mOwned,
		int myNumPoints,
		const int* localNeighborList
);

void computeDilatation
(
		const double* xOverlap,
		const double* yOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		double* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints
);

void computeInternalForceLinearElastic
(
		const double* xOverlapPtr,
		const double* yOverlapPtr,
		const double* mOwned,
		const double* volumeOverlapPtr,
		const double* dilatationOwned,
		const double* bondDamage,
		double* fInternalOverlapPtr,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS
);

void computeInternalForceIsotropicElasticPlastic
(
		const double* xOverlap,
		const double* yOverlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* deviatoricPlasticExtensionStateN,
		double* deviatoricPlasticExtensionStateNp1,
		double* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
		double HORIZON,
		double yieldStress
);

void updateGeometry
(
		const double* xOverlap,
		const double* uOverlap,
		const double* velocityOverlap,
		double* yOverlap,
		int overLapLength,
		double dt
);

/**
 * Computes norm of deviatoric force state at a particular point
 * @param numNeigh -- number of neighbors at point
 * @param theta    -- dilatation at point
 * @param neighPtr -- list of neighbors at point
 * @param bondDamage     -- damage parameter for each bond at point
 * @param X              -- original coordinates of point
 * @param Y              -- current coordinates of point
 * @param xOverlap       -- pointer to overlap vector of original coordinates; use this to get neighbor original coordinates
 * @param yOverlap       -- pointer to overlap vector of current coordinates; use this to get neighbor current coordinates
 * @param volumeOverlap  -- pointer to volume overlap vector; use this to get volume of neighboring points
 * @param alpha          -- material property (alpha = 15 mu / m
 * @param OMEGA          -- weight function at point
 */
double computeDeviatoricForceStateNorm
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

}

#endif /* PDMATERIALUTILITIES_H_ */
