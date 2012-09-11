/*! \file ordinary_utilities.h */

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


#ifndef ORDINARY_UTILITIES_H_
#define ORDINARY_UTILITIES_H_

#include <cstdlib>

class Bond_Volume_Calculator;

namespace MATERIAL_EVALUATION {

using std::size_t;

enum PURE_SHEAR { XY=0, YZ, ZX };

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
);

void computeWeightedVolume
(
		const double* xOverlap,
		const double* volumeOverlap,
		double *mOwned,
		int myNumPoints,
		const int* localNeighborList
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
		double weightedVolume
);

/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 * Y is the deformed configuration of 'X' but because
 * this is a probe, Y=X
 * This function returns the norm of the
 * deviatoric extension state squared.
 */
double compute_norm_2_deviatoric_extension
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		const double *Y,
		const double *yOverlap,
		const double *volumeOverlap,
		double weighted_volume
);

void computeShearCorrectionFactor
(
		int numOwnedPoints,
		int lengthYOverlap,
		const double *xOverlap,
		double *yOverlap_scratch_required_work_space,
		const double *volumeOverlap,
		const double *owned_weighted_volume,
		const int*  localNeighborList,
		double horizon,
		double *shearCorrectionFactorOwned
);

/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 */
double computeWeightedVolume
(
		const double *X,
		const double *xOverlap,
		const double* volumeOverlap,
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
		int numOwnedPoints
);

template<typename ScalarT>
void computeDilatationAD
(
		const double* xOverlap,
		const ScalarT* yOverlap,
		const double *mOwned,
		const double* volumeOverlap,
		const double* bondDamage,
		ScalarT* dilatationOwned,
		const int* localNeighborList,
		int numOwnedPoints
 );

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
		const int* localNeighborList
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
		const double *bondVolume,
		double weightedVolume
);


/**
 * Call this function on a single point 'X'
 * NOTE: neighPtr to should point to 'numNeigh' for 'X'
 * and thus describe the neighborhood list as usual
 * Y is the deformed configuration of 'X' but because
 * this is a probe, Y=X
 * This function returns the norm of the
 * deviatoric extension state squared.
 */
double compute_norm_2_deviatoric_extension
(
		const int *neighPtr,
		const double *X,
		const double *xOverlap,
		const double *Y,
		const double *yOverlap,
		const double *bondVolume,
		double weighted_volume
);

}


}


#endif /* ORDINARY_UTILITIES_H_ */
