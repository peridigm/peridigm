/*! \file ordinary_iso_hard_plastic.h */

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
#ifndef ORDINARY_ISO_HARD_PLASTIC_H_
#define ORDINARY_ISO_HARD_PLASTIC_H_

namespace MATERIAL_EVALUATION {

double updateLambdaNP1
(
		double tdNorm,
		const double lambdaN,
        double pointWiseYieldValue,
	    double alpha,
        double HARD_MODULUS
);

//template <typename ScalarT>
//void updateLambdaNP1AD
//(
		//double tdNorm,
		//const double lambdaN,
        //double pointWiseYieldValue,
		//double alpha,
        //double HARD_MODULUS,
        //double dt
//);

void computeInternalForceIsotropicHardeningPlastic
(
		const double* xOverlap,
		const double *yNP1Overlap,
		const double* mOwned,
		const double* volumeOverlap,
		const double* dilatationOwned,
		const double* bondDamage,
		const double* dsfOwned_,
		const double* deviatoricPlasticExtensionStateN,
		double* deviatoricPlasticExtensionStateNp1,
		const double* lambdaN_,
		double* lambdaNP1_,
		double* fInternalOverlap,
		const int*  localNeighborList,
		int numOwnedPoints,
		double BULK_MODULUS,
		double SHEAR_MODULUS,
		double HORIZON,
		double yieldStress,
		double HARD_MODULUS
);

//template<typename ScalarT>
//void computeInternalForceIsotropicElasticPlasticAD
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
		//const int* localNeighborList,
		//int numOwnedPoints,
		//double BULK_MODULUS,
		//double SHEAR_MODULUS,
		//double HORIZON,
		//double yieldStress,
		//double HARD_MODULUS,
		//double dt
//);

}

#endif /* ORDINARY_ISO_HARD_PLASTIC_H_ */
