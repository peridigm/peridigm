/*! \file ordinary_std_linear_visco_solid.cxx */

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
#include <iostream>
#include "viscoelastic.h"
using std::cout;
using std::endl;
namespace MATERIAL_EVALUATION {

void computeInternalForceViscoelasticStandardLinearSolid
  (
   double delta_t,
   const double *xOverlap,
   const double *yNOverlap,
   const double *yNP1Overlap,
   const double *mOwned,
   const double* volumeOverlap,
   const double* dilatationOwnedN,
   const double* dilatationOwnedNp1,
   const double* bondDamage,
   const double *edbN,
   double *edbNP1,
   double *fInternalOverlap,
   const int*  localNeighborList,
   int numOwnedPoints,
   double BULK_MODULUS,
   double SHEAR_MODULUS,
   double m_lambda_i,
   double m_tau_b_i
)
{

	double c1 = m_tau_b_i / delta_t;
	double decay = exp(-1.0/c1);
	double beta_i=1.-c1*(1.-decay);

	/*
	 * Compute processor local contribution to internal force
	 */
	double K = BULK_MODULUS;
	double MU = SHEAR_MODULUS;
	double OMEGA=1.0;

	const double *xOwned = xOverlap;
	const double *yNOwned = yNOverlap;
	const double *yNP1Owned = yNP1Overlap;
	const double *m = mOwned;
	const double *v = volumeOverlap;
	const double *thetaN = dilatationOwnedN;
	const double *thetaNp1 = dilatationOwnedNp1;
	double *fOwned = fInternalOverlap;

	const int *neighPtr = localNeighborList;
	double cellVolume, dx, dy, dz, zeta, dYN, dYNp1, t, ti, td, edN, edNp1, delta_ed;
	for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yNOwned +=3, yNP1Owned +=3, fOwned+=3, m++, thetaN++, thetaNp1++){

		int numNeigh = *neighPtr; neighPtr++;
		const double *X = xOwned;
		const double *YN = yNOwned;
		const double *YNP1 = yNP1Owned;
		double weightedVolume = *m;
		double dilatationN   = *thetaN;
		double dilatationNp1 = *thetaNp1;
		double alpha = 15.0*MU/weightedVolume;
		double selfCellVolume = v[p];
		double c = 3.0 * K * dilatationNp1 / weightedVolume;
		for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,edbN++,edbNP1++){
			int localId = *neighPtr;
			cellVolume = v[localId];
			const double *XP    = &xOverlap[3*localId];
			const double *YPN   = &yNOverlap[3*localId];
			const double *YPNP1 = &yNP1Overlap[3*localId];
			dx = XP[0]-X[0];
			dy = XP[1]-X[1];
			dz = XP[2]-X[2];
			zeta = sqrt(dx*dx+dy*dy+dz*dz);

			/*
			 * JAM:damage state
			 * Some additional analysis (pencil and paper) may be required with
			 * regard to handling damage.  The question is where to
			 * apply the damage.  The approach taken here is to
			 * apply damage to the incoming (computed) deviatoric extension state
			 * and evolve the back extension state with the damaged
			 * deviatoric extension state.
			 *
			 * Also note that this assumes the
			 * damage on step N is the same as step NP1;
			 * ERROR: this needs to be fixed.
			 */
			double damageN = (1.0-*bondDamage);
			double damageNp1 = (1.0-*bondDamage);

			/*
			 * volumetric scalar state
			 */
			double eiN   = dilatationN * zeta / 3.0;
			double eiNp1 = dilatationNp1 * zeta / 3.0;

			/*
			 * COMPUTE edN
			 */
			dx = YPN[0]-YN[0];
			dy = YPN[1]-YN[1];
			dz = YPN[2]-YN[2];
			dYN = sqrt(dx*dx+dy*dy+dz*dz);
			/*
			 */
			edN = damageN * (dYN - zeta) - eiN;

			/*
			 * COMPUTE edNp1
			 */
			dx = YPNP1[0]-YNP1[0];
			dy = YPNP1[1]-YNP1[1];
			dz = YPNP1[2]-YNP1[2];
			dYNp1 = sqrt(dx*dx+dy*dy+dz*dz);
			edNp1 = damageNp1 * (dYNp1 - zeta) - eiNp1;

			/*
			 * Increment to deviatoric extension state
			 */
			delta_ed = edNp1-edN;
			/*
			 * Integrate back extension state forward in time
			 */
			*edbNP1 = edN * (1-decay) + (*edbN)*decay  + beta_i * delta_ed;

			/*
			 * Compute deviatoric force state
			 */
			td = (1.0-m_lambda_i) * alpha * OMEGA * edNp1 + m_lambda_i * alpha * OMEGA * ( edNp1 - *edbNP1 );

			/*
			 * Compute volumetric force state
			 */
			ti = c * OMEGA * zeta;

			/*
			 * Note that damage has already been applied once to 'td' (through ed) above.
			 */
			t = damageNp1 * (ti + td);
			double fx = t * dx / dYNp1;
			double fy = t * dy / dYNp1;
			double fz = t * dz / dYNp1;

			*(fOwned+0) += fx*cellVolume;
			*(fOwned+1) += fy*cellVolume;
			*(fOwned+2) += fz*cellVolume;
			fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
			fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
			fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
		}
	}

}

}

