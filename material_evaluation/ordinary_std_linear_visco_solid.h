/*! \file ordinary_std_linear_visco_solid.h */

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
#ifndef ORDINARY_STD_LINEAR_VISCO_SOLID_H_
#define ORDINARY_STD_LINEAR_VISCO_SOLID_H_


namespace MATERIAL_EVALUATION {

/**
 * Internal force calculator for viscoelastic standard linear solid.
 * Integrates the deviatoric back strain forward in time.
 * Output:
 *   * force
 *   * edbPN1 -- deviatoric back strain at end of step
 */
void computeInternalForceViscoelasticStandardLinearSolid
  (double delta_t,
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
   double m_bulkModulus,
   double m_shearModulus,
   double m_tau,
   double m_tau_b
);

}

#endif /* ORDINARY_STD_LINEAR_VISCO_SOLID_H_ */
