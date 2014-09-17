//! \file linear_lps_pv.cxx

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
#include "linear_lps_pv.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeDilatationLinearLPS
(
 const double* xOverlapPtr,
 const ScalarT* yOverlapPtr,
 const double* volumeOverlapPtr,
 const double* weightedVolumePtr,
 double horizon,
 const FunctionPointer influenceFunction,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 ScalarT* dilatationOwnedPtr,
 const int* localNeighborList,
 int numOwnedPoints
)
{
  const double *x = xOverlapPtr;
  const ScalarT *y = yOverlapPtr;
  const double *m = weightedVolumePtr;
  ScalarT *theta = dilatationOwnedPtr;
  const int *neighborlist = localNeighborList;

  const double *xNeighbor;
  const ScalarT *yNeighbor;
  ScalarT u[3], uNeighbor[3], dotProduct;
  double zeta[3], volNeighbor, normZeta, omega;
  int neighborId;

  for(int p=0; p<numOwnedPoints; p++, x+=3, y+=3, m++, theta++){
    *theta = 0.0;
    int numNeighbors = *neighborlist;
    neighborlist++;
    for(int n=0; n<numNeighbors; n++, neighborlist++){
      neighborId = *neighborlist;
      xNeighbor = &xOverlapPtr[3*neighborId];
      yNeighbor = &yOverlapPtr[3*neighborId];
      volNeighbor = volumeOverlapPtr[neighborId];
      for(int i=0 ; i<3 ; ++i){
        zeta[i] = xNeighbor[i] - x[i];
        u[i] = y[i] - x[i];
        uNeighbor[i] = yNeighbor[i] - xNeighbor[i];
      }
      normZeta = std::sqrt(zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2]);
      omega = influenceFunction(normZeta, horizon);
      dotProduct = zeta[0]*(uNeighbor[0]-u[0]) + zeta[1]*(uNeighbor[1]-u[1]) + zeta[2]*(uNeighbor[2]-u[2]);
      *theta += omega*dotProduct*volNeighbor;
    } 
    *theta /= *m;
  }
}

template<typename ScalarT>
void computeInternalForceLinearLPS
(
 const double* xOverlapPtr,
 const ScalarT* yOverlapPtr,
 const double* mOwned,
 const double* volumeOverlapPtr,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 ScalarT* forceOverlapPtr,
 const int*  localNeighborList,
 int numOwnedPoints,
 double BULK_MODULUS,
 double SHEAR_MODULUS,
 double horizon,
 ScalarT* partialStressOverlap
)
{
  // std::cout << "DEBUGGING " << numOwnedPoints << std::endl;

  // double K = BULK_MODULUS;
  // double MU = SHEAR_MODULUS;

  // const double *xOwned = xOverlapPtr;
  // const ScalarT *yOwned = yOverlapPtr;
  // const double *m = mOwned;
  // const double *v = volumeOverlapPtr;
  // ScalarT *fOwned = forceOverlapPtr;
  // ScalarT *psOwned = partialStressOverlap;

  // const int *neighPtr = localNeighborList;
  // double cellVolume, X_dx, X_dy, X_dz, zeta, omega, alpha;
  // ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e, c1;
  // for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, psOwned+=9, m++){

  //   int numNeigh = *neighPtr; neighPtr++;
  //   const double *X = xOwned;
  //   const ScalarT *Y = yOwned;
  //   alpha = 15.0*MU/(*m);
  //   double selfCellVolume = v[p];
  //   for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++){
  //     int localId = *neighPtr;
  //     cellVolume = v[localId];
  //     const double *XP = &xOverlapPtr[3*localId];
  //     const ScalarT *YP = &yOverlapPtr[3*localId];
  //     X_dx = XP[0]-X[0];
  //     X_dy = XP[1]-X[1];
  //     X_dz = XP[2]-X[2];
  //     zeta = sqrt(X_dx*X_dx+X_dy*X_dy+X_dz*X_dz);
  //     Y_dx = YP[0]-Y[0];
  //     Y_dy = YP[1]-Y[1];
  //     Y_dz = YP[2]-Y[2];
  //     dY = sqrt(Y_dx*Y_dx+Y_dy*Y_dy+Y_dz*Y_dz);
  //     e = dY - zeta;

  //     omega = scalarInfluenceFunction(zeta,horizon);

  //     // c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
  //     // t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
  //     // fx = t * Y_dx / dY;
  //     // fy = t * Y_dy / dY;
  //     // fz = t * Y_dz / dY;

  //     // *(fOwned+0) += fx*cellVolume;
  //     // *(fOwned+1) += fy*cellVolume;
  //     // *(fOwned+2) += fz*cellVolume;
  //     // forceOverlapPtr[3*localId+0] -= fx*selfCellVolume;
  //     // forceOverlapPtr[3*localId+1] -= fy*selfCellVolume;
  //     // forceOverlapPtr[3*localId+2] -= fz*selfCellVolume;

  //     // if(partialStressOverlap != 0){
  //     //   *(psOwned+0) += fx*X_dx*cellVolume;
  //     //   *(psOwned+1) += fx*X_dy*cellVolume;
  //     //   *(psOwned+2) += fx*X_dz*cellVolume;
  //     //   *(psOwned+3) += fy*X_dx*cellVolume;
  //     //   *(psOwned+4) += fy*X_dy*cellVolume;
  //     //   *(psOwned+5) += fy*X_dz*cellVolume;
  //     //   *(psOwned+6) += fz*X_dx*cellVolume;
  //     //   *(psOwned+7) += fz*X_dy*cellVolume;
  //     //   *(psOwned+8) += fz*X_dz*cellVolume;
  //     // }
  //   }

  // }
}

/** Explicit template instantiation for double. */
template void computeDilatationLinearLPS<double>
(
 const double* xOverlapPtr,
 const double* yOverlapPtr,
 const double* volumeOverlapPtr,
 const double* weightedVolumePtr,
 double horizon,
 const FunctionPointer influenceFunction,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 double* dilatationOwnedPtr,
 const int* localNeighborList,
 int numOwnedPoints
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeDilatationLinearLPS<Sacado::Fad::DFad<double> >
(
 const double* xOverlapPtr,
 const Sacado::Fad::DFad<double>* yOverlapPtr,
 const double* volumeOverlapPtr,
 const double* weightedVolumePtr,
 double horizon,
 const FunctionPointer influenceFunction,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 Sacado::Fad::DFad<double>* dilatationOwnedPtr,
 const int* localNeighborList,
 int numOwnedPoints
 );

/** Explicit template instantiation for double. */
template void computeInternalForceLinearLPS<double>
(
 const double* xOverlapPtr,
 const double* yOverlapPtr,
 const double* mOwned,
 const double* volumeOverlapPtr,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 double* forceOverlapPtr,
 const int*  localNeighborList,
 int numOwnedPoints,
 double BULK_MODULUS,
 double SHEAR_MODULUS,
 double horizon,
 double* partialStressOverlap = NULL
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceLinearLPS<Sacado::Fad::DFad<double> >
(
 const double* xOverlapPtr,
 const Sacado::Fad::DFad<double>* yOverlapPtr,
 const double* mOwned,
 const double* volumeOverlapPtr,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 Sacado::Fad::DFad<double>* forceOverlapPtr,
 const int*  localNeighborList,
 int numOwnedPoints,
 double BULK_MODULUS,
 double SHEAR_MODULUS,
 double horizon,
 Sacado::Fad::DFad<double>* partialStressOverlap = NULL
);

}
