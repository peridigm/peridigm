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
  const double *damage = bondDamage;
  const double *neighborVolume = neighborVolumePtr;
  const int *neighborlist = localNeighborList;

  const double *xNeighbor;
  const ScalarT *yNeighbor;
  ScalarT u[3], uNeighbor[3], dotProduct;
  double zeta[3], volNeighbor, normZeta, omega;
  int i, p, n, numNeighbors, neighborId;

  for(p=0; p<numOwnedPoints; p++, x+=3, y+=3, m++, theta++){
    *theta = 0.0;
    numNeighbors = *neighborlist;
    neighborlist++;
    for(n=0; n<numNeighbors; n++, neighborlist++, damage++, neighborVolume++){
      neighborId = *neighborlist;
      xNeighbor = &xOverlapPtr[3*neighborId];
      yNeighbor = &yOverlapPtr[3*neighborId];
      if(neighborVolumePtr != 0)
        volNeighbor = *neighborVolume;
      else
        volNeighbor = volumeOverlapPtr[neighborId];
      for(i=0 ; i<3 ; ++i){
        zeta[i] = xNeighbor[i] - x[i];
        u[i] = y[i] - x[i];
        uNeighbor[i] = yNeighbor[i] - xNeighbor[i];
      }
      normZeta = std::sqrt(zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2]);
      omega = influenceFunction(normZeta, horizon);
      dotProduct = zeta[0]*(uNeighbor[0]-u[0]) + zeta[1]*(uNeighbor[1]-u[1]) + zeta[2]*(uNeighbor[2]-u[2]);
      *theta += omega*(1.0 - *damage)*dotProduct*volNeighbor;
    }
    if(numNeighbors > 0){
      *theta *= 3.0/(*m);
    }
  }
}

template<typename ScalarT>
void computeInternalForceLinearLPS
(
 const double* xOverlapPtr,
 const ScalarT* yOverlapPtr,
 const double* volumeOverlapPtr,
 const double* weightedVolumePtr,
 const ScalarT* dilatationPtr,
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
 ScalarT* forceOverlapPtr,
 const int* localNeighborList,
 int numOwnedPoints,
 double bulkModulus,
 double shearModulus
)
{
  const double *x = xOverlapPtr;
  const ScalarT *y = yOverlapPtr;
  const double *m = weightedVolumePtr;
  const ScalarT *theta = dilatationPtr;
  const double *damage = bondDamage;
  const double *selfVolume = selfVolumePtr;
  const double *neighborVolume = neighborVolumePtr;
  ScalarT *force = forceOverlapPtr;
  const int *neighborlist = localNeighborList;

  const double *xNeighbor;
  const ScalarT *yNeighbor;
  ScalarT u[3], uNeighbor[3], temp1, matVec[3], fx, fy, fz;
  double zeta[3], volSelf, volNeighbor, normZeta, omega, temp2, dyadicProduct[3][3];
  int i, j, p, n, numNeighbors, neighborId;

  for(p=0; p<numOwnedPoints; p++, x+=3, y+=3, m++, theta++, force+=3){
    numNeighbors = *neighborlist;
    neighborlist++;
    for(n=0; n<numNeighbors; n++, neighborlist++, damage++, selfVolume++, neighborVolume++){
      neighborId = *neighborlist;
      xNeighbor = &xOverlapPtr[3*neighborId];
      yNeighbor = &yOverlapPtr[3*neighborId];
      if(neighborVolumePtr != 0){
        volSelf = *selfVolume;
        volNeighbor = *neighborVolume;
      }
      else{
        volSelf = volumeOverlapPtr[p];
        volNeighbor = volumeOverlapPtr[neighborId];
      }
      for(i=0 ; i<3 ; ++i){
        zeta[i] = xNeighbor[i] - x[i];
        u[i] = y[i] - x[i];
        uNeighbor[i] = yNeighbor[i] - xNeighbor[i];
      }
      normZeta = std::sqrt(zeta[0]*zeta[0] + zeta[1]*zeta[1] + zeta[2]*zeta[2]);
      omega = influenceFunction(normZeta, horizon);
      temp1 = (9.0*bulkModulus - 15.0*shearModulus)*omega*(*theta)/(3.0*(*m));
      temp2 = 15.0*shearModulus*omega/((*m)*normZeta*normZeta);
      for(i=0 ; i<3 ; i++){
        for(j=0 ; j<3 ; j++){
          dyadicProduct[i][j] = zeta[i]*zeta[j];
        }
      }
      for(i=0 ; i<3 ; ++i){
        matVec[i] = dyadicProduct[i][0]*(uNeighbor[0]-u[0]) + dyadicProduct[i][1]*(uNeighbor[1]-u[1]) + dyadicProduct[i][2]*(uNeighbor[2]-u[2]);
      }
      fx = (1.0 - *damage)*(temp1*zeta[0] + temp2*matVec[0]);
      fy = (1.0 - *damage)*(temp1*zeta[1] + temp2*matVec[1]);
      fz = (1.0 - *damage)*(temp1*zeta[2] + temp2*matVec[2]);
      *(force)   += fx*volNeighbor;
      *(force+1) += fy*volNeighbor;
      *(force+2) += fz*volNeighbor;
      forceOverlapPtr[3*neighborId]   -= fx*volSelf;
      forceOverlapPtr[3*neighborId+1] -= fy*volSelf;
      forceOverlapPtr[3*neighborId+2] -= fz*volSelf;
    } 
  }
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
 const double* volumeOverlapPtr,
 const double* weightedVolumePtr,
 const double* dilatationPtr,
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
 double* forceOverlapPtr,
 const int* localNeighborList,
 int numOwnedPoints,
 double bulkModulus,
 double shearModulus
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceLinearLPS<Sacado::Fad::DFad<double> >
(
 const double* xOverlapPtr,
 const Sacado::Fad::DFad<double>* yOverlapPtr,
 const double* volumeOverlapPtr,
 const double* weightedVolumePtr,
 const Sacado::Fad::DFad<double>* dilatationPtr,
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
 Sacado::Fad::DFad<double>* forceOverlapPtr,
 const int* localNeighborList,
 int numOwnedPoints,
 double bulkModulus,
 double shearModulus
);

}
