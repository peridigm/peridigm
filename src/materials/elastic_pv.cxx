//! \file elastic_pv.cxx

//@HEADER
// ************************************************************************
//
//                             Peridigm
//                 Copyright (2014) Sandia Corporation
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
#include "elastic.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

double computeWeightedVolumePV
(
 const double* X,
 const double* xOverlap,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const int* localNeighborList,
 double horizon,
 const FunctionPointer omega
){
  double neighborVolume;
  double neighborX, neighborY, neighborZ;
  double selfX, selfY, selfZ;
  double m = 0.0;
  const int *neighPtr = localNeighborList;
  int numNeigh = *neighPtr; neighPtr++;
  for(int n=0 ; n<numNeigh ; n++,neighPtr++){
    int localId = *neighPtr;

    if(neighborVolumePtr != 0){
      neighborVolume = neighborVolumePtr[n];
    }
    else{
      neighborVolume = volumeOverlap[localId];
    }

    if(selfCentroidXPtr != 0){
      selfX = X[0];//selfCentroidXPtr[n];
      selfY = X[1];//selfCentroidYPtr[n];
      selfZ = X[2];//selfCentroidZPtr[n];
      neighborX = neighborCentroidXPtr[n];
      neighborY = neighborCentroidYPtr[n];
      neighborZ = neighborCentroidZPtr[n];
    }
    else{
      selfX = X[0];
      selfY = X[1];
      selfZ = X[2];
      neighborX = xOverlap[3*localId];
      neighborY = xOverlap[3*localId+1];
      neighborZ = xOverlap[3*localId+2];
    }

    double dx = neighborX - selfX;
    double dy = neighborY - selfY;
    double dz = neighborZ - selfZ;
    double zetaSquared = dx*dx+dy*dy+dz*dz;
    double d = sqrt(zetaSquared);
    m += omega(d,horizon)*zetaSquared*neighborVolume;
  }
  return m;
}

void computeWeightedVolumePV
(
 const double* xOverlap,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 double *mOwned,
 int myNumPoints,
 const int* localNeighborList,
 double horizon,
 const FunctionPointer omega
){
  double *m = mOwned;
  const double *xOwned = xOverlap;
  const int *neighPtr = localNeighborList;
  const double *selfVolume = selfVolumePtr;
  const double *selfCentroidX = selfCentroidXPtr;
  const double *selfCentroidY = selfCentroidYPtr;
  const double *selfCentroidZ = selfCentroidZPtr;
  const double *neighborVolume = neighborVolumePtr;
  const double *neighborCentroidX = neighborCentroidXPtr;
  const double *neighborCentroidY = neighborCentroidYPtr;
  const double *neighborCentroidZ = neighborCentroidZPtr;

  bool usePartialVolume = false;
  if(selfVolume != 0 && neighborVolume != 0)
    usePartialVolume = true;
  bool usePartialCentroid = false;
  if(selfCentroidX != 0 && selfCentroidY != 0 && selfCentroidZ != 0 && neighborCentroidX != 0 && neighborCentroidY != 0 && neighborCentroidZ != 0)
    usePartialCentroid = true;

  for(int p=0;p<myNumPoints;p++, xOwned+=3, m++){
    int numNeigh = *neighPtr;
    const double *X = xOwned;
    *m=MATERIAL_EVALUATION::computeWeightedVolumePV(X,
						    xOverlap,
						    volumeOverlap,
						    selfVolume,
						    selfCentroidX,
						    selfCentroidY,
						    selfCentroidZ,
						    neighborVolume,
						    neighborCentroidX,
						    neighborCentroidY,
						    neighborCentroidZ,
						    neighPtr,
						    horizon,
						    omega);
    neighPtr+=(numNeigh+1);
    if(usePartialVolume){
      selfVolume+=numNeigh;
      neighborVolume+=numNeigh;
    }
    if(usePartialCentroid){
      selfCentroidX+=numNeigh;
      selfCentroidY+=numNeigh;
      selfCentroidZ+=numNeigh;
      neighborCentroidX+=numNeigh;
      neighborCentroidY+=numNeigh;
      neighborCentroidZ+=numNeigh;
    }
  }
}

template<typename ScalarT>
void computeDilatationPV
(
 const double* xOverlap,
 const ScalarT* yOverlap,
 const double *mOwned,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 ScalarT* dilatationOwned,
 const int* localNeighborList,
 int numOwnedPoints,
 double horizon,
 const FunctionPointer omega,
 double thermalExpansionCoefficient,
 const double* deltaTemperature
)
{
  const double *xOwned = xOverlap;
  const ScalarT *yOwned = yOverlap;
  const double *deltaT = deltaTemperature;
  const double *m = mOwned;
  const double *v = volumeOverlap;
  const double *neighborVolume = neighborVolumePtr;
  ScalarT *theta = dilatationOwned;
  double cellVolume;
  const int *neighPtr = localNeighborList;
  for(int p=0; p<numOwnedPoints;p++, xOwned+=3, yOwned+=3, deltaT++, m++, theta++){
    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const ScalarT *Y = yOwned;
    *theta = ScalarT(0.0);
    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,neighborVolume++){
      int localId = *neighPtr;
      if(neighborVolumePtr != 0)
	cellVolume = *neighborVolume;
      else
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
      ScalarT e = sqrt(dY);
      e -= d;
      if(deltaTemperature)
	e -= thermalExpansionCoefficient*(*deltaT)*d;
      double omegaVal = omega(d,horizon);
      *theta += 3.0*omegaVal*(1.0-*bondDamage)*d*e*cellVolume/(*m);
    }
  }
}

template<typename ScalarT>
void computeInternalForceLinearElasticPV
(
 const double* xOverlap,
 const ScalarT* yOverlap,
 const double* mOwned,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const ScalarT* dilatationOwned,
 const double* bondDamage,
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
  const double *deltaT = deltaTemperature;
  const double *m = mOwned;
  const double *selfVolume = selfVolumePtr;
  const double *neighborVolume = neighborVolumePtr;
  const ScalarT *theta = dilatationOwned;
  ScalarT *fOwned = fInternalOverlap;

  const int *neighPtr = localNeighborList;
  double alpha, X_dx, X_dy, X_dz, zeta, omega, selfCellVolume, neighborCellVolume;
  ScalarT Y_dx, Y_dy, Y_dz, dY, t, fx, fy, fz, e, c1;
  for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, deltaT++, m++, theta++){

    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const ScalarT *Y = yOwned;
    alpha = 15.0*MU/(*m);

    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++,selfVolume++,neighborVolume++){
      int localId = *neighPtr;
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
      c1 = omega*(*theta)*(3.0*K/(*m)-alpha/3.0);
      t = (1.0-*bondDamage)*(c1 * zeta + (1.0-*bondDamage) * omega * alpha * e);
      fx = t * Y_dx / dY;
      fy = t * Y_dy / dY;
      fz = t * Y_dz / dY;

      if(selfVolumePtr != 0 && neighborVolumePtr != 0){
	selfCellVolume = *selfVolume;
	neighborCellVolume = *neighborVolume;
      }
      else{
	selfCellVolume = volumeOverlap[p];
	neighborCellVolume = volumeOverlap[localId];
      }

      *(fOwned+0) += fx*neighborCellVolume;
      *(fOwned+1) += fy*neighborCellVolume;
      *(fOwned+2) += fz*neighborCellVolume;
      fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
      fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
      fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
    }
  }
}

/** Explicit template instantiation for double. */
template
void computeDilatationPV<double>
(
 const double* xOverlap,
 const double* yOverlap,
 const double *mOwned,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 double* dilatationOwned,
 const int* localNeighborList,
 int numOwnedPoints,
 double horizon,
 const FunctionPointer omega,
 double thermalExpansionCoefficient,
 const double* deltaTemperature
 );

template void computeInternalForceLinearElasticPV<double>
(
 const double* xOverlap,
 const double* yOverlap,
 const double* mOwned,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* dilatationOwned,
 const double* bondDamage,
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
template
void computeDilatationPV<Sacado::Fad::DFad<double> >
(
 const double* xOverlap,
 const Sacado::Fad::DFad<double>* yOverlap,
 const double *mOwned,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* bondDamage,
 Sacado::Fad::DFad<double>* dilatationOwned,
 const int* localNeighborList,
 int numOwnedPoints,
 double horizon,
 const FunctionPointer omega,
 double thermalExpansionCoefficient,
 const double* deltaTemperature
 );

template void computeInternalForceLinearElasticPV<Sacado::Fad::DFad<double> >
(
 const double* xOverlap,
 const Sacado::Fad::DFad<double>* yOverlap,
 const double* mOwned,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const Sacado::Fad::DFad<double>* dilatationOwned,
 const double* bondDamage,
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
