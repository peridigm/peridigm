//! \file elastic_bond_based.cxx

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
#include "elastic_bond_based.h"
#include "material_utilities.h"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
void computeInternalForceElasticBondBased
(
    const double* xOverlap,
    const ScalarT* yOverlap,
    const double* volumeOverlap,
    const double* bondDamage,
    ScalarT* fInternalOverlap,
    const int* localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double horizon,
    const bool isPlaneStrainStress,
    const double height
)
{
  double volume, neighborVolume, X[3], neighborX[3], initialBondLength, damageOnBond;
  ScalarT Y[3], neighborY[3], currentBondLength, stretch, t, fx, fy, fz;
  int neighborhoodIndex(0), bondDamageIndex(0), neighborId;

  const double pi = PeridigmNS::value_of_pi();
  double constant = 0.0;
  if(isPlaneStrainStress)
    constant = 12.0*BULK_MODULUS/(pi*height*horizon*horizon*horizon); // plain strain and plain stress
  else
    constant = 18.0*BULK_MODULUS/(pi*horizon*horizon*horizon*horizon); // 3d

  for(int p=0 ; p<numOwnedPoints ; p++){

    X[0] = xOverlap[p*3];
    X[1] = xOverlap[p*3+1];
    X[2] = xOverlap[p*3+2];
    Y[0] = yOverlap[p*3];
    Y[1] = yOverlap[p*3+1];
    Y[2] = yOverlap[p*3+2];
    volume = volumeOverlap[p];

    int numNeighbors = localNeighborList[neighborhoodIndex++];
    for(int n=0; n<numNeighbors; n++){

      neighborId = localNeighborList[neighborhoodIndex++];
      neighborX[0] = xOverlap[neighborId*3];
      neighborX[1] = xOverlap[neighborId*3+1];
      neighborX[2] = xOverlap[neighborId*3+2];
      neighborY[0] = yOverlap[neighborId*3];
      neighborY[1] = yOverlap[neighborId*3+1];
      neighborY[2] = yOverlap[neighborId*3+2];
      neighborVolume = volumeOverlap[neighborId];
      
      initialBondLength = std::sqrt( (neighborX[0]-X[0])*(neighborX[0]-X[0]) + (neighborX[1]-X[1])*(neighborX[1]-X[1]) + (neighborX[2]-X[2])*(neighborX[2]-X[2]) );
      currentBondLength = std::sqrt( (neighborY[0]-Y[0])*(neighborY[0]-Y[0]) + (neighborY[1]-Y[1])*(neighborY[1]-Y[1]) + (neighborY[2]-Y[2])*(neighborY[2]-Y[2]) );
      stretch = (currentBondLength - initialBondLength)/initialBondLength;

      damageOnBond = bondDamage[bondDamageIndex++];

      t = 0.5*(1.0 - damageOnBond)*stretch*constant;

      fx = t * (neighborY[0] - Y[0]) / currentBondLength;
      fy = t * (neighborY[1] - Y[1]) / currentBondLength;
      fz = t * (neighborY[2] - Y[2]) / currentBondLength;

      fInternalOverlap[3*p+0] += fx*neighborVolume;
      fInternalOverlap[3*p+1] += fy*neighborVolume;
      fInternalOverlap[3*p+2] += fz*neighborVolume;
      fInternalOverlap[3*neighborId+0] -= fx*volume;
      fInternalOverlap[3*neighborId+1] -= fy*volume;
      fInternalOverlap[3*neighborId+2] -= fz*volume;
    }
  }
}

/** Explicit template instantiation for double. */
template void computeInternalForceElasticBondBased<double>
(
    const double* xOverlap,
    const double* yOverlap,
    const double* volumeOverlap,
    const double* bondDamage,
    double* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double horizon,
    const bool isPlaneStrainStress,
    const double height
 );

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceElasticBondBased<Sacado::Fad::DFad<double> >
(
    const double* xOverlap,
    const Sacado::Fad::DFad<double>* yOverlap,
    const double* volumeOverlap,
    const double* bondDamage,
    Sacado::Fad::DFad<double>* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double horizon,
    const bool isPlaneStrainStress,
    const double height
);

}
