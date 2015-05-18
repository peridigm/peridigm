//! \file elastic_pv.h

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
#ifndef ELASTICPV_H
#define ELASTICPV_H

#include "Peridigm_InfluenceFunction.hpp"

namespace MATERIAL_EVALUATION {

typedef PeridigmNS::InfluenceFunction::functionPointer FunctionPointer;

double computeWeightedVolumePV
(
 const double *X,
 const double *xOverlap,
 const double* volumeOverlap,
 const double* selfVolumePtr,
 const double* selfCentroidXPtr,
 const double* selfCentroidYPtr,
 const double* selfCentroidZPtr,
 const double* neighborVolumePtr,
 const double* neighborCentroidXPtr,
 const double* neighborCentroidYPtr,
 const double* neighborCentroidZPtr,
 const double* influenceFunctionValuesPtr,
 const int* localNeighborList,
 double horizon,
 const FunctionPointer omega=PeridigmNS::InfluenceFunction::self().getInfluenceFunction()
);

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
 const double* influenceFunctionValuesPtr,
 double *mOwned,
 int myNumPoints,
 const int* localNeighborList,
 double horizon,
 const FunctionPointer omega = PeridigmNS::InfluenceFunction::self().getInfluenceFunction()
);

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
 const FunctionPointer omega = PeridigmNS::InfluenceFunction::self().getInfluenceFunction(),
 double thermalExpansionCoefficient = 0,
 const double* deltaTemperature = 0
);

//! Computes contributions to the internal force resulting from owned points.
template<typename ScalarT>
void computeInternalForceLinearElasticPV
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
 const ScalarT* dilatationOwned,
 const double* bondDamage,
 ScalarT* fInternalOverlapPtr,
 const int*  localNeighborList,
 int numOwnedPoints,
 double BULK_MODULUS,
 double SHEAR_MODULUS,
 double horizon,
 double thermalExpansionCoefficient = 0,
 const double* deltaTemperature = 0
);

}

#endif // ELASTICPV_H
