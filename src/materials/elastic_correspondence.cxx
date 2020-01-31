//! \file elastic_correspondence.cxx

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

#include "elastic_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticCauchyStress
(
const ScalarT* deltaTemperatureN,
const ScalarT* deltaTemperatureNP1,
const ScalarT* unrotatedRateOfDeformation,
const ScalarT* unrotatedCauchyStressN,
ScalarT* unrotatedCauchyStressNP1,
const int numPoints,
const double bulkMod,
const double shearMod,
const double alpha,
const double dt
)
{
  const ScalarT* deltaTempN = deltaTemperatureN;
  const ScalarT* deltaTempNP1 = deltaTemperatureNP1;

  // Hooke's law
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* sigmaN = unrotatedCauchyStressN;
  ScalarT* sigmaNP1 = unrotatedCauchyStressNP1;

  ScalarT dilatationInc;
  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];

  for(int iID=0 ; iID<numPoints ; ++iID,
        rateOfDef+=9, sigmaN+=9, sigmaNP1+=9){

      //strainInc = dt * rateOfDef
      for (int i = 0; i < 9; i++) {
          strainInc[i] = *(rateOfDef+i)*dt;
          deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //thermal strains
      if (deltaTemperatureN && deltaTemperatureNP1) {
        double thermalStrainN = alpha*deltaTemperatureN[iID];
        double thermalStrainNP1 = alpha*deltaTemperatureNP1[iID];
        dilatationInc -= 3.0*(thermalStrainNP1 - thermalStrainN);
      }

      //update stress
      for (int i = 0; i < 9; i++) {
          *(sigmaNP1+i) = *(sigmaN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(sigmaNP1) += bulkMod*dilatationInc;
      *(sigmaNP1+4) += bulkMod*dilatationInc;
      *(sigmaNP1+8) += bulkMod*dilatationInc;

  }
}

// Explicit template instantiation for double
template void updateElasticCauchyStress<double>
(
const double* deltaTemperatureN,
const double* deltaTemperatureNP1,
const double* unrotatedRateOfDeformation,
const double* unrotatedCauchyStressN,
double* unrotatedCauchyStressNP1,
int numPoints,
double bulkMod,
double shearMod,
double alpha,
double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */

}
