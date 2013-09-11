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
#include "material_utilities.h"
#include <Sacado.hpp>

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticCauchyStress
(
const ScalarT* unrotatedRateOfDeformationXX, 
const ScalarT* unrotatedRateOfDeformationXY, 
const ScalarT* unrotatedRateOfDeformationXZ, 
const ScalarT* unrotatedRateOfDeformationYX, 
const ScalarT* unrotatedRateOfDeformationYY, 
const ScalarT* unrotatedRateOfDeformationYZ, 
const ScalarT* unrotatedRateOfDeformationZX, 
const ScalarT* unrotatedRateOfDeformationZY, 
const ScalarT* unrotatedRateOfDeformationZZ, 
ScalarT* cauchyStressXX, 
ScalarT* cauchyStressXY, 
ScalarT* cauchyStressXZ, 
ScalarT* cauchyStressYX, 
ScalarT* cauchyStressYY, 
ScalarT* cauchyStressYZ, 
ScalarT* cauchyStressZX, 
ScalarT* cauchyStressZY, 
ScalarT* cauchyStressZZ, 
const int numPoints, 
const double bulkMod,
const double shearMod,
const double dt
)
{
  // Hooke's law
  const ScalarT* rateOfDefXX = unrotatedRateOfDeformationXX;
  const ScalarT* rateOfDefXY = unrotatedRateOfDeformationXY;
  const ScalarT* rateOfDefXZ = unrotatedRateOfDeformationXZ;
  const ScalarT* rateOfDefYX = unrotatedRateOfDeformationYX;
  const ScalarT* rateOfDefYY = unrotatedRateOfDeformationYY;
  const ScalarT* rateOfDefYZ = unrotatedRateOfDeformationYZ;
  const ScalarT* rateOfDefZX = unrotatedRateOfDeformationZX;
  const ScalarT* rateOfDefZY = unrotatedRateOfDeformationZY;
  const ScalarT* rateOfDefZZ = unrotatedRateOfDeformationZZ;
  ScalarT* sigmaXX = cauchyStressXX;
  ScalarT* sigmaXY = cauchyStressXY;
  ScalarT* sigmaXZ = cauchyStressXZ;
  ScalarT* sigmaYX = cauchyStressYX;
  ScalarT* sigmaYY = cauchyStressYY;
  ScalarT* sigmaYZ = cauchyStressYZ;
  ScalarT* sigmaZX = cauchyStressZX;
  ScalarT* sigmaZY = cauchyStressZY;
  ScalarT* sigmaZZ = cauchyStressZZ;

  ScalarT strainIncXX, strainIncXY, strainIncXZ;
  ScalarT strainIncYX, strainIncYY, strainIncYZ;
  ScalarT strainIncZX, strainIncZY, strainIncZZ;

  ScalarT deviatoricStrainIncXX, deviatoricStrainIncXY, deviatoricStrainIncXZ;
  ScalarT deviatoricStrainIncYX, deviatoricStrainIncYY, deviatoricStrainIncYZ;
  ScalarT deviatoricStrainIncZX, deviatoricStrainIncZY, deviatoricStrainIncZZ;

  ScalarT dilatationInc;

  for(int iID=0 ; iID<numPoints ; ++iID, 
        ++rateOfDefXX, ++rateOfDefXY, ++rateOfDefXZ,
        ++rateOfDefYX, ++rateOfDefYY, ++rateOfDefYZ,
        ++rateOfDefZX, ++rateOfDefZY, ++rateOfDefZZ,
        ++sigmaXX, ++sigmaXY, ++sigmaXZ,
        ++sigmaYX, ++sigmaYY, ++sigmaYZ,
        ++sigmaZX, ++sigmaZY, ++sigmaZZ){

      strainIncXX = dt * (*rateOfDefXX);
      strainIncXY = dt * (*rateOfDefXY);
      strainIncXX = dt * (*rateOfDefXZ);
      strainIncYX = dt * (*rateOfDefYX);
      strainIncYY = dt * (*rateOfDefYY);
      strainIncYZ = dt * (*rateOfDefYZ);
      strainIncZX = dt * (*rateOfDefZX);
      strainIncZY = dt * (*rateOfDefZY);
      strainIncZZ = dt * (*rateOfDefZZ);

      dilatationInc = strainIncXX + strainIncYY + strainIncZZ;

      deviatoricStrainIncXX = strainIncXX - dilatationInc/3.0;
      deviatoricStrainIncXY = strainIncXY;
      deviatoricStrainIncXX = strainIncXZ;
      deviatoricStrainIncYX = strainIncYX;
      deviatoricStrainIncYY = strainIncYY - dilatationInc/3.0;
      deviatoricStrainIncYZ = strainIncYZ;
      deviatoricStrainIncZX = strainIncZX;
      deviatoricStrainIncZY = strainIncZY;
      deviatoricStrainIncZZ = strainIncZZ - dilatationInc/3.0;

      *sigmaXX += 2.0*shearMod*deviatoricStrainIncXX + bulkMod*dilatationInc;
      *sigmaXY += 2.0*shearMod*deviatoricStrainIncXY;
      *sigmaXZ += 2.0*shearMod*deviatoricStrainIncXZ;
      *sigmaYX += 2.0*shearMod*deviatoricStrainIncYX;
      *sigmaYY += 2.0*shearMod*deviatoricStrainIncYY + bulkMod*dilatationInc;
      *sigmaYZ += 2.0*shearMod*deviatoricStrainIncYZ;
      *sigmaZX += 2.0*shearMod*deviatoricStrainIncZX;
      *sigmaZY += 2.0*shearMod*deviatoricStrainIncZY;
      *sigmaZZ += 2.0*shearMod*deviatoricStrainIncZZ + bulkMod*dilatationInc;

  }
}

// Explicit template instantiation for double
template void updateElasticCauchyStress<double>
(
const double* unrotatedRateOfDeformationXX, 
const double* unrotatedRateOfDeformationXY, 
const double* unrotatedRateOfDeformationXZ, 
const double* unrotatedRateOfDeformationYX, 
const double* unrotatedRateOfDeformationYY, 
const double* unrotatedRateOfDeformationYZ, 
const double* unrotatedRateOfDeformationZX, 
const double* unrotatedRateOfDeformationZY, 
const double* unrotatedRateOfDeformationZZ, 
double* cauchyStressXX, 
double* cauchyStressXY, 
double* cauchyStressXZ, 
double* cauchyStressYX, 
double* cauchyStressYY, 
double* cauchyStressYZ, 
double* cauchyStressZX, 
double* cauchyStressZY, 
double* cauchyStressZZ, 
int numPoints, 
double bulkMod,
double shearMod,
double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void updateElasticCauchyStress<Sacado::Fad::DFad<double> >
(
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXX, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXY, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXZ, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYX, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYY, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYZ, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZX, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZY, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZZ, 
Sacado::Fad::DFad<double>* cauchyStressXX, 
Sacado::Fad::DFad<double>* cauchyStressXY, 
Sacado::Fad::DFad<double>* cauchyStressXZ, 
Sacado::Fad::DFad<double>* cauchyStressYX, 
Sacado::Fad::DFad<double>* cauchyStressYY, 
Sacado::Fad::DFad<double>* cauchyStressYZ, 
Sacado::Fad::DFad<double>* cauchyStressZX, 
Sacado::Fad::DFad<double>* cauchyStressZY, 
Sacado::Fad::DFad<double>* cauchyStressZZ, 
int numPoints, 
double bulkMod,
double shearMod,
double dt
);

}
