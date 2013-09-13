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
const ScalarT* cauchyStressNXX, 
const ScalarT* cauchyStressNXY, 
const ScalarT* cauchyStressNXZ, 
const ScalarT* cauchyStressNYX, 
const ScalarT* cauchyStressNYY, 
const ScalarT* cauchyStressNYZ, 
const ScalarT* cauchyStressNZX, 
const ScalarT* cauchyStressNZY, 
const ScalarT* cauchyStressNZZ, 
ScalarT* cauchyStressNP1XX, 
ScalarT* cauchyStressNP1XY, 
ScalarT* cauchyStressNP1XZ, 
ScalarT* cauchyStressNP1YX, 
ScalarT* cauchyStressNP1YY, 
ScalarT* cauchyStressNP1YZ, 
ScalarT* cauchyStressNP1ZX, 
ScalarT* cauchyStressNP1ZY, 
ScalarT* cauchyStressNP1ZZ, 
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
  const ScalarT* sigmaNXX = cauchyStressNXX;
  const ScalarT* sigmaNXY = cauchyStressNXY;
  const ScalarT* sigmaNXZ = cauchyStressNXZ;
  const ScalarT* sigmaNYX = cauchyStressNYX;
  const ScalarT* sigmaNYY = cauchyStressNYY;
  const ScalarT* sigmaNYZ = cauchyStressNYZ;
  const ScalarT* sigmaNZX = cauchyStressNZX;
  const ScalarT* sigmaNZY = cauchyStressNZY;
  const ScalarT* sigmaNZZ = cauchyStressNZZ;
  ScalarT* sigmaNP1XX = cauchyStressNP1XX;
  ScalarT* sigmaNP1XY = cauchyStressNP1XY;
  ScalarT* sigmaNP1XZ = cauchyStressNP1XZ;
  ScalarT* sigmaNP1YX = cauchyStressNP1YX;
  ScalarT* sigmaNP1YY = cauchyStressNP1YY;
  ScalarT* sigmaNP1YZ = cauchyStressNP1YZ;
  ScalarT* sigmaNP1ZX = cauchyStressNP1ZX;
  ScalarT* sigmaNP1ZY = cauchyStressNP1ZY;
  ScalarT* sigmaNP1ZZ = cauchyStressNP1ZZ;

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
        ++sigmaNXX, ++sigmaNXY, ++sigmaNXZ,
        ++sigmaNYX, ++sigmaNYY, ++sigmaNYZ,
        ++sigmaNZX, ++sigmaNZY, ++sigmaNZZ,
        ++sigmaNP1XX, ++sigmaNP1XY, ++sigmaNP1XZ,
        ++sigmaNP1YX, ++sigmaNP1YY, ++sigmaNP1YZ,
        ++sigmaNP1ZX, ++sigmaNP1ZY, ++sigmaNP1ZZ){

      strainIncXX = dt * (*rateOfDefXX);
      strainIncXY = dt * (*rateOfDefXY);
      strainIncXZ = dt * (*rateOfDefXZ);
      strainIncYX = dt * (*rateOfDefYX);
      strainIncYY = dt * (*rateOfDefYY);
      strainIncYZ = dt * (*rateOfDefYZ);
      strainIncZX = dt * (*rateOfDefZX);
      strainIncZY = dt * (*rateOfDefZY);
      strainIncZZ = dt * (*rateOfDefZZ);

      dilatationInc = strainIncXX + strainIncYY + strainIncZZ;

      deviatoricStrainIncXX = strainIncXX - dilatationInc/3.0;
      deviatoricStrainIncXY = strainIncXY;
      deviatoricStrainIncXZ = strainIncXZ;
      deviatoricStrainIncYX = strainIncYX;
      deviatoricStrainIncYY = strainIncYY - dilatationInc/3.0;
      deviatoricStrainIncYZ = strainIncYZ;
      deviatoricStrainIncZX = strainIncZX;
      deviatoricStrainIncZY = strainIncZY;
      deviatoricStrainIncZZ = strainIncZZ - dilatationInc/3.0;

      *sigmaNP1XX = (*sigmaNXX) + 2.0*shearMod*deviatoricStrainIncXX + bulkMod*dilatationInc;
      *sigmaNP1XY = (*sigmaNXY) + 2.0*shearMod*deviatoricStrainIncXY;
      *sigmaNP1XZ = (*sigmaNXZ) + 2.0*shearMod*deviatoricStrainIncXZ;
      *sigmaNP1YX = (*sigmaNYX) + 2.0*shearMod*deviatoricStrainIncYX;
      *sigmaNP1YY = (*sigmaNYY) + 2.0*shearMod*deviatoricStrainIncYY + bulkMod*dilatationInc;
      *sigmaNP1YZ = (*sigmaNYZ) + 2.0*shearMod*deviatoricStrainIncYZ;
      *sigmaNP1ZX = (*sigmaNZX) + 2.0*shearMod*deviatoricStrainIncZX;
      *sigmaNP1ZY = (*sigmaNZY) + 2.0*shearMod*deviatoricStrainIncZY;
      *sigmaNP1ZZ = (*sigmaNZZ) + 2.0*shearMod*deviatoricStrainIncZZ + bulkMod*dilatationInc;


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
const double* cauchyStressNXX, 
const double* cauchyStressNXY, 
const double* cauchyStressNXZ, 
const double* cauchyStressNYX, 
const double* cauchyStressNYY, 
const double* cauchyStressNYZ, 
const double* cauchyStressNZX, 
const double* cauchyStressNZY, 
const double* cauchyStressNZZ, 
double* cauchyStressNP1XX, 
double* cauchyStressNP1XY, 
double* cauchyStressNP1XZ, 
double* cauchyStressNP1YX, 
double* cauchyStressNP1YY, 
double* cauchyStressNP1YZ, 
double* cauchyStressNP1ZX, 
double* cauchyStressNP1ZY, 
double* cauchyStressNP1ZZ, 
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
const Sacado::Fad::DFad<double>* cauchyStressNXX, 
const Sacado::Fad::DFad<double>* cauchyStressNXY, 
const Sacado::Fad::DFad<double>* cauchyStressNXZ, 
const Sacado::Fad::DFad<double>* cauchyStressNYX, 
const Sacado::Fad::DFad<double>* cauchyStressNYY, 
const Sacado::Fad::DFad<double>* cauchyStressNYZ, 
const Sacado::Fad::DFad<double>* cauchyStressNZX, 
const Sacado::Fad::DFad<double>* cauchyStressNZY, 
const Sacado::Fad::DFad<double>* cauchyStressNZZ, 
Sacado::Fad::DFad<double>* cauchyStressNP1XX, 
Sacado::Fad::DFad<double>* cauchyStressNP1XY, 
Sacado::Fad::DFad<double>* cauchyStressNP1XZ, 
Sacado::Fad::DFad<double>* cauchyStressNP1YX, 
Sacado::Fad::DFad<double>* cauchyStressNP1YY, 
Sacado::Fad::DFad<double>* cauchyStressNP1YZ, 
Sacado::Fad::DFad<double>* cauchyStressNP1ZX, 
Sacado::Fad::DFad<double>* cauchyStressNP1ZY, 
Sacado::Fad::DFad<double>* cauchyStressNP1ZZ, 
int numPoints, 
double bulkMod,
double shearMod,
double dt
);

}
