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
const ScalarT* unrotatedCauchyStressNXX, 
const ScalarT* unrotatedCauchyStressNXY, 
const ScalarT* unrotatedCauchyStressNXZ, 
const ScalarT* unrotatedCauchyStressNYX, 
const ScalarT* unrotatedCauchyStressNYY, 
const ScalarT* unrotatedCauchyStressNYZ, 
const ScalarT* unrotatedCauchyStressNZX, 
const ScalarT* unrotatedCauchyStressNZY, 
const ScalarT* unrotatedCauchyStressNZZ, 
ScalarT* unrotatedCauchyStressNP1XX, 
ScalarT* unrotatedCauchyStressNP1XY, 
ScalarT* unrotatedCauchyStressNP1XZ, 
ScalarT* unrotatedCauchyStressNP1YX, 
ScalarT* unrotatedCauchyStressNP1YY, 
ScalarT* unrotatedCauchyStressNP1YZ, 
ScalarT* unrotatedCauchyStressNP1ZX, 
ScalarT* unrotatedCauchyStressNP1ZY, 
ScalarT* unrotatedCauchyStressNP1ZZ, 
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
  const ScalarT* sigmaNXX = unrotatedCauchyStressNXX;
  const ScalarT* sigmaNXY = unrotatedCauchyStressNXY;
  const ScalarT* sigmaNXZ = unrotatedCauchyStressNXZ;
  const ScalarT* sigmaNYX = unrotatedCauchyStressNYX;
  const ScalarT* sigmaNYY = unrotatedCauchyStressNYY;
  const ScalarT* sigmaNYZ = unrotatedCauchyStressNYZ;
  const ScalarT* sigmaNZX = unrotatedCauchyStressNZX;
  const ScalarT* sigmaNZY = unrotatedCauchyStressNZY;
  const ScalarT* sigmaNZZ = unrotatedCauchyStressNZZ;
  ScalarT* sigmaNP1XX = unrotatedCauchyStressNP1XX;
  ScalarT* sigmaNP1XY = unrotatedCauchyStressNP1XY;
  ScalarT* sigmaNP1XZ = unrotatedCauchyStressNP1XZ;
  ScalarT* sigmaNP1YX = unrotatedCauchyStressNP1YX;
  ScalarT* sigmaNP1YY = unrotatedCauchyStressNP1YY;
  ScalarT* sigmaNP1YZ = unrotatedCauchyStressNP1YZ;
  ScalarT* sigmaNP1ZX = unrotatedCauchyStressNP1ZX;
  ScalarT* sigmaNP1ZY = unrotatedCauchyStressNP1ZY;
  ScalarT* sigmaNP1ZZ = unrotatedCauchyStressNP1ZZ;

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
const double* unrotatedCauchyStressNXX, 
const double* unrotatedCauchyStressNXY, 
const double* unrotatedCauchyStressNXZ, 
const double* unrotatedCauchyStressNYX, 
const double* unrotatedCauchyStressNYY, 
const double* unrotatedCauchyStressNYZ, 
const double* unrotatedCauchyStressNZX, 
const double* unrotatedCauchyStressNZY, 
const double* unrotatedCauchyStressNZZ, 
double* unrotatedCauchyStressNP1XX, 
double* unrotatedCauchyStressNP1XY, 
double* unrotatedCauchyStressNP1XZ, 
double* unrotatedCauchyStressNP1YX, 
double* unrotatedCauchyStressNP1YY, 
double* unrotatedCauchyStressNP1YZ, 
double* unrotatedCauchyStressNP1ZX, 
double* unrotatedCauchyStressNP1ZY, 
double* unrotatedCauchyStressNP1ZZ, 
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
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNXX, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNXY, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNXZ, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNYX, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNYY, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNYZ, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNZX, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNZY, 
const Sacado::Fad::DFad<double>* unrotatedCauchyStressNZZ, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1XX, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1XY, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1XZ, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1YX, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1YY, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1YZ, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1ZX, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1ZY, 
Sacado::Fad::DFad<double>* unrotatedCauchyStressNP1ZZ, 
int numPoints, 
double bulkMod,
double shearMod,
double dt
);

}
