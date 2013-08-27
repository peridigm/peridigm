//! \file correspondence.cxx

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
#include "correspondence.h"
#include "material_utilities.h"

namespace CORRESPONDENCE {

template<typename ScalarT>
void computeGreenLagrangeStrain
(
  const ScalarT* deformationGradientXX,
  const ScalarT* deformationGradientXY,
  const ScalarT* deformationGradientXZ,
  const ScalarT* deformationGradientYX,
  const ScalarT* deformationGradientYY,
  const ScalarT* deformationGradientYZ,
  const ScalarT* deformationGradientZX,
  const ScalarT* deformationGradientZY,
  const ScalarT* deformationGradientZZ,
  ScalarT* greenLagrangeStrainXX,
  ScalarT* greenLagrangeStrainXY,
  ScalarT* greenLagrangeStrainXZ,
  ScalarT* greenLagrangeStrainYX,
  ScalarT* greenLagrangeStrainYY,
  ScalarT* greenLagrangeStrainYZ,
  ScalarT* greenLagrangeStrainZX,
  ScalarT* greenLagrangeStrainZY,
  ScalarT* greenLagrangeStrainZZ,
  int numOwnedPoints
)
{
  // Green-Lagrange Strain E = 0.5*(F^T F - I)

  const ScalarT* defGradXX = deformationGradientXX;
  const ScalarT* defGradXY = deformationGradientXY;
  const ScalarT* defGradXZ = deformationGradientXZ;
  const ScalarT* defGradYX = deformationGradientYX;
  const ScalarT* defGradYY = deformationGradientYY;
  const ScalarT* defGradYZ = deformationGradientYZ;
  const ScalarT* defGradZX = deformationGradientZX;
  const ScalarT* defGradZY = deformationGradientZY;
  const ScalarT* defGradZZ = deformationGradientZZ;
  ScalarT* strainXX = greenLagrangeStrainXX;
  ScalarT* strainXY = greenLagrangeStrainXY;
  ScalarT* strainXZ = greenLagrangeStrainXZ;
  ScalarT* strainYX = greenLagrangeStrainYX;
  ScalarT* strainYY = greenLagrangeStrainYY;
  ScalarT* strainYZ = greenLagrangeStrainYZ;
  ScalarT* strainZX = greenLagrangeStrainZX;
  ScalarT* strainZY = greenLagrangeStrainZY;
  ScalarT* strainZZ = greenLagrangeStrainZZ;

  for(int iID=0 ; iID<numOwnedPoints ; ++iID, 
        ++defGradXX, ++defGradXY, ++defGradXZ,
        ++defGradYX, ++defGradYY, ++defGradYZ,
        ++defGradZX, ++defGradZY, ++defGradZZ,
        ++strainXX, ++strainXY, ++strainXZ,
        ++strainYX, ++strainYY, ++strainYZ,
        ++strainZX, ++strainZY, ++strainZZ){

    *strainXX = 0.5 * ( (*defGradXX) * (*defGradXX) + (*defGradYX) * (*defGradYX) + (*defGradZX) * (*defGradZX) - 1.0 );
    *strainXY = 0.5 * ( (*defGradXX) * (*defGradXY) + (*defGradYX) * (*defGradYY) + (*defGradZX) * (*defGradZY) );
    *strainXZ = 0.5 * ( (*defGradXX) * (*defGradXZ) + (*defGradYX) * (*defGradYZ) + (*defGradZX) * (*defGradZZ) );
    *strainYX = 0.5 * ( (*defGradXY) * (*defGradXX) + (*defGradYY) * (*defGradYX) + (*defGradZY) * (*defGradZX) );
    *strainYY = 0.5 * ( (*defGradXY) * (*defGradXY) + (*defGradYY) * (*defGradYY) + (*defGradZY) * (*defGradZY) - 1.0 );
    *strainYZ = 0.5 * ( (*defGradXY) * (*defGradXZ) + (*defGradYY) * (*defGradYZ) + (*defGradZY) * (*defGradZZ) );
    *strainZX = 0.5 * ( (*defGradXZ) * (*defGradXX) + (*defGradYZ) * (*defGradYX) + (*defGradZZ) * (*defGradZX) );
    *strainZY = 0.5 * ( (*defGradXZ) * (*defGradXY) + (*defGradYZ) * (*defGradYY) + (*defGradZZ) * (*defGradZY) );
    *strainZZ = 0.5 * ( (*defGradXZ) * (*defGradXZ) + (*defGradYZ) * (*defGradYZ) + (*defGradZZ) * (*defGradZZ) - 1.0 );
  }
}

/** Explicit template instantiation for double. */
template void computeGreenLagrangeStrain<double>
(
  const double* deformationGradientXX,
  const double* deformationGradientXY,
  const double* deformationGradientXZ,
  const double* deformationGradientYX,
  const double* deformationGradientYY,
  const double* deformationGradientYZ,
  const double* deformationGradientZX,
  const double* deformationGradientZY,
  const double* deformationGradientZZ,
  double* greenLagrangeStrainXX,
  double* greenLagrangeStrainXY,
  double* greenLagrangeStrainXZ,
  double* greenLagrangeStrainYX,
  double* greenLagrangeStrainYY,
  double* greenLagrangeStrainYZ,
  double* greenLagrangeStrainZX,
  double* greenLagrangeStrainZY,
  double* greenLagrangeStrainZZ,
  int numOwnedPoints
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeGreenLagrangeStrain<Sacado::Fad::DFad<double>>
(
  const Sacado::Fad::DFad<double>* deformationGradientXX,
  const Sacado::Fad::DFad<double>* deformationGradientXY,
  const Sacado::Fad::DFad<double>* deformationGradientXZ,
  const Sacado::Fad::DFad<double>* deformationGradientYX,
  const Sacado::Fad::DFad<double>* deformationGradientYY,
  const Sacado::Fad::DFad<double>* deformationGradientYZ,
  const Sacado::Fad::DFad<double>* deformationGradientZX,
  const Sacado::Fad::DFad<double>* deformationGradientZY,
  const Sacado::Fad::DFad<double>* deformationGradientZZ,
  Sacado::Fad::DFad<double>* greenLagrangeStrainXX,
  Sacado::Fad::DFad<double>* greenLagrangeStrainXY,
  Sacado::Fad::DFad<double>* greenLagrangeStrainXZ,
  Sacado::Fad::DFad<double>* greenLagrangeStrainYX,
  Sacado::Fad::DFad<double>* greenLagrangeStrainYY,
  Sacado::Fad::DFad<double>* greenLagrangeStrainYZ,
  Sacado::Fad::DFad<double>* greenLagrangeStrainZX,
  Sacado::Fad::DFad<double>* greenLagrangeStrainZY,
  Sacado::Fad::DFad<double>* greenLagrangeStrainZZ,
  int numOwnedPoints
);

}
