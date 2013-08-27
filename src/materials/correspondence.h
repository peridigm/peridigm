//! \file correspondence.h

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
#ifndef CORRESPONDENCE_H
#define CORRESPONDENCE_H

namespace CORRESPONDENCE {
  
//! Green-Lagrange Strain E = 0.5*(F^T F - I).
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
  int numPoints
);

//! Classical elastic material model stress calculation (Hooke's law).
template<typename ScalarT>
void computeClassicalElasticStress
(
 const ScalarT* strainXX,
 const ScalarT* strainXY,
 const ScalarT* strainXZ,
 const ScalarT* strainYX,
 const ScalarT* strainYY,
 const ScalarT* strainYZ,
 const ScalarT* strainZX,
 const ScalarT* strainZY,
 const ScalarT* strainZZ,
 ScalarT* cauchyStressXX,
 ScalarT* cauchyStressXY,
 ScalarT* cauchyStressXZ,
 ScalarT* cauchyStressYX,
 ScalarT* cauchyStressYY,
 ScalarT* cauchyStressYZ,
 ScalarT* cauchyStressZX,
 ScalarT* cauchyStressZY,
 ScalarT* cauchyStressZZ,
 int numPoints,
 double youngsModulus,
 double poissonsRatio
);

}

#endif // CORRESPONDENCE_H
