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

//! Invert a single 3-by-3 matrix; returns zero of successful, one if not successful (e.g., singular matrix).
template<typename ScalarT>
int invert3by3Matrix
(
 const ScalarT& matrixXX,
 const ScalarT& matrixXY,
 const ScalarT& matrixXZ,
 const ScalarT& matrixYX,
 const ScalarT& matrixYY,
 const ScalarT& matrixYZ,
 const ScalarT& matrixZX,
 const ScalarT& matrixZY,
 const ScalarT& matrixZZ,
 ScalarT& inverseXX,
 ScalarT& inverseXY,
 ScalarT& inverseXZ,
 ScalarT& inverseYX,
 ScalarT& inverseYY,
 ScalarT& inverseYZ,
 ScalarT& inverseZX,
 ScalarT& inverseZY,
 ScalarT& inverseZZ
);

//! Inner product of two 3-by-3 matrices.
template<typename ScalarT>
void MatrixMultiply
(
 const ScalarT& aXX,
 const ScalarT& aXY,
 const ScalarT& aXZ,
 const ScalarT& aYX,
 const ScalarT& aYY,
 const ScalarT& aYZ,
 const ScalarT& aZX,
 const ScalarT& aZY,
 const ScalarT& aZZ,
 const ScalarT& bXX,
 const ScalarT& bXY,
 const ScalarT& bXZ,
 const ScalarT& bYX,
 const ScalarT& bYY,
 const ScalarT& bYZ,
 const ScalarT& bZX,
 const ScalarT& bZY,
 const ScalarT& bZZ,
 ScalarT& resultXX,
 ScalarT& resultXY,
 ScalarT& resultXZ,
 ScalarT& resultYX,
 ScalarT& resultYY,
 ScalarT& resultYZ,
 ScalarT& resultZX,
 ScalarT& resultZY,
 ScalarT& resultZZ
);

template<typename ScalarT>
void MatrixUpdate
(
 const ScalarT& alpha,
 const ScalarT& beta,
 const ScalarT& aXX,
 const ScalarT& aXY,
 const ScalarT& aXZ,
 const ScalarT& aYX,
 const ScalarT& aYY,
 const ScalarT& aYZ,
 const ScalarT& aZX,
 const ScalarT& aZY,
 const ScalarT& aZZ,
 const ScalarT& bXX,
 const ScalarT& bXY,
 const ScalarT& bXZ,
 const ScalarT& bYX,
 const ScalarT& bYY,
 const ScalarT& bYZ,
 const ScalarT& bZX,
 const ScalarT& bZY,
 const ScalarT& bZZ,
 ScalarT& resultXX,
 ScalarT& resultXY,
 ScalarT& resultXZ,
 ScalarT& resultYX,
 ScalarT& resultYY,
 ScalarT& resultYZ,
 ScalarT& resultZX,
 ScalarT& resultZY,
 ScalarT& resultZZ
);

template<typename ScalarT>
void UnrotateTensor
(
 const ScalarT& aXX,
 const ScalarT& aXY,
 const ScalarT& aXZ,
 const ScalarT& aYX,
 const ScalarT& aYY,
 const ScalarT& aYZ,
 const ScalarT& aZX,
 const ScalarT& aZY,
 const ScalarT& aZZ,
 const ScalarT& rXX,
 const ScalarT& rXY,
 const ScalarT& rXZ,
 const ScalarT& rYX,
 const ScalarT& rYY,
 const ScalarT& rYZ,
 const ScalarT& rZX,
 const ScalarT& rZY,
 const ScalarT& rZZ,
 ScalarT& resultXX,
 ScalarT& resultXY,
 ScalarT& resultXZ,
 ScalarT& resultYX,
 ScalarT& resultYY,
 ScalarT& resultYZ,
 ScalarT& resultZX,
 ScalarT& resultZY,
 ScalarT& resultZZ
 );

template<typename ScalarT>
void RotateTensor
(
 const ScalarT& aXX,
 const ScalarT& aXY,
 const ScalarT& aXZ,
 const ScalarT& aYX,
 const ScalarT& aYY,
 const ScalarT& aYZ,
 const ScalarT& aZX,
 const ScalarT& aZY,
 const ScalarT& aZZ,
 const ScalarT& rXX,
 const ScalarT& rXY,
 const ScalarT& rXZ,
 const ScalarT& rYX,
 const ScalarT& rYY,
 const ScalarT& rYZ,
 const ScalarT& rZX,
 const ScalarT& rZY,
 const ScalarT& rZZ,
 ScalarT& resultXX,
 ScalarT& resultXY,
 ScalarT& resultXZ,
 ScalarT& resultYX,
 ScalarT& resultYY,
 ScalarT& resultYZ,
 ScalarT& resultZX,
 ScalarT& resultZY,
 ScalarT& resultZZ
);

/* template<typename ScalarT> */
/* void unrotateCauchyStress */
/* ( */
/*  const ScalarT* rotationTensorXX, */
/*  const ScalarT* rotationTensorXY, */
/*  const ScalarT* rotationTensorXZ, */
/*  const ScalarT* rotationTensorYX, */
/*  const ScalarT* rotationTensorYY, */
/*  const ScalarT* rotationTensorYZ, */
/*  const ScalarT* rotationTensorZX, */
/*  const ScalarT* rotationTensorZY, */
/*  const ScalarT* rotationTensorZZ, */
/*  const ScalarT* rotatedCauchyStressXX, */
/*  const ScalarT* rotatedCauchyStressXY, */
/*  const ScalarT* rotatedCauchyStressXZ, */
/*  const ScalarT* rotatedCauchyStressYX, */
/*  const ScalarT* rotatedCauchyStressYY, */
/*  const ScalarT* rotatedCauchyStressYZ, */
/*  const ScalarT* rotatedCauchyStressZX, */
/*  const ScalarT* rotatedCauchyStressZY, */
/*  const ScalarT* rotatedCauchyStressZZ, */
/*  ScalarT* unrotatedCauchyStressXX, */
/*  ScalarT* unrotatedCauchyStressXY, */
/*  ScalarT* unrotatedCauchyStressXZ, */
/*  ScalarT* unrotatedCauchyStressYX, */
/*  ScalarT* unrotatedCauchyStressYY, */
/*  ScalarT* unrotatedCauchyStressYZ, */
/*  ScalarT* unrotatedCauchyStressZX, */
/*  ScalarT* unrotatedCauchyStressZY, */
/*  ScalarT* unrotatedCauchyStressZZ, */
/*  int numPoints */
/* ); */

template<typename ScalarT>
void rotateCauchyStress
(
 const ScalarT* rotationTensorXX,
 const ScalarT* rotationTensorXY,
 const ScalarT* rotationTensorXZ,
 const ScalarT* rotationTensorYX,
 const ScalarT* rotationTensorYY,
 const ScalarT* rotationTensorYZ,
 const ScalarT* rotationTensorZX,
 const ScalarT* rotationTensorZY,
 const ScalarT* rotationTensorZZ,
 const ScalarT* unrotatedCauchyStressXX,
 const ScalarT* unrotatedCauchyStressXY,
 const ScalarT* unrotatedCauchyStressXZ,
 const ScalarT* unrotatedCauchyStressYX,
 const ScalarT* unrotatedCauchyStressYY,
 const ScalarT* unrotatedCauchyStressYZ,
 const ScalarT* unrotatedCauchyStressZX,
 const ScalarT* unrotatedCauchyStressZY,
 const ScalarT* unrotatedCauchyStressZZ,
 ScalarT* rotatedCauchyStressXX,
 ScalarT* rotatedCauchyStressXY,
 ScalarT* rotatedCauchyStressXZ,
 ScalarT* rotatedCauchyStressYX,
 ScalarT* rotatedCauchyStressYY,
 ScalarT* rotatedCauchyStressYZ,
 ScalarT* rotatedCauchyStressZX,
 ScalarT* rotatedCauchyStressZY,
 ScalarT* rotatedCauchyStressZZ,
 int numPoints
);

template<typename ScalarT>
int computeShapeTensorInverseAndApproximateDeformationGradient
(
const double* volume,
const double* modelCoordinates,
const ScalarT* coordinates,
ScalarT* shapeTensorInverseXX,
ScalarT* shapeTensorInverseXY,
ScalarT* shapeTensorInverseXZ,
ScalarT* shapeTensorInverseYX,
ScalarT* shapeTensorInverseYY,
ScalarT* shapeTensorInverseYZ,
ScalarT* shapeTensorInverseZX,
ScalarT* shapeTensorInverseZY,
ScalarT* shapeTensorInverseZZ,
ScalarT* deformationGradientXX,
ScalarT* deformationGradientXY,
ScalarT* deformationGradientXZ,
ScalarT* deformationGradientYX,
ScalarT* deformationGradientYY,
ScalarT* deformationGradientYZ,
ScalarT* deformationGradientZX,
ScalarT* deformationGradientZY,
ScalarT* deformationGradientZZ,
const int* neighborhoodList,
int numPoints,
double horizon
);

// Calculation of stretch rates following Flanagan & Taylor
template<typename ScalarT>
int computeUnrotatedRateOfDeformationAndRotationTensor(
const double* volume,
const double* modelCoordinates,
const ScalarT* velocities,
const ScalarT* deformationGradientXX,
const ScalarT* deformationGradientXY,
const ScalarT* deformationGradientXZ,
const ScalarT* deformationGradientYX,
const ScalarT* deformationGradientYY,
const ScalarT* deformationGradientYZ,
const ScalarT* deformationGradientZX,
const ScalarT* deformationGradientZY,
const ScalarT* deformationGradientZZ,
const ScalarT* shapeTensorInverseXX,
const ScalarT* shapeTensorInverseXY,
const ScalarT* shapeTensorInverseXZ,
const ScalarT* shapeTensorInverseYX,
const ScalarT* shapeTensorInverseYY,
const ScalarT* shapeTensorInverseYZ,
const ScalarT* shapeTensorInverseZX,
const ScalarT* shapeTensorInverseZY,
const ScalarT* shapeTensorInverseZZ,
const ScalarT* leftStretchTensorNXX,
const ScalarT* leftStretchTensorNXY,
const ScalarT* leftStretchTensorNXZ,
const ScalarT* leftStretchTensorNYX,
const ScalarT* leftStretchTensorNYY,
const ScalarT* leftStretchTensorNYZ,
const ScalarT* leftStretchTensorNZX,
const ScalarT* leftStretchTensorNZY,
const ScalarT* leftStretchTensorNZZ,
const ScalarT* rotationTensorNXX,
const ScalarT* rotationTensorNXY,
const ScalarT* rotationTensorNXZ,
const ScalarT* rotationTensorNYX,
const ScalarT* rotationTensorNYY,
const ScalarT* rotationTensorNYZ,
const ScalarT* rotationTensorNZX,
const ScalarT* rotationTensorNZY,
const ScalarT* rotationTensorNZZ,
ScalarT* leftStretchTensorNP1XX,
ScalarT* leftStretchTensorNP1XY,
ScalarT* leftStretchTensorNP1XZ,
ScalarT* leftStretchTensorNP1YX,
ScalarT* leftStretchTensorNP1YY,
ScalarT* leftStretchTensorNP1YZ,
ScalarT* leftStretchTensorNP1ZX,
ScalarT* leftStretchTensorNP1ZY,
ScalarT* leftStretchTensorNP1ZZ,
ScalarT* rotationTensorNP1XX,
ScalarT* rotationTensorNP1XY,
ScalarT* rotationTensorNP1XZ,
ScalarT* rotationTensorNP1YX,
ScalarT* rotationTensorNP1YY,
ScalarT* rotationTensorNP1YZ,
ScalarT* rotationTensorNP1ZX,
ScalarT* rotationTensorNP1ZY,
ScalarT* rotationTensorNP1ZZ,
ScalarT* unrotatedRateOfDeformationXX,
ScalarT* unrotatedRateOfDeformationXY,
ScalarT* unrotatedRateOfDeformationXZ,
ScalarT* unrotatedRateOfDeformationYX,
ScalarT* unrotatedRateOfDeformationYY,
ScalarT* unrotatedRateOfDeformationYZ,
ScalarT* unrotatedRateOfDeformationZX,
ScalarT* unrotatedRateOfDeformationZY,
ScalarT* unrotatedRateOfDeformationZZ,
const int* neighborhoodList,
int numPoints,
double horizon,
double dt
);

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

template<typename ScalarT>
void computeHourglassForce
(
const double* volume,
const double* modelCoordinates,
const ScalarT* coordinates,
const ScalarT* deformationGradientXX,
const ScalarT* deformationGradientXY,
const ScalarT* deformationGradientXZ,
const ScalarT* deformationGradientYX,
const ScalarT* deformationGradientYY,
const ScalarT* deformationGradientYZ,
const ScalarT* deformationGradientZX,
const ScalarT* deformationGradientZY,
const ScalarT* deformationGradientZZ,
ScalarT* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double horizon,
double bulkModulus,
double hourglassCoefficient
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
