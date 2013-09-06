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

#include "correspondence.h"
#include <Sacado.hpp>

namespace CORRESPONDENCE {

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
)
{
  int returnCode(0);

  ScalarT minor0 =  matrixYY * matrixZZ - matrixYZ * matrixZY;
  ScalarT minor1 =  matrixYX * matrixZZ - matrixYZ * matrixZX;
  ScalarT minor2 =  matrixYX * matrixZY - matrixYY * matrixZX;
  ScalarT minor3 =  matrixXY * matrixZZ - matrixXZ * matrixZY;
  ScalarT minor4 =  matrixXX * matrixZZ - matrixZX * matrixXZ;
  ScalarT minor5 =  matrixXX * matrixZY - matrixXY * matrixZX;
  ScalarT minor6 =  matrixXY * matrixYZ - matrixXZ * matrixYY;
  ScalarT minor7 =  matrixXX * matrixYZ - matrixXZ * matrixYX;
  ScalarT minor8 =  matrixXX * matrixYY - matrixXY * matrixYX;
  ScalarT det = matrixXX * minor0 - matrixXY * minor1 + matrixXZ * minor2;

  if(det == ScalarT(0.0)){
    returnCode = 1;
    inverseXX = 0.0;
    inverseXY = 0.0;
    inverseXZ = 0.0;
    inverseYX = 0.0;
    inverseYY = 0.0;
    inverseYZ = 0.0;
    inverseZX = 0.0;
    inverseZY = 0.0;
    inverseZZ = 0.0;
  }
  else{
    inverseXX = minor0/det;
    inverseXY = -1.0*minor3/det;
    inverseXZ = minor6/det;
    inverseYX = -1.0*minor1/det;
    inverseYY = minor4/det;
    inverseYZ = -1.0*minor7/det;
    inverseZX = minor2/det;
    inverseZY = -1.0*minor5/det;
    inverseZZ = minor8/det;
  }

  return returnCode;
}

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
)
{
  resultXX = aXX * bXX + aXY * bYX + aXZ * bZX;
  resultXY = aXX * bXY + aXY * bYY + aXZ * bZY;
  resultXZ = aXX * bXZ + aXY * bYZ + aXZ * bZZ;
  resultYX = aYX * bXX + aYY * bYX + aYZ * bZX;
  resultYY = aYX * bXY + aYY * bYY + aYZ * bZY;
  resultYZ = aYX * bXZ + aYY * bYZ + aYZ * bZZ;
  resultZX = aZX * bXX + aZY * bYX + aZZ * bZX;
  resultZY = aZX * bXY + aZY * bYY + aZZ * bZY;
  resultZZ = aZX * bXZ + aZY * bYZ + aZZ * bZZ;
}

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

  for(int iID=0 ; iID<numPoints ; ++iID, 
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
)
{
  double undeformedBond[3], undeformedBondLength;
  ScalarT deformedBond[3], deformedBondLength;
  ScalarT expectedNeighborLocation[3], hourglassVector[3], bondDamage, dot, magnitude;
  double vol, neighborVol;
  int neighborIndex;

  double constant = 18.0*hourglassCoefficient*bulkModulus/(3.1415926536*horizon*horizon*horizon*horizon);

  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

    int numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){
      neighborIndex = *neighborListPtr;

      undeformedBond[0] = modelCoordinates[3*neighborIndex]   - modelCoordinates[3*iID];
      undeformedBond[1] = modelCoordinates[3*neighborIndex+1] - modelCoordinates[3*iID+1];
      undeformedBond[2] = modelCoordinates[3*neighborIndex+2] - modelCoordinates[3*iID+2];
      undeformedBondLength = sqrt(undeformedBond[0]*undeformedBond[0] +
                                  undeformedBond[1]*undeformedBond[1] +
                                  undeformedBond[2]*undeformedBond[2]);

      deformedBond[0] = coordinates[3*neighborIndex]   - coordinates[3*iID];
      deformedBond[1] = coordinates[3*neighborIndex+1] - coordinates[3*iID+1];
      deformedBond[2] = coordinates[3*neighborIndex+2] - coordinates[3*iID+2];
      deformedBondLength = sqrt(deformedBond[0]*deformedBond[0] +
				deformedBond[1]*deformedBond[1] +
				deformedBond[2]*deformedBond[2]);

      expectedNeighborLocation[0] = coordinates[3*iID] +
	deformationGradientXX[iID]*undeformedBond[0] +
	deformationGradientXY[iID]*undeformedBond[1] +
	deformationGradientXZ[iID]*undeformedBond[2];
      expectedNeighborLocation[1] = coordinates[3*iID+1] +
	deformationGradientYX[iID]*undeformedBond[0] +
	deformationGradientYY[iID]*undeformedBond[1] +
	deformationGradientYZ[iID]*undeformedBond[2];
      expectedNeighborLocation[2] = coordinates[3*iID+2] +
	deformationGradientZX[iID]*undeformedBond[0] +
	deformationGradientZY[iID]*undeformedBond[1] +
	deformationGradientZZ[iID]*undeformedBond[2];

      // \todo Include bond damage in hourglass force calculation
      bondDamage = 0.0;

      hourglassVector[0] = expectedNeighborLocation[0] - coordinates[3*neighborIndex];
      hourglassVector[1] = expectedNeighborLocation[1] - coordinates[3*neighborIndex+1];
      hourglassVector[2] = expectedNeighborLocation[2] - coordinates[3*neighborIndex+2];

      dot = -1.0 * (hourglassVector[0]*deformedBond[0] + hourglassVector[1]*deformedBond[1] + hourglassVector[2]*deformedBond[2]);

      magnitude = (1.0-bondDamage) * constant * (dot/undeformedBondLength) * (1.0/deformedBondLength);

      vol = volume[iID];
      neighborVol = volume[neighborIndex];

      hourglassForceDensity[3*iID]   += magnitude * deformedBond[0] * neighborVol;
      hourglassForceDensity[3*iID+1] += magnitude * deformedBond[1] * neighborVol;
      hourglassForceDensity[3*iID+2] += magnitude * deformedBond[2] * neighborVol;
      hourglassForceDensity[3*neighborIndex]   -= magnitude * deformedBond[0] * vol;
      hourglassForceDensity[3*neighborIndex+1] -= magnitude * deformedBond[1] * vol;
      hourglassForceDensity[3*neighborIndex+2] -= magnitude * deformedBond[2] * vol;

    }
  }
}

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
)
{
  // Hooke's law

  const ScalarT* epsilonXX = strainXX;
  const ScalarT* epsilonXY = strainXY;
  const ScalarT* epsilonXZ = strainXZ;
  const ScalarT* epsilonYX = strainYX;
  const ScalarT* epsilonYY = strainYY;
  const ScalarT* epsilonYZ = strainYZ;
  const ScalarT* epsilonZX = strainZX;
  const ScalarT* epsilonZY = strainZY;
  const ScalarT* epsilonZZ = strainZZ;
  ScalarT* sigmaXX = cauchyStressXX;
  ScalarT* sigmaXY = cauchyStressXY;
  ScalarT* sigmaXZ = cauchyStressXZ;
  ScalarT* sigmaYX = cauchyStressYX;
  ScalarT* sigmaYY = cauchyStressYY;
  ScalarT* sigmaYZ = cauchyStressYZ;
  ScalarT* sigmaZX = cauchyStressZX;
  ScalarT* sigmaZY = cauchyStressZY;
  ScalarT* sigmaZZ = cauchyStressZZ;

  double constant = youngsModulus/((1.0 + poissonsRatio)*(1.0 - 2.0*poissonsRatio));

  for(int iID=0 ; iID<numPoints ; ++iID, 
        ++epsilonXX, ++epsilonXY, ++epsilonXZ,
        ++epsilonYX, ++epsilonYY, ++epsilonYZ,
        ++epsilonZX, ++epsilonZY, ++epsilonZZ,
        ++sigmaXX, ++sigmaXY, ++sigmaXZ,
        ++sigmaYX, ++sigmaYY, ++sigmaYZ,
        ++sigmaZX, ++sigmaZY, ++sigmaZZ){

    *sigmaXX = constant * ( (1.0 - poissonsRatio) * (*epsilonXX)  +        poissonsRatio  * (*epsilonYY) +        poissonsRatio  * (*epsilonZZ) );
    *sigmaYY = constant * (        poissonsRatio  * (*epsilonXX)  + (1.0 - poissonsRatio) * (*epsilonYY) +        poissonsRatio  * (*epsilonZZ) );
    *sigmaZZ = constant * (        poissonsRatio  * (*epsilonXX)  +        poissonsRatio  * (*epsilonYY) + (1.0 - poissonsRatio) * (*epsilonZZ) );
    *sigmaXY = constant * (1.0 - 2.0*poissonsRatio) * (*epsilonXY);
    *sigmaYZ = constant * (1.0 - 2.0*poissonsRatio) * (*epsilonYZ);
    *sigmaZX = constant * (1.0 - 2.0*poissonsRatio) * (*epsilonZX);
    *sigmaYX = *sigmaXY;
    *sigmaZY = *sigmaYZ;
    *sigmaXZ = *sigmaZX;
  }
}

/** Explicit template instantiation for double. */
template void MatrixMultiply<double>
(
 const double& aXX,
 const double& aXY,
 const double& aXZ,
 const double& aYX,
 const double& aYY,
 const double& aYZ,
 const double& aZX,
 const double& aZY,
 const double& aZZ,
 const double& bXX,
 const double& bXY,
 const double& bXZ,
 const double& bYX,
 const double& bYY,
 const double& bYZ,
 const double& bZX,
 const double& bZY,
 const double& bZZ,
 double& resultXX,
 double& resultXY,
 double& resultXZ,
 double& resultYX,
 double& resultYY,
 double& resultYZ,
 double& resultZX,
 double& resultZY,
 double& resultZZ
);

template int invert3by3Matrix<double>
(
 const double& matrixXX,
 const double& matrixXY,
 const double& matrixXZ,
 const double& matrixYX,
 const double& matrixYY,
 const double& matrixYZ,
 const double& matrixZX,
 const double& matrixZY,
 const double& matrixZZ,
 double& inverseXX,
 double& inverseXY,
 double& inverseXZ,
 double& inverseYX,
 double& inverseYY,
 double& inverseYZ,
 double& inverseZX,
 double& inverseZY,
 double& inverseZZ
);

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
  int numPoints
);

template void computeHourglassForce<double>
(
const double* volume,
const double* modelCoordinates,
const double* coordinates,
const double* deformationGradientXX,
const double* deformationGradientXY,
const double* deformationGradientXZ,
const double* deformationGradientYX,
const double* deformationGradientYY,
const double* deformationGradientYZ,
const double* deformationGradientZX,
const double* deformationGradientZY,
const double* deformationGradientZZ,
double* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double horizon,
double bulkModulus,
double hourglassCoefficient
);

template void computeClassicalElasticStress<double>
(
 const double* strainXX,
 const double* strainXY,
 const double* strainXZ,
 const double* strainYX,
 const double* strainYY,
 const double* strainYZ,
 const double* strainZX,
 const double* strainZY,
 const double* strainZZ,
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
 double youngsModulus,
 double poissonsRatio
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void MatrixMultiply<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>& aXX,
 const Sacado::Fad::DFad<double>& aXY,
 const Sacado::Fad::DFad<double>& aXZ,
 const Sacado::Fad::DFad<double>& aYX,
 const Sacado::Fad::DFad<double>& aYY,
 const Sacado::Fad::DFad<double>& aYZ,
 const Sacado::Fad::DFad<double>& aZX,
 const Sacado::Fad::DFad<double>& aZY,
 const Sacado::Fad::DFad<double>& aZZ,
 const Sacado::Fad::DFad<double>& bXX,
 const Sacado::Fad::DFad<double>& bXY,
 const Sacado::Fad::DFad<double>& bXZ,
 const Sacado::Fad::DFad<double>& bYX,
 const Sacado::Fad::DFad<double>& bYY,
 const Sacado::Fad::DFad<double>& bYZ,
 const Sacado::Fad::DFad<double>& bZX,
 const Sacado::Fad::DFad<double>& bZY,
 const Sacado::Fad::DFad<double>& bZZ,
 Sacado::Fad::DFad<double>& resultXX,
 Sacado::Fad::DFad<double>& resultXY,
 Sacado::Fad::DFad<double>& resultXZ,
 Sacado::Fad::DFad<double>& resultYX,
 Sacado::Fad::DFad<double>& resultYY,
 Sacado::Fad::DFad<double>& resultYZ,
 Sacado::Fad::DFad<double>& resultZX,
 Sacado::Fad::DFad<double>& resultZY,
 Sacado::Fad::DFad<double>& resultZZ
);

template int invert3by3Matrix<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>& matrixXX,
 const Sacado::Fad::DFad<double>& matrixXY,
 const Sacado::Fad::DFad<double>& matrixXZ,
 const Sacado::Fad::DFad<double>& matrixYX,
 const Sacado::Fad::DFad<double>& matrixYY,
 const Sacado::Fad::DFad<double>& matrixYZ,
 const Sacado::Fad::DFad<double>& matrixZX,
 const Sacado::Fad::DFad<double>& matrixZY,
 const Sacado::Fad::DFad<double>& matrixZZ,
 Sacado::Fad::DFad<double>& inverseXX,
 Sacado::Fad::DFad<double>& inverseXY,
 Sacado::Fad::DFad<double>& inverseXZ,
 Sacado::Fad::DFad<double>& inverseYX,
 Sacado::Fad::DFad<double>& inverseYY,
 Sacado::Fad::DFad<double>& inverseYZ,
 Sacado::Fad::DFad<double>& inverseZX,
 Sacado::Fad::DFad<double>& inverseZY,
 Sacado::Fad::DFad<double>& inverseZZ
);

template void computeGreenLagrangeStrain<Sacado::Fad::DFad<double> >
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
  int numPoints
);

template void computeHourglassForce<Sacado::Fad::DFad<double> >
(
const double* volume,
const double* modelCoordinates,
const Sacado::Fad::DFad<double>* coordinates,
const Sacado::Fad::DFad<double>* deformationGradientXX,
const Sacado::Fad::DFad<double>* deformationGradientXY,
const Sacado::Fad::DFad<double>* deformationGradientXZ,
const Sacado::Fad::DFad<double>* deformationGradientYX,
const Sacado::Fad::DFad<double>* deformationGradientYY,
const Sacado::Fad::DFad<double>* deformationGradientYZ,
const Sacado::Fad::DFad<double>* deformationGradientZX,
const Sacado::Fad::DFad<double>* deformationGradientZY,
const Sacado::Fad::DFad<double>* deformationGradientZZ,
Sacado::Fad::DFad<double>* hourglassForceDensity,
const int* neighborhoodList,
int numPoints,
double horizon,
double bulkModulus,
double hourglassCoefficient
);

template void computeClassicalElasticStress<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* strainXX,
 const Sacado::Fad::DFad<double>* strainXY,
 const Sacado::Fad::DFad<double>* strainXZ,
 const Sacado::Fad::DFad<double>* strainYX,
 const Sacado::Fad::DFad<double>* strainYY,
 const Sacado::Fad::DFad<double>* strainYZ,
 const Sacado::Fad::DFad<double>* strainZX,
 const Sacado::Fad::DFad<double>* strainZY,
 const Sacado::Fad::DFad<double>* strainZZ,
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
 double youngsModulus,
 double poissonsRatio
);

}
