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
#include "material_utilities.h"
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
)
{
  int returnCode = 0;

  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  ScalarT* shapeTensorInvXX = shapeTensorInverseXX;
  ScalarT* shapeTensorInvXY = shapeTensorInverseXY;
  ScalarT* shapeTensorInvXZ = shapeTensorInverseXZ;
  ScalarT* shapeTensorInvYX = shapeTensorInverseYX;
  ScalarT* shapeTensorInvYY = shapeTensorInverseYY;
  ScalarT* shapeTensorInvYZ = shapeTensorInverseYZ;
  ScalarT* shapeTensorInvZX = shapeTensorInverseZX;
  ScalarT* shapeTensorInvZY = shapeTensorInverseZY;
  ScalarT* shapeTensorInvZZ = shapeTensorInverseZZ;
  ScalarT* defGradXX = deformationGradientXX;
  ScalarT* defGradXY = deformationGradientXY;
  ScalarT* defGradXZ = deformationGradientXZ;
  ScalarT* defGradYX = deformationGradientYX;
  ScalarT* defGradYY = deformationGradientYY;
  ScalarT* defGradYZ = deformationGradientYZ;
  ScalarT* defGradZX = deformationGradientZX;
  ScalarT* defGradZY = deformationGradientZY;
  ScalarT* defGradZZ = deformationGradientZZ;

  ScalarT shapeTensorXX, shapeTensorXY, shapeTensorXZ;
  ScalarT shapeTensorYX, shapeTensorYY, shapeTensorYZ;
  ScalarT shapeTensorZX, shapeTensorZY, shapeTensorZZ;
  
  ScalarT defGradFirstTermXX, defGradFirstTermXY, defGradFirstTermXZ;
  ScalarT defGradFirstTermYX, defGradFirstTermYY, defGradFirstTermYZ;
  ScalarT defGradFirstTermZX, defGradFirstTermZY, defGradFirstTermZZ;

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  double neighborVolume, omega, temp;

  // placeholder for bond damage
  double bondDamage = 0.0;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, coord+=3,
        ++shapeTensorInvXX, ++shapeTensorInvXY, ++shapeTensorInvXZ,
        ++shapeTensorInvYX, ++shapeTensorInvYY, ++shapeTensorInvYZ,
        ++shapeTensorInvZX, ++shapeTensorInvZY, ++shapeTensorInvZZ,
        ++defGradXX, ++defGradXY, ++defGradXZ,
        ++defGradYX, ++defGradYY, ++defGradYZ,
        ++defGradZX, ++defGradZY, ++defGradZZ){

    // Zero out data
    shapeTensorXX = 0.0 ; shapeTensorXY = 0.0 ; shapeTensorXZ = 0.0 ;
    shapeTensorYX = 0.0 ; shapeTensorYY = 0.0 ; shapeTensorYZ = 0.0 ;
    shapeTensorZX = 0.0 ; shapeTensorZY = 0.0 ; shapeTensorZZ = 0.0 ;
    defGradFirstTermXX = 0.0 ; defGradFirstTermXY = 0.0 ; defGradFirstTermXZ = 0.0 ;
    defGradFirstTermYX = 0.0 ; defGradFirstTermYY = 0.0 ; defGradFirstTermYZ = 0.0 ;
    defGradFirstTermZX = 0.0 ; defGradFirstTermZY = 0.0 ; defGradFirstTermZZ = 0.0 ;
    
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = volume[neighborIndex];
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, horizon);

      temp = (1.0 - bondDamage) * omega * neighborVolume;

      shapeTensorXX += temp * undeformedBondX * undeformedBondX;
      shapeTensorXY += temp * undeformedBondX * undeformedBondY;
      shapeTensorXZ += temp * undeformedBondX * undeformedBondZ;
      shapeTensorYX += temp * undeformedBondY * undeformedBondX;
      shapeTensorYY += temp * undeformedBondY * undeformedBondY;
      shapeTensorYZ += temp * undeformedBondY * undeformedBondZ;
      shapeTensorZX += temp * undeformedBondZ * undeformedBondX;
      shapeTensorZY += temp * undeformedBondZ * undeformedBondY;
      shapeTensorZZ += temp * undeformedBondZ * undeformedBondZ;

      defGradFirstTermXX += temp * deformedBondX * undeformedBondX;
      defGradFirstTermXY += temp * deformedBondX * undeformedBondY;
      defGradFirstTermXZ += temp * deformedBondX * undeformedBondZ;
      defGradFirstTermYX += temp * deformedBondY * undeformedBondX;
      defGradFirstTermYY += temp * deformedBondY * undeformedBondY;
      defGradFirstTermYZ += temp * deformedBondY * undeformedBondZ;
      defGradFirstTermZX += temp * deformedBondZ * undeformedBondX;
      defGradFirstTermZY += temp * deformedBondZ * undeformedBondY;
      defGradFirstTermZZ += temp * deformedBondZ * undeformedBondZ;
    }

    int inversionReturnCode = invert3by3Matrix(shapeTensorXX, shapeTensorXY, shapeTensorXZ,
                                               shapeTensorYX, shapeTensorYY, shapeTensorYZ,
                                               shapeTensorZX, shapeTensorZY, shapeTensorZZ,
                                               *shapeTensorInvXX, *shapeTensorInvXY, *shapeTensorInvXZ,
                                               *shapeTensorInvYX, *shapeTensorInvYY, *shapeTensorInvYZ,
                                               *shapeTensorInvZX, *shapeTensorInvZY, *shapeTensorInvZZ);

    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    *defGradXX = defGradFirstTermXX* (*shapeTensorInvXX) + defGradFirstTermXY* (*shapeTensorInvYX) + defGradFirstTermXZ* (*shapeTensorInvZX);
    *defGradXY = defGradFirstTermXX* (*shapeTensorInvXY) + defGradFirstTermXY* (*shapeTensorInvYY) + defGradFirstTermXZ* (*shapeTensorInvZY);
    *defGradXZ = defGradFirstTermXX* (*shapeTensorInvXZ) + defGradFirstTermXY* (*shapeTensorInvYZ) + defGradFirstTermXZ* (*shapeTensorInvZZ);
    *defGradYX = defGradFirstTermYX* (*shapeTensorInvXX) + defGradFirstTermYY* (*shapeTensorInvYX) + defGradFirstTermYZ* (*shapeTensorInvZX);
    *defGradYY = defGradFirstTermYX* (*shapeTensorInvXY) + defGradFirstTermYY* (*shapeTensorInvYY) + defGradFirstTermYZ* (*shapeTensorInvZY);
    *defGradYZ = defGradFirstTermYX* (*shapeTensorInvXZ) + defGradFirstTermYY* (*shapeTensorInvYZ) + defGradFirstTermYZ* (*shapeTensorInvZZ);
    *defGradZX = defGradFirstTermZX* (*shapeTensorInvXX) + defGradFirstTermZY* (*shapeTensorInvYX) + defGradFirstTermZZ* (*shapeTensorInvZX);
    *defGradZY = defGradFirstTermZX* (*shapeTensorInvXY) + defGradFirstTermZY* (*shapeTensorInvYY) + defGradFirstTermZZ* (*shapeTensorInvZY);
    *defGradZZ = defGradFirstTermZX* (*shapeTensorInvXZ) + defGradFirstTermZY* (*shapeTensorInvYZ) + defGradFirstTermZZ* (*shapeTensorInvZZ);
  }

  return returnCode;
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
  double vol, neighborVol;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  ScalarT expectedNeighborLocationX, expectedNeighborLocationY, expectedNeighborLocationZ;
  ScalarT hourglassVectorX, hourglassVectorY, hourglassVectorZ;
  ScalarT dot, magnitude;
  int neighborIndex, numNeighbors;

  const ScalarT* defGradXX = deformationGradientXX;
  const ScalarT* defGradXY = deformationGradientXY;
  const ScalarT* defGradXZ = deformationGradientXZ;
  const ScalarT* defGradYX = deformationGradientYX;
  const ScalarT* defGradYY = deformationGradientYY;
  const ScalarT* defGradYZ = deformationGradientYZ;
  const ScalarT* defGradZX = deformationGradientZX;
  const ScalarT* defGradZY = deformationGradientZY;
  const ScalarT* defGradZZ = deformationGradientZZ;

  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  ScalarT* hourglassForceDensityPtr = hourglassForceDensity;
  ScalarT* neighborHourglassForceDensityPtr;

  // placeholder for inclusion of bond damage
  double bondDamage = 0.0;

  double constant = 18.0*hourglassCoefficient*bulkModulus/(3.1415926536*horizon*horizon*horizon*horizon);

  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, coord+=3,
        ++defGradXX, ++defGradXY, ++defGradXZ,
        ++defGradYX, ++defGradYY, ++defGradYZ,
        ++defGradZX, ++defGradZY, ++defGradZZ,
        hourglassForceDensityPtr+=3){

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){
      neighborIndex = *neighborListPtr;
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      expectedNeighborLocationX = *(coord) +
        *defGradXX * undeformedBondX +
        *defGradXY * undeformedBondY +
        *defGradXZ * undeformedBondZ;
      expectedNeighborLocationY = *(coord+1) +
        *defGradYX * undeformedBondX +
        *defGradYY * undeformedBondY +
        *defGradYZ * undeformedBondZ;
      expectedNeighborLocationZ = *(coord+2) +
        *defGradZX * undeformedBondX +
        *defGradZY * undeformedBondY +
        *defGradZZ * undeformedBondZ;

      hourglassVectorX = expectedNeighborLocationX - *(neighborCoord);
      hourglassVectorY = expectedNeighborLocationY - *(neighborCoord+1);
      hourglassVectorZ = expectedNeighborLocationZ - *(neighborCoord+2);

      dot = hourglassVectorX*deformedBondX + hourglassVectorY*deformedBondY + hourglassVectorZ*deformedBondZ;
      dot *= -1.0;

      magnitude = (1.0-bondDamage) * constant * (dot/undeformedBondLength) * (1.0/deformedBondLength);

      vol = volume[iID];
      neighborVol = volume[neighborIndex];
      neighborHourglassForceDensityPtr = hourglassForceDensity + 3*neighborIndex;

      *(hourglassForceDensityPtr)   += magnitude * deformedBondX * neighborVol;
      *(hourglassForceDensityPtr+1) += magnitude * deformedBondY * neighborVol;
      *(hourglassForceDensityPtr+2) += magnitude * deformedBondZ * neighborVol;
      *(neighborHourglassForceDensityPtr)   -= magnitude * deformedBondX * vol;
      *(neighborHourglassForceDensityPtr+1) -= magnitude * deformedBondY * vol;
      *(neighborHourglassForceDensityPtr+2) -= magnitude * deformedBondZ * vol;
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

template<typename ScalarT>
int computeApproximateVelocityGradient
(
const double* volume,
const double* modelCoordinates,
const ScalarT* coordinates,
const ScalarT* velocities,
ScalarT* shapeTensorInverseXX,
ScalarT* shapeTensorInverseXY,
ScalarT* shapeTensorInverseXZ,
ScalarT* shapeTensorInverseYX,
ScalarT* shapeTensorInverseYY,
ScalarT* shapeTensorInverseYZ,
ScalarT* shapeTensorInverseZX,
ScalarT* shapeTensorInverseZY,
ScalarT* shapeTensorInverseZZ,
ScalarT* velocityGradientXX,
ScalarT* velocityGradientXY,
ScalarT* velocityGradientXZ,
ScalarT* velocityGradientYX,
ScalarT* velocityGradientYY,
ScalarT* velocityGradientYZ,
ScalarT* velocityGradientZX,
ScalarT* velocityGradientZY,
ScalarT* velocityGradientZZ,
const int* neighborhoodList,
int numPoints,
double horizon
)
{
  int returnCode = 0;

  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  ScalarT* shapeTensorInvXX = shapeTensorInverseXX;
  ScalarT* shapeTensorInvXY = shapeTensorInverseXY;
  ScalarT* shapeTensorInvXZ = shapeTensorInverseXZ;
  ScalarT* shapeTensorInvYX = shapeTensorInverseYX;
  ScalarT* shapeTensorInvYY = shapeTensorInverseYY;
  ScalarT* shapeTensorInvYZ = shapeTensorInverseYZ;
  ScalarT* shapeTensorInvZX = shapeTensorInverseZX;
  ScalarT* shapeTensorInvZY = shapeTensorInverseZY;
  ScalarT* shapeTensorInvZZ = shapeTensorInverseZZ;
  ScalarT* velGradXX = velocityGradientXX;
  ScalarT* velGradXY = velocityGradientXY;
  ScalarT* velGradXZ = velocityGradientXZ;
  ScalarT* velGradYX = velocityGradientYX;
  ScalarT* velGradYY = velocityGradientYY;
  ScalarT* velGradYZ = velocityGradientYZ;
  ScalarT* velGradZX = velocityGradientZX;
  ScalarT* velGradZY = velocityGradientZY;
  ScalarT* velGradZZ = velocityGradientZZ;

  ScalarT velGradFirstTermXX, velGradFirstTermXY, velGradFirstTermXZ;
  ScalarT velGradFirstTermYX, velGradFirstTermYY, velGradFirstTermYZ;
  ScalarT velGradFirstTermZX, velGradFirstTermZY, velGradFirstTermZZ;

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  double neighborVolume, omega, temp;
  ScalarT velStateX, velStateY, velStateZ;

  // placeholder for bond damage
  double bondDamage = 0.0;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, coord+=3,
        ++shapeTensorInvXX, ++shapeTensorInvXY, ++shapeTensorInvXZ,
        ++shapeTensorInvYX, ++shapeTensorInvYY, ++shapeTensorInvYZ,
        ++shapeTensorInvZX, ++shapeTensorInvZY, ++shapeTensorInvZZ,
        ++velGradXX, ++velGradXY, ++velGradXZ,
        ++velGradYX, ++velGradYY, ++velGradYZ,
        ++velGradZX, ++velGradZY, ++velGradZZ){

    // Zero out data
    velGradFirstTermXX = 0.0 ; velGradFirstTermXY = 0.0 ; velGradFirstTermXZ = 0.0 ;
    velGradFirstTermYX = 0.0 ; velGradFirstTermYY = 0.0 ; velGradFirstTermYZ = 0.0 ;
    velGradFirstTermZX = 0.0 ; velGradFirstTermZY = 0.0 ; velGradFirstTermZZ = 0.0 ;
    
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = volume[neighborIndex];
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborCoord = coordinates + 3*neighborIndex;
      neighborVel = velocities + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);

      // The velState is the relative difference in velocities of the nodes at
      // each end of a bond. i.e., v_j - v_i
      velStateX = *(neighborVel) - *(vel);
      velStateY = *(neighborVel+1) - *(vel+1);
      velStateZ = *(neighborVel+2) - *(vel+2);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, horizon);

      temp = (1.0 - bondDamage) * omega * neighborVolume;

      velGradFirstTermXX += temp * velStateX * undeformedBondX;
      velGradFirstTermXY += temp * velStateX * undeformedBondY;
      velGradFirstTermXZ += temp * velStateX * undeformedBondZ;
      velGradFirstTermYX += temp * velStateY * undeformedBondX;
      velGradFirstTermYY += temp * velStateY * undeformedBondY;
      velGradFirstTermYZ += temp * velStateY * undeformedBondZ;
      velGradFirstTermZX += temp * velStateZ * undeformedBondX;
      velGradFirstTermZY += temp * velStateZ * undeformedBondY;
      velGradFirstTermZZ += temp * velStateZ * undeformedBondZ;
    }


    *velGradXX = velGradFirstTermXX* (*shapeTensorInvXX) + velGradFirstTermXY* (*shapeTensorInvYX) + velGradFirstTermXZ* (*shapeTensorInvZX);
    *velGradXY = velGradFirstTermXX* (*shapeTensorInvXY) + velGradFirstTermXY* (*shapeTensorInvYY) + velGradFirstTermXZ* (*shapeTensorInvZY);
    *velGradXZ = velGradFirstTermXX* (*shapeTensorInvXZ) + velGradFirstTermXY* (*shapeTensorInvYZ) + velGradFirstTermXZ* (*shapeTensorInvZZ);
    *velGradYX = velGradFirstTermYX* (*shapeTensorInvXX) + velGradFirstTermYY* (*shapeTensorInvYX) + velGradFirstTermYZ* (*shapeTensorInvZX);
    *velGradYY = velGradFirstTermYX* (*shapeTensorInvXY) + velGradFirstTermYY* (*shapeTensorInvYY) + velGradFirstTermYZ* (*shapeTensorInvZY);
    *velGradYZ = velGradFirstTermYX* (*shapeTensorInvXZ) + velGradFirstTermYY* (*shapeTensorInvYZ) + velGradFirstTermYZ* (*shapeTensorInvZZ);
    *velGradZX = velGradFirstTermZX* (*shapeTensorInvXX) + velGradFirstTermZY* (*shapeTensorInvYX) + velGradFirstTermZZ* (*shapeTensorInvZX);
    *velGradZY = velGradFirstTermZX* (*shapeTensorInvXY) + velGradFirstTermZY* (*shapeTensorInvYY) + velGradFirstTermZZ* (*shapeTensorInvZY);
    *velGradZZ = velGradFirstTermZX* (*shapeTensorInvXZ) + velGradFirstTermZY* (*shapeTensorInvYZ) + velGradFirstTermZZ* (*shapeTensorInvZZ);
  }

  // Placeholder for errors due to node being too damaged
  return returnCode;
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

template int computeShapeTensorInverseAndApproximateDeformationGradient<double>
(
const double* volume,
const double* modelCoordinates,
const double* coordinates,
double* shapeTensorInverseXX,
double* shapeTensorInverseXY,
double* shapeTensorInverseXZ,
double* shapeTensorInverseYX,
double* shapeTensorInverseYY,
double* shapeTensorInverseYZ,
double* shapeTensorInverseZX,
double* shapeTensorInverseZY,
double* shapeTensorInverseZZ,
double* deformationGradientXX,
double* deformationGradientXY,
double* deformationGradientXZ,
double* deformationGradientYX,
double* deformationGradientYY,
double* deformationGradientYZ,
double* deformationGradientZX,
double* deformationGradientZY,
double* deformationGradientZZ,
const int* neighborhoodList,
int numPoints,
double horizon
);

template int computeApproximateVelocityGradient<double>
(
const double* volume,
const double* modelCoordinates,
const double* coordinates,
const double* velocities,
double* shapeTensorInverseXX,
double* shapeTensorInverseXY,
double* shapeTensorInverseXZ,
double* shapeTensorInverseYX,
double* shapeTensorInverseYY,
double* shapeTensorInverseYZ,
double* shapeTensorInverseZX,
double* shapeTensorInverseZY,
double* shapeTensorInverseZZ,
double* velocityGradientXX,
double* velocityGradientXY,
double* velocityGradientXZ,
double* velocityGradientYX,
double* velocityGradientYY,
double* velocityGradientYZ,
double* velocityGradientZX,
double* velocityGradientZY,
double* velocityGradientZZ,
const int* neighborhoodList,
int numPoints,
double horizon
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

template int computeShapeTensorInverseAndApproximateDeformationGradient<Sacado::Fad::DFad<double> >
(
const double* volume,
const double* modelCoordinates,
const Sacado::Fad::DFad<double>* coordinates,
Sacado::Fad::DFad<double>* shapeTensorInverseXX,
Sacado::Fad::DFad<double>* shapeTensorInverseXY,
Sacado::Fad::DFad<double>* shapeTensorInverseXZ,
Sacado::Fad::DFad<double>* shapeTensorInverseYX,
Sacado::Fad::DFad<double>* shapeTensorInverseYY,
Sacado::Fad::DFad<double>* shapeTensorInverseYZ,
Sacado::Fad::DFad<double>* shapeTensorInverseZX,
Sacado::Fad::DFad<double>* shapeTensorInverseZY,
Sacado::Fad::DFad<double>* shapeTensorInverseZZ,
Sacado::Fad::DFad<double>* deformationGradientXX,
Sacado::Fad::DFad<double>* deformationGradientXY,
Sacado::Fad::DFad<double>* deformationGradientXZ,
Sacado::Fad::DFad<double>* deformationGradientYX,
Sacado::Fad::DFad<double>* deformationGradientYY,
Sacado::Fad::DFad<double>* deformationGradientYZ,
Sacado::Fad::DFad<double>* deformationGradientZX,
Sacado::Fad::DFad<double>* deformationGradientZY,
Sacado::Fad::DFad<double>* deformationGradientZZ,
const int* neighborhoodList,
int numPoints,
double horizon
);

template int computeApproximateVelocityGradient<Sacado::Fad::DFad<double> >
(
const double* volume,
const double* modelCoordinates,
const Sacado::Fad::DFad<double>* coordinates,
const Sacado::Fad::DFad<double>* velocities,
Sacado::Fad::DFad<double>* shapeTensorInverseXX,
Sacado::Fad::DFad<double>* shapeTensorInverseXY,
Sacado::Fad::DFad<double>* shapeTensorInverseXZ,
Sacado::Fad::DFad<double>* shapeTensorInverseYX,
Sacado::Fad::DFad<double>* shapeTensorInverseYY,
Sacado::Fad::DFad<double>* shapeTensorInverseYZ,
Sacado::Fad::DFad<double>* shapeTensorInverseZX,
Sacado::Fad::DFad<double>* shapeTensorInverseZY,
Sacado::Fad::DFad<double>* shapeTensorInverseZZ,
Sacado::Fad::DFad<double>* velocityGradientXX,
Sacado::Fad::DFad<double>* velocityGradientXY,
Sacado::Fad::DFad<double>* velocityGradientXZ,
Sacado::Fad::DFad<double>* velocityGradientYX,
Sacado::Fad::DFad<double>* velocityGradientYY,
Sacado::Fad::DFad<double>* velocityGradientYZ,
Sacado::Fad::DFad<double>* velocityGradientZX,
Sacado::Fad::DFad<double>* velocityGradientZY,
Sacado::Fad::DFad<double>* velocityGradientZZ,
const int* neighborhoodList,
int numPoints,
double horizon
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
