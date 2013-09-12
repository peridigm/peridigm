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

// Compute result = Rtranspose * A * R
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
)
{
  resultXX = rXX*(aXX*rXX + aYX*rYX + aZX*rZX) + rYX*(aXY*rXX + aYY*rYX + aZY*rZX) + 
   rZX*(aXZ*rXX + aYZ*rYX + aZZ*rZX);
  resultXY = aXX*rXX*rXY + aYX*rXY*rYX + aXY*rXX*rYY + aYY*rYX*rYY + aZX*rXY*rZX + 
   aZY*rYY*rZX + aXZ*rXX*rZY + aYZ*rYX*rZY + aZZ*rZX*rZY;
  resultXZ = aXX*rXX*rXZ + aYX*rXZ*rYX + aXY*rXX*rYZ + aYY*rYX*rYZ + aZX*rXZ*rZX + 
   aZY*rYZ*rZX + aXZ*rXX*rZZ + aYZ*rYX*rZZ + aZZ*rZX*rZZ;
  resultYX = aXX*rXX*rXY + aXY*rXY*rYX + aYX*rXX*rYY + aYY*rYX*rYY + aXZ*rXY*rZX + 
   aYZ*rYY*rZX + aZX*rXX*rZY + aZY*rYX*rZY + aZZ*rZX*rZY;
  resultYY = rXY*(aXX*rXY + aYX*rYY + aZX*rZY) + rYY*(aXY*rXY + aYY*rYY + aZY*rZY) + 
   rZY*(aXZ*rXY + aYZ*rYY + aZZ*rZY);
  resultYZ = aXX*rXY*rXZ + aYX*rXZ*rYY + aXY*rXY*rYZ + aYY*rYY*rYZ + aZX*rXZ*rZY + 
   aZY*rYZ*rZY + aXZ*rXY*rZZ + aYZ*rYY*rZZ + aZZ*rZY*rZZ;
  resultZX = aXX*rXX*rXZ + aXY*rXZ*rYX + aYX*rXX*rYZ + aYY*rYX*rYZ + aXZ*rXZ*rZX + 
   aYZ*rYZ*rZX + aZX*rXX*rZZ + aZY*rYX*rZZ + aZZ*rZX*rZZ;
  resultZY = aXX*rXY*rXZ + aXY*rXZ*rYY + aYX*rXY*rYZ + aYY*rYY*rYZ + aXZ*rXZ*rZY + 
   aYZ*rYZ*rZY + aZX*rXY*rZZ + aZY*rYY*rZZ + aZZ*rZY*rZZ;
  resultZZ = rXZ*(aXX*rXZ + aYX*rYZ + aZX*rZZ) + rYZ*(aXY*rXZ + aYY*rYZ + aZY*rZZ) + 
   rZZ*(aXZ*rXZ + aYZ*rYZ + aZZ*rZZ);
}

// Compute result = R * A * Rtranspose
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
)
{
  resultXX = rXX*(aXX*rXX + aYX*rXY + aZX*rXZ) + rXY*(aXY*rXX + aYY*rXY + aZY*rXZ) + 
   rXZ*(aXZ*rXX + aYZ*rXY + aZZ*rXZ);
  resultXY = aXX*rXX*rYX + aYX*rXY*rYX + aZX*rXZ*rYX + aXY*rXX*rYY + aYY*rXY*rYY + 
   aZY*rXZ*rYY + aXZ*rXX*rYZ + aYZ*rXY*rYZ + aZZ*rXZ*rYZ;
  resultXZ = aXX*rXX*rZX + aYX*rXY*rZX + aZX*rXZ*rZX + aXY*rXX*rZY + aYY*rXY*rZY + 
   aZY*rXZ*rZY + aXZ*rXX*rZZ + aYZ*rXY*rZZ + aZZ*rXZ*rZZ;
  resultYX = aXX*rXX*rYX + aXY*rXY*rYX + aXZ*rXZ*rYX + aYX*rXX*rYY + aYY*rXY*rYY + 
   aYZ*rXZ*rYY + aZX*rXX*rYZ + aZY*rXY*rYZ + aZZ*rXZ*rYZ;
  resultYY = rYX*(aXX*rYX + aYX*rYY + aZX*rYZ) + rYY*(aXY*rYX + aYY*rYY + aZY*rYZ) + 
   rYZ*(aXZ*rYX + aYZ*rYY + aZZ*rYZ);
  resultYZ = aXX*rYX*rZX + aYX*rYY*rZX + aZX*rYZ*rZX + aXY*rYX*rZY + aYY*rYY*rZY + 
   aZY*rYZ*rZY + aXZ*rYX*rZZ + aYZ*rYY*rZZ + aZZ*rYZ*rZZ;
  resultZX = aXX*rXX*rZX + aXY*rXY*rZX + aXZ*rXZ*rZX + aYX*rXX*rZY + aYY*rXY*rZY + 
   aYZ*rXZ*rZY + aZX*rXX*rZZ + aZY*rXY*rZZ + aZZ*rXZ*rZZ;
  resultZY = aXX*rYX*rZX + aXY*rYY*rZX + aXZ*rYZ*rZX + aYX*rYX*rZY + aYY*rYY*rZY + 
   aYZ*rYZ*rZY + aZX*rYX*rZZ + aZY*rYY*rZZ + aZZ*rYZ*rZZ;
  resultZZ = rZX*(aXX*rZX + aYX*rZY + aZX*rZZ) + rZY*(aXY*rZX + aYY*rZY + aZY*rZZ) + 
   rZZ*(aXZ*rZX + aYZ*rZY + aZZ*rZZ);
}

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
)
{
  //compute result = alpha * a + beta * b
  resultXX = alpha * aXX + beta * bXX;
  resultXY = alpha * aXY + beta * bXY;
  resultXZ = alpha * aXZ + beta * bXZ;
  resultYX = alpha * aYX + beta * bYX;
  resultYY = alpha * aYY + beta * bYY;
  resultYZ = alpha * aYZ + beta * bYZ;
  resultZX = alpha * aZX + beta * bZX;
  resultZY = alpha * aZY + beta * bZY;
  resultZZ = alpha * aZZ + beta * bZZ;
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

//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
template<typename ScalarT>
int computeUnrotatedRateOfDeformationAndRotationTensor(
const double* volume,
const double* modelCoordinates,
const ScalarT* coordinates,
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
ScalarT* leftStretchTensorXX,
ScalarT* leftStretchTensorXY,
ScalarT* leftStretchTensorXZ,
ScalarT* leftStretchTensorYX,
ScalarT* leftStretchTensorYY,
ScalarT* leftStretchTensorYZ,
ScalarT* leftStretchTensorZX,
ScalarT* leftStretchTensorZY,
ScalarT* leftStretchTensorZZ,
ScalarT* rotationTensorXX,
ScalarT* rotationTensorXY,
ScalarT* rotationTensorXZ,
ScalarT* rotationTensorYX,
ScalarT* rotationTensorYY,
ScalarT* rotationTensorYZ,
ScalarT* rotationTensorZX,
ScalarT* rotationTensorZY,
ScalarT* rotationTensorZZ,
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
)
{
  int returnCode = 0;

  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;

  const ScalarT* defGradXX = deformationGradientXX;
  const ScalarT* defGradXY = deformationGradientXY;
  const ScalarT* defGradXZ = deformationGradientXZ;
  const ScalarT* defGradYX = deformationGradientYX;
  const ScalarT* defGradYY = deformationGradientYY;
  const ScalarT* defGradYZ = deformationGradientYZ;
  const ScalarT* defGradZX = deformationGradientZX;
  const ScalarT* defGradZY = deformationGradientZY;
  const ScalarT* defGradZZ = deformationGradientZZ;

  const ScalarT* shapeTensorInvXX = shapeTensorInverseXX;
  const ScalarT* shapeTensorInvXY = shapeTensorInverseXY;
  const ScalarT* shapeTensorInvXZ = shapeTensorInverseXZ;
  const ScalarT* shapeTensorInvYX = shapeTensorInverseYX;
  const ScalarT* shapeTensorInvYY = shapeTensorInverseYY;
  const ScalarT* shapeTensorInvYZ = shapeTensorInverseYZ;
  const ScalarT* shapeTensorInvZX = shapeTensorInverseZX;
  const ScalarT* shapeTensorInvZY = shapeTensorInverseZY;
  const ScalarT* shapeTensorInvZZ = shapeTensorInverseZZ;

  ScalarT* leftStretchXX = leftStretchTensorXX;
  ScalarT* leftStretchXY = leftStretchTensorXY;
  ScalarT* leftStretchXZ = leftStretchTensorXZ;
  ScalarT* leftStretchYX = leftStretchTensorYX;
  ScalarT* leftStretchYY = leftStretchTensorYY;
  ScalarT* leftStretchYZ = leftStretchTensorYZ;
  ScalarT* leftStretchZX = leftStretchTensorZX;
  ScalarT* leftStretchZY = leftStretchTensorZY;
  ScalarT* leftStretchZZ = leftStretchTensorZZ;

  ScalarT* rotTensorXX = rotationTensorXX;
  ScalarT* rotTensorXY = rotationTensorXY;
  ScalarT* rotTensorXZ = rotationTensorXZ;
  ScalarT* rotTensorYX = rotationTensorYX;
  ScalarT* rotTensorYY = rotationTensorYY;
  ScalarT* rotTensorYZ = rotationTensorYZ;
  ScalarT* rotTensorZX = rotationTensorZX;
  ScalarT* rotTensorZY = rotationTensorZY;
  ScalarT* rotTensorZZ = rotationTensorZZ;

  ScalarT* unrotRateOfDefXX = unrotatedRateOfDeformationXX;
  ScalarT* unrotRateOfDefXY = unrotatedRateOfDeformationXY;
  ScalarT* unrotRateOfDefXZ = unrotatedRateOfDeformationXZ;
  ScalarT* unrotRateOfDefYX = unrotatedRateOfDeformationYX;
  ScalarT* unrotRateOfDefYY = unrotatedRateOfDeformationYY;
  ScalarT* unrotRateOfDefYZ = unrotatedRateOfDeformationYZ;
  ScalarT* unrotRateOfDefZX = unrotatedRateOfDeformationZX;
  ScalarT* unrotRateOfDefZY = unrotatedRateOfDeformationZY;
  ScalarT* unrotRateOfDefZZ = unrotatedRateOfDeformationZZ;

  ScalarT FdotFirstTermXX, FdotFirstTermXY, FdotFirstTermXZ;
  ScalarT FdotFirstTermYX, FdotFirstTermYY, FdotFirstTermYZ;
  ScalarT FdotFirstTermZX, FdotFirstTermZY, FdotFirstTermZZ;

  ScalarT FdotXX, FdotXY, FdotXZ;
  ScalarT FdotYX, FdotYY, FdotYZ;
  ScalarT FdotZX, FdotZY, FdotZZ;

  ScalarT eulerianVelGradXX, eulerianVelGradXY, eulerianVelGradXZ;
  ScalarT eulerianVelGradYX, eulerianVelGradYY, eulerianVelGradYZ;
  ScalarT eulerianVelGradZX, eulerianVelGradZY, eulerianVelGradZZ;

  ScalarT defGradInvXX, defGradInvXY, defGradInvXZ;
  ScalarT defGradInvYX, defGradInvYY, defGradInvYZ;
  ScalarT defGradInvZX, defGradInvZY, defGradInvZZ;

  ScalarT rateOfDefXX, rateOfDefXY, rateOfDefXZ;
  ScalarT rateOfDefYX, rateOfDefYY, rateOfDefYZ;
  ScalarT rateOfDefZX, rateOfDefZY, rateOfDefZZ;

  ScalarT spinXX, spinXY, spinXZ;
  ScalarT spinYX, spinYY, spinYZ;
  ScalarT spinZX, spinZY, spinZZ;

  ScalarT dvXX, dvXY, dvXZ;
  ScalarT dvYX, dvYY, dvYZ;
  ScalarT dvZX, dvZY, dvZZ;

  ScalarT tempXX, tempXY, tempXZ;
  ScalarT tempYX, tempYY, tempYZ;
  ScalarT tempZX, tempZY, tempZZ;

  ScalarT tempInvXX, tempInvXY, tempInvXZ;
  ScalarT tempInvYX, tempInvYY, tempInvYZ;
  ScalarT tempInvZX, tempInvZY, tempInvZZ;

  ScalarT tempAXX, tempAXY, tempAXZ;
  ScalarT tempAYX, tempAYY, tempAYZ;
  ScalarT tempAZX, tempAZY, tempAZZ;
  
  ScalarT tempBXX, tempBXY, tempBXZ;
  ScalarT tempBYX, tempBYY, tempBYZ;
  ScalarT tempBZX, tempBZY, tempBZZ;

  ScalarT tempBInvXX, tempBInvXY, tempBInvXZ;
  ScalarT tempBInvYX, tempBInvYY, tempBInvYZ;
  ScalarT tempBInvZX, tempBInvZY, tempBInvZZ;

  ScalarT omegaTensorXX, omegaTensorXY, omegaTensorXZ;
  ScalarT omegaTensorYX, omegaTensorYY, omegaTensorYZ;
  ScalarT omegaTensorZX, omegaTensorZY, omegaTensorZZ;

  ScalarT rotTensorOldXX, rotTensorOldXY, rotTensorOldXZ;
  ScalarT rotTensorOldYX, rotTensorOldYY, rotTensorOldYZ;
  ScalarT rotTensorOldZX, rotTensorOldZY, rotTensorOldZZ;

  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
 
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT traceV, scaleFactor;

  ScalarT identityMatrixXX, identityMatrixXY, identityMatrixXZ;
  ScalarT identityMatrixYX, identityMatrixYY, identityMatrixYZ;
  ScalarT identityMatrixZX, identityMatrixZY, identityMatrixZZ;

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double neighborVolume, omega, scalarTemp; 

  // placeholder for bond damage
  double bondDamage = 0.0;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, coord+=3, vel+=3,
        ++shapeTensorInvXX, ++shapeTensorInvXY, ++shapeTensorInvXZ,
        ++shapeTensorInvYX, ++shapeTensorInvYY, ++shapeTensorInvYZ,
        ++shapeTensorInvZX, ++shapeTensorInvZY, ++shapeTensorInvZZ,
        ++rotTensorXX, ++rotTensorXY, ++rotTensorXZ,
        ++rotTensorYX, ++rotTensorYY, ++rotTensorYZ,
        ++rotTensorZX, ++rotTensorZY, ++rotTensorZZ,
        ++leftStretchXX, ++leftStretchXY, ++leftStretchXZ,
        ++leftStretchYX, ++leftStretchYY, ++leftStretchYZ,
        ++leftStretchZX, ++leftStretchZY, ++leftStretchZZ,
        ++unrotRateOfDefXX, ++unrotRateOfDefXY, ++unrotRateOfDefXZ,
        ++unrotRateOfDefYX, ++unrotRateOfDefYY, ++unrotRateOfDefYZ,
        ++unrotRateOfDefZX, ++unrotRateOfDefZY, ++unrotRateOfDefZZ,
        ++defGradXX, ++defGradXY, ++defGradXZ,
        ++defGradYX, ++defGradYY, ++defGradYZ,
        ++defGradZX, ++defGradZY, ++defGradZZ){

    // Initialize data
    FdotFirstTermXX = 0.0, FdotFirstTermXY = 0.0, FdotFirstTermXZ = 0.0;
    FdotFirstTermYX = 0.0, FdotFirstTermYY = 0.0, FdotFirstTermYZ = 0.0;
    FdotFirstTermZX = 0.0, FdotFirstTermZY = 0.0, FdotFirstTermZZ = 0.0;
    FdotXX = 0.0, FdotXY = 0.0, FdotXZ = 0.0;
    FdotYX = 0.0, FdotYY = 0.0, FdotYZ = 0.0;
    FdotZX = 0.0, FdotZY = 0.0, FdotZZ = 0.0;
    eulerianVelGradXX = 0.0, eulerianVelGradXY = 0.0, eulerianVelGradXZ = 0.0;
    eulerianVelGradYX = 0.0, eulerianVelGradYY = 0.0, eulerianVelGradYZ = 0.0;
    eulerianVelGradZX = 0.0, eulerianVelGradZY = 0.0, eulerianVelGradZZ = 0.0;
    defGradInvXX = 0.0, defGradInvXY = 0.0, defGradInvXZ = 0.0;
    defGradInvYX = 0.0, defGradInvYY = 0.0, defGradInvYZ = 0.0;
    defGradInvZX = 0.0, defGradInvZY = 0.0, defGradInvZZ = 0.0;
    rateOfDefXX = 0.0, rateOfDefXY = 0.0, rateOfDefXZ = 0.0;
    rateOfDefYX = 0.0, rateOfDefYY = 0.0, rateOfDefYZ = 0.0;
    rateOfDefZX = 0.0, rateOfDefZY = 0.0, rateOfDefZZ = 0.0;
    spinXX = 0.0, spinXY = 0.0, spinXZ = 0.0;
    spinYX = 0.0, spinYY = 0.0, spinYZ = 0.0;
    spinZX = 0.0, spinZY = 0.0, spinZZ = 0.0;
    dvXX = 0.0, dvXY = 0.0, dvXZ = 0.0;
    dvYX = 0.0, dvYY = 0.0, dvYZ = 0.0;
    dvZX = 0.0, dvZY = 0.0, dvZZ = 0.0;
    tempXX = 0.0, tempXY = 0.0, tempXZ = 0.0;
    tempYX = 0.0, tempYY = 0.0, tempYZ = 0.0;
    tempZX = 0.0, tempZY = 0.0, tempZZ = 0.0;
    tempInvXX = 0.0, tempInvXY = 0.0, tempInvXZ = 0.0;
    tempInvYX = 0.0, tempInvYY = 0.0, tempInvYZ = 0.0;
    tempInvZX = 0.0, tempInvZY = 0.0, tempInvZZ = 0.0;
    tempAXX = 0.0, tempAXY = 0.0, tempAXZ = 0.0;
    tempAYX = 0.0, tempAYY = 0.0, tempAYZ = 0.0;
    tempAZX = 0.0, tempAZY = 0.0, tempAZZ = 0.0;
    tempBXX = 0.0, tempBXY = 0.0, tempBXZ = 0.0;
    tempBYX = 0.0, tempBYY = 0.0, tempBYZ = 0.0;
    tempBZX = 0.0, tempBZY = 0.0, tempBZZ = 0.0;
    tempBInvXX = 0.0, tempBInvXY = 0.0, tempBInvXZ = 0.0;
    tempBInvYX = 0.0, tempBInvYY = 0.0, tempBInvYZ = 0.0;
    tempBInvZX = 0.0, tempBInvZY = 0.0, tempBInvZZ = 0.0;
    omegaTensorXX = 0.0, omegaTensorXY = 0.0, omegaTensorXZ = 0.0;
    omegaTensorYX = 0.0, omegaTensorYY = 0.0, omegaTensorYZ = 0.0;
    omegaTensorZX = 0.0, omegaTensorZY = 0.0, omegaTensorZZ = 0.0;
    omegaX = 0.0, omegaY = 0.0, omegaZ = 0.0;
    zX = 0.0, zY = 0.0, zZ = 0.0;
    wX = 0.0, wY = 0.0, wZ = 0.0;
    identityMatrixXX = 1.0, identityMatrixXY = 0.0, identityMatrixXZ = 0.0;
    identityMatrixYX = 0.0, identityMatrixYY = 1.0, identityMatrixYZ = 0.0;
    identityMatrixZX = 0.0, identityMatrixZY = 0.0, identityMatrixZZ = 1.0;
    deformedBondX = 0.0, deformedBondY = 0.0, deformedBondZ = 0.0;
    velStateX = 0.0, velStateY = 0.0, velStateZ = 0.0;
    undeformedBondX = 0.0, undeformedBondY = 0.0, undeformedBondZ = 0.0; 
    undeformedBondLength = 0.0;
    neighborVolume = 0.0, omega = 0.0, scalarTemp = 0.0, scaleFactor = 0.0;
    traceV = 0.0;
    rotTensorOldXX = *rotTensorXX, rotTensorOldXY = *rotTensorXY, rotTensorOldXZ = *rotTensorXZ;
    rotTensorOldYX = *rotTensorYX, rotTensorOldYY = *rotTensorYY, rotTensorOldYZ = *rotTensorYZ;
    rotTensorOldZX = *rotTensorZX, rotTensorOldZY = *rotTensorZY, rotTensorOldZZ = *rotTensorZZ;
    
    //Compute Fdot
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

      scalarTemp = (1.0 - bondDamage) * omega * neighborVolume;

      FdotFirstTermXX += scalarTemp * velStateX * undeformedBondX;
      FdotFirstTermXY += scalarTemp * velStateX * undeformedBondY;
      FdotFirstTermXZ += scalarTemp * velStateX * undeformedBondZ;
      FdotFirstTermYX += scalarTemp * velStateY * undeformedBondX;
      FdotFirstTermYY += scalarTemp * velStateY * undeformedBondY;
      FdotFirstTermYZ += scalarTemp * velStateY * undeformedBondZ;
      FdotFirstTermZX += scalarTemp * velStateZ * undeformedBondX;
      FdotFirstTermZY += scalarTemp * velStateZ * undeformedBondY;
      FdotFirstTermZZ += scalarTemp * velStateZ * undeformedBondZ;
    }

    FdotXX = FdotFirstTermXX* (*shapeTensorInvXX) + FdotFirstTermXY* (*shapeTensorInvYX) + FdotFirstTermXZ* (*shapeTensorInvZX);
    FdotXY = FdotFirstTermXX* (*shapeTensorInvXY) + FdotFirstTermXY* (*shapeTensorInvYY) + FdotFirstTermXZ* (*shapeTensorInvZY);
    FdotXZ = FdotFirstTermXX* (*shapeTensorInvXZ) + FdotFirstTermXY* (*shapeTensorInvYZ) + FdotFirstTermXZ* (*shapeTensorInvZZ);
    FdotYX = FdotFirstTermYX* (*shapeTensorInvXX) + FdotFirstTermYY* (*shapeTensorInvYX) + FdotFirstTermYZ* (*shapeTensorInvZX);
    FdotYY = FdotFirstTermYX* (*shapeTensorInvXY) + FdotFirstTermYY* (*shapeTensorInvYY) + FdotFirstTermYZ* (*shapeTensorInvZY);
    FdotYZ = FdotFirstTermYX* (*shapeTensorInvXZ) + FdotFirstTermYY* (*shapeTensorInvYZ) + FdotFirstTermYZ* (*shapeTensorInvZZ);
    FdotZX = FdotFirstTermZX* (*shapeTensorInvXX) + FdotFirstTermZY* (*shapeTensorInvYX) + FdotFirstTermZZ* (*shapeTensorInvZX);
    FdotZY = FdotFirstTermZX* (*shapeTensorInvXY) + FdotFirstTermZY* (*shapeTensorInvYY) + FdotFirstTermZZ* (*shapeTensorInvZY);
    FdotZZ = FdotFirstTermZX* (*shapeTensorInvXZ) + FdotFirstTermZY* (*shapeTensorInvYZ) + FdotFirstTermZZ* (*shapeTensorInvZZ);


    // Compute the inverse of the deformation gradient, Finv
    int inversionReturnCode = invert3by3Matrix(*defGradXX,*defGradXY,*defGradXZ,
                                               *defGradYX,*defGradYY,*defGradYZ,
                                               *defGradZX,*defGradZY,*defGradZZ,
                                               defGradInvXX, defGradInvXY, defGradInvXZ,
                                               defGradInvYX, defGradInvYY, defGradInvYZ,
                                               defGradInvZX, defGradInvZY, defGradInvZZ);
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    // Compute the Eulerian velocity gradient L = Fdot * Finv
    MatrixMultiply(FdotXX, FdotXY, FdotXZ,
                   FdotYX, FdotYY, FdotYZ,
                   FdotZX, FdotZY, FdotZZ,
                   defGradInvXX, defGradInvXY, defGradInvXZ,
                   defGradInvYX, defGradInvYY, defGradInvYZ,
                   defGradInvZX, defGradInvZY, defGradInvZZ,
                   eulerianVelGradXX, eulerianVelGradXY, eulerianVelGradXZ,
                   eulerianVelGradYX, eulerianVelGradYY, eulerianVelGradYZ,
                   eulerianVelGradZX, eulerianVelGradZY, eulerianVelGradZZ);


    // Compute rate-of-deformation tensor, D
    MatrixUpdate(ScalarT(0.5), ScalarT(0.5),
                 eulerianVelGradXX, eulerianVelGradXY, eulerianVelGradXZ,
                 eulerianVelGradYX, eulerianVelGradYY, eulerianVelGradYZ,
                 eulerianVelGradZX, eulerianVelGradZY, eulerianVelGradZZ,
                 eulerianVelGradXX, eulerianVelGradYX, eulerianVelGradZX,
                 eulerianVelGradYX, eulerianVelGradYY, eulerianVelGradZY,
                 eulerianVelGradXZ, eulerianVelGradYZ, eulerianVelGradZZ,
                 rateOfDefXX, rateOfDefXY, rateOfDefXZ,
                 rateOfDefYX, rateOfDefYY, rateOfDefYZ,
                 rateOfDefZX, rateOfDefZY, rateOfDefZZ);

    // Compute spin tensor, W
    MatrixUpdate(ScalarT(0.5), ScalarT(-0.5),
                 eulerianVelGradXX, eulerianVelGradXY, eulerianVelGradXZ,
                 eulerianVelGradYX, eulerianVelGradYY, eulerianVelGradYZ,
                 eulerianVelGradZX, eulerianVelGradZY, eulerianVelGradZZ,
                 eulerianVelGradXX, eulerianVelGradYX, eulerianVelGradZX,
                 eulerianVelGradYX, eulerianVelGradYY, eulerianVelGradZY,
                 eulerianVelGradXZ, eulerianVelGradYZ, eulerianVelGradZZ,
                 spinXX, spinXY, spinXZ,
                 spinYX, spinYY, spinYZ,
                 spinZX, spinZY, spinZZ);

    //Following Flanagan & Taylor (T&F) find the matrix product D*V
    MatrixMultiply(rateOfDefXX, rateOfDefXY, rateOfDefXZ,
                   rateOfDefYX, rateOfDefYY, rateOfDefYZ,
                   rateOfDefZX, rateOfDefZY, rateOfDefZZ,
                   *leftStretchXX, *leftStretchXY, *leftStretchXZ,
                   *leftStretchYX, *leftStretchYY, *leftStretchYZ,
                   *leftStretchZX, *leftStretchZY, *leftStretchZZ,
                   dvXX, dvXY, dvXZ,
                   dvYX, dvYY, dvYZ,
                   dvZX, dvZY, dvZZ);

    //Find the vector z, the vector dual of DV - VD = 2*dual(DV)
    zX = dvZY - dvYZ;
    zY = -dvZX + dvXZ;
    zZ = dvYX - dvXY;

    //Find the vector w, the vector dual of W
    wX = -0.5 * ( spinYZ - spinZY);
    wY = -0.5 * ( spinXZ - spinZX);
    wZ = -0.5 * ( spinXY - spinYX);

    //Find trace(V)
    traceV = (*leftStretchXX) + (*leftStretchYY) + (*leftStretchZZ);

    // Compute trace(V) I - V store in temp
    MatrixUpdate(traceV, ScalarT(-1.0),
                 identityMatrixXX, identityMatrixXY, identityMatrixXZ,
                 identityMatrixYX, identityMatrixYY, identityMatrixYZ,
                 identityMatrixZX, identityMatrixZY, identityMatrixZZ,
                 *leftStretchXX, *leftStretchYX, *leftStretchZX,
                 *leftStretchYX, *leftStretchYY, *leftStretchZY,
                 *leftStretchXZ, *leftStretchYZ, *leftStretchZZ,
                 tempXX, tempXY, tempXZ,
                 tempYX, tempYY, tempYZ,
                 tempZX, tempZY, tempZZ);

    // Compute the inverse of the temp matrix
    inversionReturnCode = invert3by3Matrix(tempXX, tempXY, tempXZ,
                                           tempYX, tempYY, tempYZ,
                                           tempZX, tempZY, tempZZ,
                                           tempInvXX, tempInvXY, tempInvXZ,
                                           tempInvYX, tempInvYY, tempInvYZ,
                                           tempInvZX, tempInvZY, tempInvZZ);
    // Placeholder for more sophisticated error checking
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    //Find omega vector, i.e. omega = w + (inverse of (trace(V) I - V))
    omegaX =  wX + tempInvXX*zX + tempInvXY*zY + tempInvXZ*zZ;
    omegaY =  wY + tempInvYX*zX + tempInvYY*zY + tempInvYZ*zZ;
    omegaZ =  wZ + tempInvZX*zX + tempInvZY*zY + tempInvZZ*zZ;

    //Find the omega tensor from it's vector dual
    omegaTensorXX = 0.0;
    omegaTensorXY = -omegaZ;
    omegaTensorXZ = omegaY;
    omegaTensorYX = omegaZ;
    omegaTensorYY = 0.0;
    omegaTensorYZ = -omegaX;
    omegaTensorZX = -omegaY;
    omegaTensorZY = omegaX;
    omegaTensorXX = 0.0;

    //Increment R
    //
    // Compute tempA = I + 0.5 * dt * omega
    scaleFactor = 0.5*dt;
    MatrixUpdate(ScalarT(1.0), scaleFactor,
                 identityMatrixXX, identityMatrixXY, identityMatrixXZ,
                 identityMatrixYX, identityMatrixYY, identityMatrixYZ,
                 identityMatrixZX, identityMatrixZY, identityMatrixZZ,
                 omegaTensorXX, omegaTensorYX, omegaTensorZX,
                 omegaTensorYX, omegaTensorYY, omegaTensorZY,
                 omegaTensorXZ, omegaTensorYZ, omegaTensorZZ,
                 tempAXX, tempAXY, tempAXZ,
                 tempAYX, tempAYY, tempAYZ,
                 tempAZX, tempAZY, tempAZZ);

    // Compute tempB = I - 0.5 * dt * omega
    scaleFactor = -0.5*dt;
    MatrixUpdate(ScalarT(1.0), scaleFactor,
                 identityMatrixXX, identityMatrixXY, identityMatrixXZ,
                 identityMatrixYX, identityMatrixYY, identityMatrixYZ,
                 identityMatrixZX, identityMatrixZY, identityMatrixZZ,
                 omegaTensorXX, omegaTensorYX, omegaTensorZX,
                 omegaTensorYX, omegaTensorYY, omegaTensorZY,
                 omegaTensorXZ, omegaTensorYZ, omegaTensorZZ,
                 tempBXX, tempBXY, tempBXZ,
                 tempBYX, tempBYY, tempBYZ,
                 tempBZX, tempBZY, tempBZZ);

    // Compute the inverse of the tempB matrix
    inversionReturnCode = invert3by3Matrix(tempBXX, tempBXY, tempBXZ,
                                           tempBYX, tempBYY, tempBYZ,
                                           tempBZX, tempBZY, tempBZZ,
                                           tempBInvXX, tempBInvXY, tempBInvXZ,
                                           tempBInvYX, tempBInvYY, tempBInvYZ,
                                           tempBInvZX, tempBInvZY, tempBInvZZ);
    // Placeholder for more sophisticated error checking
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    // Compute temp = tempBinv * tempA
    MatrixMultiply(tempBInvXX, tempBInvXY, tempBInvXZ,
                   tempBInvYX, tempBInvYY, tempBInvYZ,
                   tempBInvZX, tempBInvZY, tempBInvZZ,
                   tempAXX, tempAXY, tempAXZ,
                   tempAYX, tempAYY, tempAYZ,
                   tempAZX, tempAZY, tempAZZ,
                   tempXX, tempXY, tempXZ,
                   tempYX, tempYY, tempYZ,
                   tempZX, tempZY, tempZZ);

    // Compute R = temp * Rold
    MatrixMultiply(tempInvXX, tempInvXY, tempInvXZ,
                   tempInvYX, tempInvYY, tempInvYZ,
                   tempInvZX, tempInvZY, tempInvZZ,
                   rotTensorOldXX, rotTensorOldXY, rotTensorOldXZ,
                   rotTensorOldYX, rotTensorOldYY, rotTensorOldYZ,
                   rotTensorOldZX, rotTensorOldZY, rotTensorOldZZ,
                   *rotTensorXX, *rotTensorXY, *rotTensorXZ,
                   *rotTensorYX, *rotTensorYY, *rotTensorYZ,
                   *rotTensorZX, *rotTensorZY, *rotTensorZZ);

    // Compute the unrotated rate-of-deformation, d, i.e., temp = Rt * D * R
    UnrotateTensor(rateOfDefXX, rateOfDefXY, rateOfDefXZ,
                 rateOfDefYX, rateOfDefYY, rateOfDefYZ,
                 rateOfDefZX, rateOfDefZY, rateOfDefZZ,
                 *rotTensorXX, *rotTensorXY, *rotTensorXZ,
                 *rotTensorYX, *rotTensorYY, *rotTensorYZ,
                 *rotTensorZX, *rotTensorZY, *rotTensorZZ,
                 *unrotRateOfDefXX, *unrotRateOfDefXY, *unrotRateOfDefXZ,
                 *unrotRateOfDefYX, *unrotRateOfDefYY, *unrotRateOfDefYZ,
                 *unrotRateOfDefZX, *unrotRateOfDefZY, *unrotRateOfDefZZ);
    

    // Find V = F * Rt
    MatrixMultiply(*defGradXX, *defGradYX, *defGradZX,
                   *defGradXY, *defGradYY, *defGradZY,
                   *defGradXZ, *defGradYZ, *defGradZZ,
                   *rotTensorXX, *rotTensorYX, *rotTensorZX,
                   *rotTensorXY, *rotTensorYY, *rotTensorZY,
                   *rotTensorXZ, *rotTensorYZ, *rotTensorZZ,
                   *leftStretchXX, *leftStretchXY, *leftStretchXZ,
                   *leftStretchYX, *leftStretchYY, *leftStretchYZ,
                   *leftStretchZX, *leftStretchZY, *leftStretchZZ);
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
 ScalarT* cauchyStressXX,
 ScalarT* cauchyStressXY,
 ScalarT* cauchyStressXZ,
 ScalarT* cauchyStressYX,
 ScalarT* cauchyStressYY,
 ScalarT* cauchyStressYZ,
 ScalarT* cauchyStressZX,
 ScalarT* cauchyStressZY,
 ScalarT* cauchyStressZZ,
 int numPoints
)
{
  const ScalarT* rotTensorXX = rotationTensorXX;
  const ScalarT* rotTensorXY = rotationTensorXY;
  const ScalarT* rotTensorXZ = rotationTensorXZ;
  const ScalarT* rotTensorYX = rotationTensorYX;
  const ScalarT* rotTensorYY = rotationTensorYY;
  const ScalarT* rotTensorYZ = rotationTensorYZ;
  const ScalarT* rotTensorZX = rotationTensorZX;
  const ScalarT* rotTensorZY = rotationTensorZY;
  const ScalarT* rotTensorZZ = rotationTensorZZ;
  ScalarT* rotatedStressXX = cauchyStressXX;
  ScalarT* rotatedStressXY = cauchyStressXY;
  ScalarT* rotatedStressXZ = cauchyStressXZ;
  ScalarT* rotatedStressYX = cauchyStressYX;
  ScalarT* rotatedStressYY = cauchyStressYY;
  ScalarT* rotatedStressYZ = cauchyStressYZ;
  ScalarT* rotatedStressZX = cauchyStressZX;
  ScalarT* rotatedStressZY = cauchyStressZY;
  ScalarT* rotatedStressZZ = cauchyStressZZ;

  ScalarT unrotatedStressXX;
  ScalarT unrotatedStressXY;
  ScalarT unrotatedStressXZ;
  ScalarT unrotatedStressYX;
  ScalarT unrotatedStressYY;
  ScalarT unrotatedStressYZ;
  ScalarT unrotatedStressZX;
  ScalarT unrotatedStressZY;
  ScalarT unrotatedStressZZ;

  for(int iID=0 ; iID<numPoints ; ++iID, 
        ++rotTensorXX, ++rotTensorXY, ++rotTensorXZ,
        ++rotTensorYX, ++rotTensorYY, ++rotTensorYZ,
        ++rotTensorZX, ++rotTensorZY, ++rotTensorZZ,
        ++rotatedStressXX, ++rotatedStressXY, ++rotatedStressXZ,
        ++rotatedStressYX, ++rotatedStressYY, ++rotatedStressYZ,
        ++rotatedStressZX, ++rotatedStressZY, ++rotatedStressZZ){ 
      
      unrotatedStressXX = *rotatedStressXX;
      unrotatedStressXY = *rotatedStressXY;
      unrotatedStressXZ = *rotatedStressXZ;
      unrotatedStressYX = *rotatedStressYX;
      unrotatedStressYY = *rotatedStressYY;
      unrotatedStressYZ = *rotatedStressYZ;
      unrotatedStressZX = *rotatedStressZX;
      unrotatedStressZY = *rotatedStressZY;
      unrotatedStressZZ = *rotatedStressZZ;
      
      RotateTensor(unrotatedStressXX, unrotatedStressXY, unrotatedStressXZ, 
                   unrotatedStressYX, unrotatedStressYY, unrotatedStressYZ, 
                   unrotatedStressZX, unrotatedStressZY, unrotatedStressZZ,
                   *rotTensorXX, *rotTensorXY, *rotTensorXZ,
                   *rotTensorYX, *rotTensorYY, *rotTensorYZ,
                   *rotTensorZX, *rotTensorZY, *rotTensorZZ,
                   *rotatedStressXX, *rotatedStressXY, *rotatedStressXZ,
                   *rotatedStressYX, *rotatedStressYY, *rotatedStressYZ,
                   *rotatedStressZX, *rotatedStressZY, *rotatedStressZZ);
  }
}

template<typename ScalarT>
void unrotateCauchyStress
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
 ScalarT* cauchyStressXX,
 ScalarT* cauchyStressXY,
 ScalarT* cauchyStressXZ,
 ScalarT* cauchyStressYX,
 ScalarT* cauchyStressYY,
 ScalarT* cauchyStressYZ,
 ScalarT* cauchyStressZX,
 ScalarT* cauchyStressZY,
 ScalarT* cauchyStressZZ,
 int numPoints
)
{
  const ScalarT* rotTensorXX = rotationTensorXX;
  const ScalarT* rotTensorXY = rotationTensorXY;
  const ScalarT* rotTensorXZ = rotationTensorXZ;
  const ScalarT* rotTensorYX = rotationTensorYX;
  const ScalarT* rotTensorYY = rotationTensorYY;
  const ScalarT* rotTensorYZ = rotationTensorYZ;
  const ScalarT* rotTensorZX = rotationTensorZX;
  const ScalarT* rotTensorZY = rotationTensorZY;
  const ScalarT* rotTensorZZ = rotationTensorZZ;
  ScalarT* unrotatedStressXX = cauchyStressXX;
  ScalarT* unrotatedStressXY = cauchyStressXY;
  ScalarT* unrotatedStressXZ = cauchyStressXZ;
  ScalarT* unrotatedStressYX = cauchyStressYX;
  ScalarT* unrotatedStressYY = cauchyStressYY;
  ScalarT* unrotatedStressYZ = cauchyStressYZ;
  ScalarT* unrotatedStressZX = cauchyStressZX;
  ScalarT* unrotatedStressZY = cauchyStressZY;
  ScalarT* unrotatedStressZZ = cauchyStressZZ;

  ScalarT rotatedStressXX;
  ScalarT rotatedStressXY;
  ScalarT rotatedStressXZ;
  ScalarT rotatedStressYX;
  ScalarT rotatedStressYY;
  ScalarT rotatedStressYZ;
  ScalarT rotatedStressZX;
  ScalarT rotatedStressZY;
  ScalarT rotatedStressZZ;

  for(int iID=0 ; iID<numPoints ; ++iID, 
        ++rotTensorXX, ++rotTensorXY, ++rotTensorXZ,
        ++rotTensorYX, ++rotTensorYY, ++rotTensorYZ,
        ++rotTensorZX, ++rotTensorZY, ++rotTensorZZ,
        ++unrotatedStressXX, ++unrotatedStressXY, ++unrotatedStressXZ,
        ++unrotatedStressYX, ++unrotatedStressYY, ++unrotatedStressYZ,
        ++unrotatedStressZX, ++unrotatedStressZY, ++unrotatedStressZZ){ 
      
      rotatedStressXX = *unrotatedStressXX;
      rotatedStressXY = *unrotatedStressXY;
      rotatedStressXZ = *unrotatedStressXZ;
      rotatedStressYX = *unrotatedStressYX;
      rotatedStressYY = *unrotatedStressYY;
      rotatedStressYZ = *unrotatedStressYZ;
      rotatedStressZX = *unrotatedStressZX;
      rotatedStressZY = *unrotatedStressZY;
      rotatedStressZZ = *unrotatedStressZZ;
      
      UnrotateTensor(rotatedStressXX, rotatedStressXY, rotatedStressXZ, 
                     rotatedStressYX, rotatedStressYY, rotatedStressYZ, 
                     rotatedStressZX, rotatedStressZY, rotatedStressZZ,
                     *rotTensorXX, *rotTensorXY, *rotTensorXZ,
                     *rotTensorYX, *rotTensorYY, *rotTensorYZ,
                     *rotTensorZX, *rotTensorZY, *rotTensorZZ,
                     *unrotatedStressXX, *unrotatedStressXY, *unrotatedStressXZ,
                     *unrotatedStressYX, *unrotatedStressYY, *unrotatedStressYZ,
                     *unrotatedStressZX, *unrotatedStressZY, *unrotatedStressZZ);
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

template void MatrixUpdate<double>
(
 const double& alpha,
 const double& beta,
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

template void UnrotateTensor<double>
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
 const double& rXX,
 const double& rXY,
 const double& rXZ,
 const double& rYX,
 const double& rYY,
 const double& rYZ,
 const double& rZX,
 const double& rZY,
 const double& rZZ,
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

template void RotateTensor<double>
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
 const double& rXX,
 const double& rXY,
 const double& rXZ,
 const double& rYX,
 const double& rYY,
 const double& rYZ,
 const double& rZX,
 const double& rZY,
 const double& rZZ,
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

template void unrotateCauchyStress<double>
(
 const double* rotationTensorXX,
 const double* rotationTensorXY,
 const double* rotationTensorXZ,
 const double* rotationTensorYX,
 const double* rotationTensorYY,
 const double* rotationTensorYZ,
 const double* rotationTensorZX,
 const double* rotationTensorZY,
 const double* rotationTensorZZ,
 double* cauchyStressXX,
 double* cauchyStressXY,
 double* cauchyStressXZ,
 double* cauchyStressYX,
 double* cauchyStressYY,
 double* cauchyStressYZ,
 double* cauchyStressZX,
 double* cauchyStressZY,
 double* cauchyStressZZ,
 int numPoints
);

template void rotateCauchyStress<double>
(
 const double* rotationTensorXX,
 const double* rotationTensorXY,
 const double* rotationTensorXZ,
 const double* rotationTensorYX,
 const double* rotationTensorYY,
 const double* rotationTensorYZ,
 const double* rotationTensorZX,
 const double* rotationTensorZY,
 const double* rotationTensorZZ,
 double* cauchyStressXX,
 double* cauchyStressXY,
 double* cauchyStressXZ,
 double* cauchyStressYX,
 double* cauchyStressYY,
 double* cauchyStressYZ,
 double* cauchyStressZX,
 double* cauchyStressZY,
 double* cauchyStressZZ,
 int numPoints
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

template int computeUnrotatedRateOfDeformationAndRotationTensor<double>
(
const double* volume,
const double* modelCoordinates,
const double* coordinates,
const double* velocities,
const double* deformationGradientXX,
const double* deformationGradientXY,
const double* deformationGradientXZ,
const double* deformationGradientYX,
const double* deformationGradientYY,
const double* deformationGradientYZ,
const double* deformationGradientZX,
const double* deformationGradientZY,
const double* deformationGradientZZ,
const double* shapeTensorInverseXX,
const double* shapeTensorInverseXY,
const double* shapeTensorInverseXZ,
const double* shapeTensorInverseYX,
const double* shapeTensorInverseYY,
const double* shapeTensorInverseYZ,
const double* shapeTensorInverseZX,
const double* shapeTensorInverseZY,
const double* shapeTensorInverseZZ,
double* leftStretchTensorXX,
double* leftStretchTensorXY,
double* leftStretchTensorXZ,
double* leftStretchTensorYX,
double* leftStretchTensorYY,
double* leftStretchTensorYZ,
double* leftStretchTensorZX,
double* leftStretchTensorZY,
double* leftStretchTensorZZ,
double* rotationTensorXX,
double* rotationTensorXY,
double* rotationTensorXZ,
double* rotationTensorYX,
double* rotationTensorYY,
double* rotationTensorYZ,
double* rotationTensorZX,
double* rotationTensorZY,
double* rotationTensorZZ,
double* unrotatedRateOfDeformationXX,
double* unrotatedRateOfDeformationXY,
double* unrotatedRateOfDeformationXZ,
double* unrotatedRateOfDeformationYX,
double* unrotatedRateOfDeformationYY,
double* unrotatedRateOfDeformationYZ,
double* unrotatedRateOfDeformationZX,
double* unrotatedRateOfDeformationZY,
double* unrotatedRateOfDeformationZZ,
const int* neighborhoodList,
int numPoints,
double horizon,
double dt
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

template void UnrotateTensor<Sacado::Fad::DFad<double> >
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
 const Sacado::Fad::DFad<double>& rXX,
 const Sacado::Fad::DFad<double>& rXY,
 const Sacado::Fad::DFad<double>& rXZ,
 const Sacado::Fad::DFad<double>& rYX,
 const Sacado::Fad::DFad<double>& rYY,
 const Sacado::Fad::DFad<double>& rYZ,
 const Sacado::Fad::DFad<double>& rZX,
 const Sacado::Fad::DFad<double>& rZY,
 const Sacado::Fad::DFad<double>& rZZ,
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

template void unrotateCauchyStress<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* rotationTensorXX,
 const Sacado::Fad::DFad<double>* rotationTensorXY,
 const Sacado::Fad::DFad<double>* rotationTensorXZ,
 const Sacado::Fad::DFad<double>* rotationTensorYX,
 const Sacado::Fad::DFad<double>* rotationTensorYY,
 const Sacado::Fad::DFad<double>* rotationTensorYZ,
 const Sacado::Fad::DFad<double>* rotationTensorZX,
 const Sacado::Fad::DFad<double>* rotationTensorZY,
 const Sacado::Fad::DFad<double>* rotationTensorZZ,
 Sacado::Fad::DFad<double>* cauchyStressXX,
 Sacado::Fad::DFad<double>* cauchyStressXY,
 Sacado::Fad::DFad<double>* cauchyStressXZ,
 Sacado::Fad::DFad<double>* cauchyStressYX,
 Sacado::Fad::DFad<double>* cauchyStressYY,
 Sacado::Fad::DFad<double>* cauchyStressYZ,
 Sacado::Fad::DFad<double>* cauchyStressZX,
 Sacado::Fad::DFad<double>* cauchyStressZY,
 Sacado::Fad::DFad<double>* cauchyStressZZ,
 int numPoints
);

template void rotateCauchyStress<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>* rotationTensorXX,
 const Sacado::Fad::DFad<double>* rotationTensorXY,
 const Sacado::Fad::DFad<double>* rotationTensorXZ,
 const Sacado::Fad::DFad<double>* rotationTensorYX,
 const Sacado::Fad::DFad<double>* rotationTensorYY,
 const Sacado::Fad::DFad<double>* rotationTensorYZ,
 const Sacado::Fad::DFad<double>* rotationTensorZX,
 const Sacado::Fad::DFad<double>* rotationTensorZY,
 const Sacado::Fad::DFad<double>* rotationTensorZZ,
 Sacado::Fad::DFad<double>* cauchyStressXX,
 Sacado::Fad::DFad<double>* cauchyStressXY,
 Sacado::Fad::DFad<double>* cauchyStressXZ,
 Sacado::Fad::DFad<double>* cauchyStressYX,
 Sacado::Fad::DFad<double>* cauchyStressYY,
 Sacado::Fad::DFad<double>* cauchyStressYZ,
 Sacado::Fad::DFad<double>* cauchyStressZX,
 Sacado::Fad::DFad<double>* cauchyStressZY,
 Sacado::Fad::DFad<double>* cauchyStressZZ,
 int numPoints
);

template void RotateTensor<Sacado::Fad::DFad<double> >
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
 const Sacado::Fad::DFad<double>& rXX,
 const Sacado::Fad::DFad<double>& rXY,
 const Sacado::Fad::DFad<double>& rXZ,
 const Sacado::Fad::DFad<double>& rYX,
 const Sacado::Fad::DFad<double>& rYY,
 const Sacado::Fad::DFad<double>& rYZ,
 const Sacado::Fad::DFad<double>& rZX,
 const Sacado::Fad::DFad<double>& rZY,
 const Sacado::Fad::DFad<double>& rZZ,
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

template void MatrixUpdate<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double>& alpha,
 const Sacado::Fad::DFad<double>& beta,
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


template int computeUnrotatedRateOfDeformationAndRotationTensor<Sacado::Fad::DFad<double> >
(
const double* volume,
const double* modelCoordinates,
const Sacado::Fad::DFad<double>* coordinates,
const Sacado::Fad::DFad<double>* velocities,
const Sacado::Fad::DFad<double>* deformationGradientXX,
const Sacado::Fad::DFad<double>* deformationGradientXY,
const Sacado::Fad::DFad<double>* deformationGradientXZ,
const Sacado::Fad::DFad<double>* deformationGradientYX,
const Sacado::Fad::DFad<double>* deformationGradientYY,
const Sacado::Fad::DFad<double>* deformationGradientYZ,
const Sacado::Fad::DFad<double>* deformationGradientZX,
const Sacado::Fad::DFad<double>* deformationGradientZY,
const Sacado::Fad::DFad<double>* deformationGradientZZ,
const Sacado::Fad::DFad<double>* shapeTensorInverseXX,
const Sacado::Fad::DFad<double>* shapeTensorInverseXY,
const Sacado::Fad::DFad<double>* shapeTensorInverseXZ,
const Sacado::Fad::DFad<double>* shapeTensorInverseYX,
const Sacado::Fad::DFad<double>* shapeTensorInverseYY,
const Sacado::Fad::DFad<double>* shapeTensorInverseYZ,
const Sacado::Fad::DFad<double>* shapeTensorInverseZX,
const Sacado::Fad::DFad<double>* shapeTensorInverseZY,
const Sacado::Fad::DFad<double>* shapeTensorInverseZZ,
Sacado::Fad::DFad<double>* leftStretchTensorXX,
Sacado::Fad::DFad<double>* leftStretchTensorXY,
Sacado::Fad::DFad<double>* leftStretchTensorXZ,
Sacado::Fad::DFad<double>* leftStretchTensorYX,
Sacado::Fad::DFad<double>* leftStretchTensorYY,
Sacado::Fad::DFad<double>* leftStretchTensorYZ,
Sacado::Fad::DFad<double>* leftStretchTensorZX,
Sacado::Fad::DFad<double>* leftStretchTensorZY,
Sacado::Fad::DFad<double>* leftStretchTensorZZ,
Sacado::Fad::DFad<double>* rotationTensorXX,
Sacado::Fad::DFad<double>* rotationTensorXY,
Sacado::Fad::DFad<double>* rotationTensorXZ,
Sacado::Fad::DFad<double>* rotationTensorYX,
Sacado::Fad::DFad<double>* rotationTensorYY,
Sacado::Fad::DFad<double>* rotationTensorYZ,
Sacado::Fad::DFad<double>* rotationTensorZX,
Sacado::Fad::DFad<double>* rotationTensorZY,
Sacado::Fad::DFad<double>* rotationTensorZZ,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXX,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXY,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXZ,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYX,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYY,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYZ,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZX,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZY,
Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZZ,
const int* neighborhoodList,
int numPoints,
double horizon,
double dt
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
