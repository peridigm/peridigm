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
#include <Teuchos_ScalarTraits.hpp>
#include <math.h>
#include <functional>
#include <vector>

namespace CORRESPONDENCE {

template<typename ScalarT>
void setOnesOnDiagonalFullTensor(ScalarT* tensor, int numPoints){

  ScalarT *tens = tensor;

  for(int iID=0; iID<numPoints; ++iID, tens+=9){
    *(tens) = 1.0;
    *(tens+4) = 1.0;
    *(tens+8) = 1.0;
  };
}

template<typename ScalarT>
int Invert3by3Matrix
(
    const ScalarT* matrix,
    ScalarT& determinant,
    ScalarT* inverse
)
{
  int returnCode(0);

  ScalarT minor0 =  *(matrix+4) * *(matrix+8) - *(matrix+5) * *(matrix+7);
  ScalarT minor1 =  *(matrix+3) * *(matrix+8) - *(matrix+5) * *(matrix+6);
  ScalarT minor2 =  *(matrix+3) * *(matrix+7) - *(matrix+4) * *(matrix+6);
  ScalarT minor3 =  *(matrix+1) * *(matrix+8) - *(matrix+2) * *(matrix+7);
  ScalarT minor4 =  *(matrix)   * *(matrix+8) - *(matrix+6) * *(matrix+2);
  ScalarT minor5 =  *(matrix)   * *(matrix+7) - *(matrix+1) * *(matrix+6);
  ScalarT minor6 =  *(matrix+1) * *(matrix+5) - *(matrix+2) * *(matrix+4);
  ScalarT minor7 =  *(matrix)   * *(matrix+5) - *(matrix+2) * *(matrix+3);
  ScalarT minor8 =  *(matrix)   * *(matrix+4) - *(matrix+1) * *(matrix+3);
  determinant = *(matrix) * minor0 - *(matrix+1) * minor1 + *(matrix+2) * minor2;

  if(determinant == ScalarT(0.0)){
    returnCode = 1;
    *(inverse) = 0.0;
    *(inverse+1) = 0.0;
    *(inverse+2) = 0.0;
    *(inverse+3) = 0.0;
    *(inverse+4) = 0.0;
    *(inverse+5) = 0.0;
    *(inverse+6) = 0.0;
    *(inverse+7) = 0.0;
    *(inverse+8) = 0.0;
  }
  else{
    *(inverse) = minor0/determinant;
    *(inverse+1) = -1.0*minor3/determinant;
    *(inverse+2) = minor6/determinant;
    *(inverse+3) = -1.0*minor1/determinant;
    *(inverse+4) = minor4/determinant;
    *(inverse+5) = -1.0*minor7/determinant;
    *(inverse+6) = minor2/determinant;
    *(inverse+7) = -1.0*minor5/determinant;
    *(inverse+8) = minor8/determinant;
  }

  return returnCode;
}

template<typename ScalarT>
void TransposeMatrix
(
 const ScalarT* matrix,
 ScalarT* transpose
)
{
  // Store some values so that the matrix and transpose can be the
  // same matrix (i.e., transpose in place)
  ScalarT temp_xy( *(matrix+1) );
  ScalarT temp_xz( *(matrix+2) );
  ScalarT temp_yz( *(matrix+5) );

  *(transpose)   = *(matrix);
  *(transpose+1) = *(matrix+3);
  *(transpose+2) = *(matrix+6);
  *(transpose+3) = temp_xy;
  *(transpose+4) = *(matrix+4);
  *(transpose+5) = *(matrix+7);
  *(transpose+6) = temp_xz;
  *(transpose+7) = temp_yz;
  *(transpose+8) = *(matrix+8);
}

template<typename ScalarT>
void MatrixMultiply
(
    bool transA,
    bool transB,
    ScalarT alpha,
    const ScalarT* a,
    const ScalarT* b,
    ScalarT* result
)
{
  // This function computes result = alpha * a * b
  // where alpha is a scalar and a and b are 3x3 matrices
  // The arguments transA and transB denote whether or not
  // to use the transpose of a and b, respectively.

  // The default ordering is row-major:
  //
  // XX(0) XY(1) XZ(2)
  // YX(3) YY(4) YZ(5)
  // ZX(6) ZY(7) ZZ(8)

  if(!transA && !transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+1) * *(b+3) + *(a+2) * *(b+6);
    *(result+1) = *(a+0) * *(b+1) + *(a+1) * *(b+4) + *(a+2) * *(b+7);
    *(result+2) = *(a+0) * *(b+2) + *(a+1) * *(b+5) + *(a+2) * *(b+8);
    *(result+3) = *(a+3) * *(b+0) + *(a+4) * *(b+3) + *(a+5) * *(b+6);
    *(result+4) = *(a+3) * *(b+1) + *(a+4) * *(b+4) + *(a+5) * *(b+7);
    *(result+5) = *(a+3) * *(b+2) + *(a+4) * *(b+5) + *(a+5) * *(b+8);
    *(result+6) = *(a+6) * *(b+0) + *(a+7) * *(b+3) + *(a+8) * *(b+6);
    *(result+7) = *(a+6) * *(b+1) + *(a+7) * *(b+4) + *(a+8) * *(b+7);
    *(result+8) = *(a+6) * *(b+2) + *(a+7) * *(b+5) + *(a+8) * *(b+8);
  }
  else if(transA && !transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+3) * *(b+3) + *(a+6) * *(b+6);
    *(result+1) = *(a+0) * *(b+1) + *(a+3) * *(b+4) + *(a+6) * *(b+7);
    *(result+2) = *(a+0) * *(b+2) + *(a+3) * *(b+5) + *(a+6) * *(b+8);
    *(result+3) = *(a+1) * *(b+0) + *(a+4) * *(b+3) + *(a+7) * *(b+6);
    *(result+4) = *(a+1) * *(b+1) + *(a+4) * *(b+4) + *(a+7) * *(b+7);
    *(result+5) = *(a+1) * *(b+2) + *(a+4) * *(b+5) + *(a+7) * *(b+8);
    *(result+6) = *(a+2) * *(b+0) + *(a+5) * *(b+3) + *(a+8) * *(b+6);
    *(result+7) = *(a+2) * *(b+1) + *(a+5) * *(b+4) + *(a+8) * *(b+7);
    *(result+8) = *(a+2) * *(b+2) + *(a+5) * *(b+5) + *(a+8) * *(b+8);
  }
  else if(!transA && transB){
    *(result+0) = *(a+0) * *(b+0) + *(a+1) * *(b+1) + *(a+2) * *(b+2);
    *(result+1) = *(a+0) * *(b+3) + *(a+1) * *(b+4) + *(a+2) * *(b+5);
    *(result+2) = *(a+0) * *(b+6) + *(a+1) * *(b+7) + *(a+2) * *(b+8);
    *(result+3) = *(a+3) * *(b+0) + *(a+4) * *(b+1) + *(a+5) * *(b+2);
    *(result+4) = *(a+3) * *(b+3) + *(a+4) * *(b+4) + *(a+5) * *(b+5);
    *(result+5) = *(a+3) * *(b+6) + *(a+4) * *(b+7) + *(a+5) * *(b+8);
    *(result+6) = *(a+6) * *(b+0) + *(a+7) * *(b+1) + *(a+8) * *(b+2);
    *(result+7) = *(a+6) * *(b+3) + *(a+7) * *(b+4) + *(a+8) * *(b+5);
    *(result+8) = *(a+6) * *(b+6) + *(a+7) * *(b+7) + *(a+8) * *(b+8);
  }
  else{
    *(result+0) = *(a+0) * *(b+0) + *(a+3) * *(b+1) + *(a+6) * *(b+2);
    *(result+1) = *(a+0) * *(b+3) + *(a+3) * *(b+4) + *(a+6) * *(b+5);
    *(result+2) = *(a+0) * *(b+6) + *(a+3) * *(b+7) + *(a+6) * *(b+8);
    *(result+3) = *(a+1) * *(b+0) + *(a+4) * *(b+1) + *(a+7) * *(b+2);
    *(result+4) = *(a+1) * *(b+3) + *(a+4) * *(b+4) + *(a+7) * *(b+5);
    *(result+5) = *(a+1) * *(b+6) + *(a+4) * *(b+7) + *(a+7) * *(b+8);
    *(result+6) = *(a+2) * *(b+0) + *(a+5) * *(b+1) + *(a+8) * *(b+2);
    *(result+7) = *(a+2) * *(b+3) + *(a+5) * *(b+4) + *(a+8) * *(b+5);
    *(result+8) = *(a+2) * *(b+6) + *(a+5) * *(b+7) + *(a+8) * *(b+8);
  }

  if(alpha != 1.0){
    for(int i=0 ; i<9 ; ++i)
      *(result+i) *= alpha;
  }
}

template<typename ScalarT>
int computeShapeTensorInverseAndApproximateDeformationGradient
(
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const ScalarT* coordinates,
    ScalarT* shapeTensorInverse,
    ScalarT* deformationGradient,
    const int* neighborhoodList,
    int numPoints
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* defGrad = deformationGradient;

  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  ScalarT deformedBondX, deformedBondY, deformedBondZ;
  double neighborVolume, omega, temp;

  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> defGradFirstTermVector(9);
  ScalarT* defGradFirstTerm = &defGradFirstTermVector[0];

  // placeholder for bond damage
  double bondDamage = 0.0;

  int inversionReturnCode(0);

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
      shapeTensorInv+=9, defGrad+=9){

    // Zero out data
    *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
    *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
    *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
    *(defGradFirstTerm)   = 0.0 ; *(defGradFirstTerm+1) = 0.0 ; *(defGradFirstTerm+2) = 0.0 ;
    *(defGradFirstTerm+3) = 0.0 ; *(defGradFirstTerm+4) = 0.0 ; *(defGradFirstTerm+5) = 0.0 ;
    *(defGradFirstTerm+6) = 0.0 ; *(defGradFirstTerm+7) = 0.0 ; *(defGradFirstTerm+8) = 0.0 ;

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

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);

      temp = (1.0 - bondDamage) * omega * neighborVolume;

      *(shapeTensor)   += temp * undeformedBondX * undeformedBondX;
      *(shapeTensor+1) += temp * undeformedBondX * undeformedBondY;
      *(shapeTensor+2) += temp * undeformedBondX * undeformedBondZ;
      *(shapeTensor+3) += temp * undeformedBondY * undeformedBondX;
      *(shapeTensor+4) += temp * undeformedBondY * undeformedBondY;
      *(shapeTensor+5) += temp * undeformedBondY * undeformedBondZ;
      *(shapeTensor+6) += temp * undeformedBondZ * undeformedBondX;
      *(shapeTensor+7) += temp * undeformedBondZ * undeformedBondY;
      *(shapeTensor+8) += temp * undeformedBondZ * undeformedBondZ;

      *(defGradFirstTerm)   += temp * deformedBondX * undeformedBondX;
      *(defGradFirstTerm+1) += temp * deformedBondX * undeformedBondY;
      *(defGradFirstTerm+2) += temp * deformedBondX * undeformedBondZ;
      *(defGradFirstTerm+3) += temp * deformedBondY * undeformedBondX;
      *(defGradFirstTerm+4) += temp * deformedBondY * undeformedBondY;
      *(defGradFirstTerm+5) += temp * deformedBondY * undeformedBondZ;
      *(defGradFirstTerm+6) += temp * deformedBondZ * undeformedBondX;
      *(defGradFirstTerm+7) += temp * deformedBondZ * undeformedBondY;
      *(defGradFirstTerm+8) += temp * deformedBondZ * undeformedBondZ;
    }
    
    inversionReturnCode = Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    // Matrix multiply the first term and the shape tensor inverse to compute
    // the deformation gradient
    MatrixMultiply(false, false, 1.0, defGradFirstTerm, shapeTensorInv, defGrad);
  }

  return returnCode;
}

//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
template<typename ScalarT>
int computeUnrotatedRateOfDeformationAndRotationTensor(
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const ScalarT* velocities,
    const ScalarT* deformationGradient,
    const ScalarT* shapeTensorInverse,
    const ScalarT* leftStretchTensorN,
    const ScalarT* rotationTensorN,
    ScalarT* leftStretchTensorNP1,
    ScalarT* rotationTensorNP1,
    ScalarT* unrotatedRateOfDeformation,
    const int* neighborhoodList,
    int numPoints,
    double dt
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* defGrad = deformationGradient;
  const ScalarT* shapeTensorInv = shapeTensorInverse;
  const ScalarT* leftStretchN = leftStretchTensorN;
  const ScalarT* rotTensorN = rotationTensorN;

  ScalarT* leftStretchNP1 = leftStretchTensorNP1;
  ScalarT* rotTensorNP1 = rotationTensorNP1;
  ScalarT* unrotRateOfDef = unrotatedRateOfDeformation;

  std::vector<ScalarT> FdotFirstTermVector(9) ; ScalarT* FdotFirstTerm = &FdotFirstTermVector[0];
  std::vector<ScalarT> FdotVector(9) ; ScalarT* Fdot = &FdotVector[0];
  std::vector<ScalarT> FinverseVector(9) ; ScalarT* Finverse = &FinverseVector[0];
  std::vector<ScalarT> eulerianVelGradVector(9) ; ScalarT* eulerianVelGrad = &eulerianVelGradVector[0];
  std::vector<ScalarT> rateOfDefVector(9) ; ScalarT* rateOfDef = &rateOfDefVector[0];
  std::vector<ScalarT> spinVector(9) ; ScalarT* spin = &spinVector[0];
  std::vector<ScalarT> tempVector(9) ; ScalarT* temp = &tempVector[0];
  std::vector<ScalarT> tempInvVector(9) ; ScalarT* tempInv = &tempInvVector[0];
  std::vector<ScalarT> OmegaTensorVector(9) ; ScalarT* OmegaTensor = &OmegaTensorVector[0];
  std::vector<ScalarT> QMatrixVector(9) ; ScalarT* QMatrix = &QMatrixVector[0];
  std::vector<ScalarT> OmegaTensorSqVector(9) ; ScalarT* OmegaTensorSq = &OmegaTensorSqVector[0];
  std::vector<ScalarT> tempAVector(9) ; ScalarT* tempA = &tempAVector[0];
  std::vector<ScalarT> tempBVector(9) ; ScalarT* tempB = &tempBVector[0];
  std::vector<ScalarT> rateOfStretchVector(9) ; ScalarT* rateOfStretch = &rateOfStretchVector[0];

  ScalarT determinant;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  double undeformedBondX, undeformedBondY, undeformedBondZ, undeformedBondLength;
  double neighborVolume, omega, scalarTemp; 
  int inversionReturnCode(0);

  // placeholder for bond damage
  double bondDamage = 0.0;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, vel+=3,
      shapeTensorInv+=9, rotTensorN+=9, rotTensorNP1+=9, leftStretchNP1+=9, leftStretchN+=9,
      unrotRateOfDef+=9, defGrad+=9){

    // Initialize data
    *(FdotFirstTerm)   = 0.0 ; *(FdotFirstTerm+1) = 0.0 ;  *(FdotFirstTerm+2) = 0.0;
    *(FdotFirstTerm+3) = 0.0 ; *(FdotFirstTerm+4) = 0.0 ;  *(FdotFirstTerm+5) = 0.0;
    *(FdotFirstTerm+6) = 0.0 ; *(FdotFirstTerm+7) = 0.0 ;  *(FdotFirstTerm+8) = 0.0;
    
    //Compute Fdot
    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = volume[neighborIndex];
      neighborModelCoord = modelCoordinates + 3*neighborIndex;
      neighborVel = velocities + 3*neighborIndex;

      undeformedBondX = *(neighborModelCoord)   - *(modelCoord);
      undeformedBondY = *(neighborModelCoord+1) - *(modelCoord+1);
      undeformedBondZ = *(neighborModelCoord+2) - *(modelCoord+2);
      undeformedBondLength = sqrt(undeformedBondX*undeformedBondX +
                                  undeformedBondY*undeformedBondY +
                                  undeformedBondZ*undeformedBondZ);

      // The velState is the relative difference in velocities of the nodes at
      // each end of a bond. i.e., v_j - v_i
      velStateX = *(neighborVel)   - *(vel);
      velStateY = *(neighborVel+1) - *(vel+1);
      velStateZ = *(neighborVel+2) - *(vel+2);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(undeformedBondLength, *delta);

      scalarTemp = (1.0 - bondDamage) * omega * neighborVolume;

      *(FdotFirstTerm)   += scalarTemp * velStateX * undeformedBondX;
      *(FdotFirstTerm+1) += scalarTemp * velStateX * undeformedBondY;
      *(FdotFirstTerm+2) += scalarTemp * velStateX * undeformedBondZ;
      *(FdotFirstTerm+3) += scalarTemp * velStateY * undeformedBondX;
      *(FdotFirstTerm+4) += scalarTemp * velStateY * undeformedBondY;
      *(FdotFirstTerm+5) += scalarTemp * velStateY * undeformedBondZ;
      *(FdotFirstTerm+6) += scalarTemp * velStateZ * undeformedBondX;
      *(FdotFirstTerm+7) += scalarTemp * velStateZ * undeformedBondY;
      *(FdotFirstTerm+8) += scalarTemp * velStateZ * undeformedBondZ;
    }

    // Compute Fdot
    MatrixMultiply(false, false, 1.0, FdotFirstTerm, shapeTensorInv, Fdot);

    // Compute the inverse of the deformation gradient, Finverse
    inversionReturnCode = Invert3by3Matrix(defGrad, determinant, Finverse);
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    // Compute the Eulerian velocity gradient L = Fdot * Finv
    MatrixMultiply(false, false, 1.0, Fdot, Finverse, eulerianVelGrad);

    // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
    *(rateOfDef)   = *(eulerianVelGrad);
    *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
    *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
    *(rateOfDef+3) = *(rateOfDef+1);
    *(rateOfDef+4) = *(eulerianVelGrad+4);
    *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
    *(rateOfDef+6) = *(rateOfDef+2);
    *(rateOfDef+7) = *(rateOfDef+5);
    *(rateOfDef+8) = *(eulerianVelGrad+8);

    // Compute spin tensor, W = 1/2 * (L - Lt)
    *(spin)   = 0.0;
    *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
    *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
    *(spin+3) = -1.0 * *(spin+1);
    *(spin+4) = 0.0;
    *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
    *(spin+6) = -1.0 * *(spin+2);
    *(spin+7) = -1.0 * *(spin+5);
    *(spin+8) = 0.0;
   
    //Following Flanagan & Taylor (T&F) 
    //
    //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
    //
    //where \epsilon_{ikj} is the alternator tensor.
    //
    //Components below copied from computer algebra solution to the expansion
    //above
    
    
    zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
           *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
           *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
    zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
           *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
           *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
    zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
           *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
           *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

    //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
    wX = 0.5 * ( *(spin+7) - *(spin+5) );
    wY = 0.5 * ( *(spin+2) - *(spin+6) );
    wZ = 0.5 * ( *(spin+3) - *(spin+1) );

    //Find trace(V)
    traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

    // Compute (trace(V) * I - V) store in temp
    *(temp)   = traceV - *(leftStretchN);
    *(temp+1) = - *(leftStretchN+1);
    *(temp+2) = - *(leftStretchN+2);
    *(temp+3) = - *(leftStretchN+3);
    *(temp+4) = traceV - *(leftStretchN+4);
    *(temp+5) = - *(leftStretchN+5);
    *(temp+6) = - *(leftStretchN+6);
    *(temp+7) = - *(leftStretchN+7);
    *(temp+8) = traceV - *(leftStretchN+8);

    // Compute the inverse of the temp matrix
    Invert3by3Matrix(temp, determinant, tempInv);
    if(inversionReturnCode > 0)
      returnCode = inversionReturnCode;

    //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
    omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
    omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
    omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

    //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
    *(OmegaTensor) = 0.0;
    *(OmegaTensor+1) = -omegaZ;
    *(OmegaTensor+2) = omegaY;
    *(OmegaTensor+3) = omegaZ;
    *(OmegaTensor+4) = 0.0;
    *(OmegaTensor+5) = -omegaX;
    *(OmegaTensor+6) = -omegaY;
    *(OmegaTensor+7) = omegaX;
    *(OmegaTensor+8) = 0.0;

    //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
    //is desirable for accuracy in implicit solves and has no effect on
    //explicit solves (other than a slight decrease in speed).
    //
    // Compute Q with (T&F Eq. 44)
    //
    // Omega^2 = w_i * w_i (T&F Eq. 42)
    OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
    // Omega = \sqrt{OmegaSq}
    Omega = sqrt(OmegaSq);

    // Avoid a potential divide-by-zero
    if( OmegaSq > 1.e-30){

      // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
      //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
      scaleFactor1 = sin(dt*Omega) / Omega;
      scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
      MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
      *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
      *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
      *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
      *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
      *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
      *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
      *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
      *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
      *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

    } else {
      *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
      *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
      *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
    };

    // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
    MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

    // Compute rate of stretch, Vdot = L*V - V*Omega
    // First tempA = L*V, 
    MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

    // tempB = V*Omega
    MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

    //Vdot = tempA - tempB
    for(int i=0 ; i<9 ; ++i)
      *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

    //V_STEP_NP1 = V_STEP_N + dt*Vdot
    for(int i=0 ; i<9 ; ++i)
      *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

    // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
    MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

    // d = Rt * temp
    MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);
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

    *strainXX = 0.5 * ( *(defGradXX) * *(defGradXX) + *(defGradYX) * *(defGradYX) + *(defGradZX) * *(defGradZX) - 1.0 );
    *strainXY = 0.5 * ( *(defGradXX) * *(defGradXY) + *(defGradYX) * *(defGradYY) + *(defGradZX) * *(defGradZY) );
    *strainXZ = 0.5 * ( *(defGradXX) * *(defGradXZ) + *(defGradYX) * *(defGradYZ) + *(defGradZX) * *(defGradZZ) );
    *strainYX = 0.5 * ( *(defGradXY) * *(defGradXX) + *(defGradYY) * *(defGradYX) + *(defGradZY) * *(defGradZX) );
    *strainYY = 0.5 * ( *(defGradXY) * *(defGradXY) + *(defGradYY) * *(defGradYY) + *(defGradZY) * *(defGradZY) - 1.0 );
    *strainYZ = 0.5 * ( *(defGradXY) * *(defGradXZ) + *(defGradYY) * *(defGradYZ) + *(defGradZY) * *(defGradZZ) );
    *strainZX = 0.5 * ( *(defGradXZ) * *(defGradXX) + *(defGradYZ) * *(defGradYX) + *(defGradZZ) * *(defGradZX) );
    *strainZY = 0.5 * ( *(defGradXZ) * *(defGradXY) + *(defGradYZ) * *(defGradYY) + *(defGradZZ) * *(defGradZY) );
    *strainZZ = 0.5 * ( *(defGradXZ) * *(defGradXZ) + *(defGradYZ) * *(defGradYZ) + *(defGradZZ) * *(defGradZZ) - 1.0 );
  }
}

template<typename ScalarT>
void computeHourglassForce
(
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const ScalarT* coordinates,
    const ScalarT* deformationGradient,
    ScalarT* hourglassForceDensity,
    const int* neighborhoodList,
    int numPoints,
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

  const ScalarT* defGrad = deformationGradient;

  const double* delta = horizon;
  const double* modelCoord = modelCoordinates;
  const double* neighborModelCoord;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  ScalarT* hourglassForceDensityPtr = hourglassForceDensity;
  ScalarT* neighborHourglassForceDensityPtr;

  // placeholder for inclusion of bond damage
  double bondDamage = 0.0;

  const double pi = PeridigmNS::value_of_pi();
  double firstPartOfConstant = 18.0*hourglassCoefficient*bulkModulus/pi;
  double constant;

  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, modelCoord+=3, coord+=3,
      defGrad+=9, hourglassForceDensityPtr+=3){

    constant = firstPartOfConstant/( (*delta)*(*delta)*(*delta)*(*delta) );

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
        *(defGrad) * undeformedBondX +
        *(defGrad+1) * undeformedBondY +
        *(defGrad+2) * undeformedBondZ;
      expectedNeighborLocationY = *(coord+1) +
        *(defGrad+3) * undeformedBondX +
        *(defGrad+4) * undeformedBondY +
        *(defGrad+5) * undeformedBondZ;
      expectedNeighborLocationZ = *(coord+2) +
        *(defGrad+6) * undeformedBondX +
        *(defGrad+7) * undeformedBondY +
        *(defGrad+8) * undeformedBondZ;

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
    const ScalarT* rotationTensor,
    const ScalarT* unrotatedCauchyStress,
    ScalarT* rotatedCauchyStress,
    int numPoints
)
{
  const ScalarT* rotTensor = rotationTensor;
  const ScalarT* unrotatedStress = unrotatedCauchyStress;
  ScalarT* rotatedStress = rotatedCauchyStress;
  ScalarT temp[9];

  for(int iID=0 ; iID<numPoints ; ++iID, 
      rotTensor+=9, unrotatedStress+=9, rotatedStress+=9){ 

    // temp = \sigma_unrot * Rt
    CORRESPONDENCE::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
    // \sigma_rot = R * temp
    CORRESPONDENCE::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);
  }
}

template<typename ScalarT>
void computeUndamagedWeightedVolume
(
    const double* volume,
    double* weightedVolume,
    const ScalarT* jacobianDeterminant,
    const double* horizon,
    const ScalarT* coordinates,
    const int* neighborhoodList,
    int numPoints
)
{
  const double* delta = horizon;
  double* w0 = weightedVolume;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  double neighborVolume, omega;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, w0++){

    // Zero out the weighted volume
    *w0 = 0.0;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
      neighborCoord = coordinates + 3*neighborIndex;

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

      *w0 += omega * neighborVolume;
    }
  }
}

template<typename ScalarT>
void computeWeightedVolume
(
    const double* volume,
    double* weightedVolume,
    const ScalarT* jacobianDeterminant,
    const double* horizon,
    const ScalarT* coordinates,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
)
{
  const double* delta = horizon;
  double* w0 = weightedVolume;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength;
  double neighborVolume, omega;

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, w0++, flyingPointFlg++){

    // Zero out the weighted volume
    *w0 = 0.0;

    numNeighbors = *neighborListPtr; neighborListPtr++;
    for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

      neighborIndex = *neighborListPtr;
      neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
      neighborCoord = coordinates + 3*neighborIndex;

      deformedBondX = *(neighborCoord)   - *(coord);
      deformedBondY = *(neighborCoord+1) - *(coord+1);
      deformedBondZ = *(neighborCoord+2) - *(coord+2);
      deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                deformedBondY*deformedBondY +
                                deformedBondZ*deformedBondZ);

      omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

      *w0 += (1.0 - *bondDamagePtr) * omega * neighborVolume;
    }
  }
}

//This function computes the node-level velocity gradient
template<typename ScalarT>
int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient
(
    const double* volume,
    const ScalarT* jacobianDeterminantN,
    ScalarT* jacobianDeterminantNP1,
    const double* horizon,
    const ScalarT* coordinates,
    const ScalarT* velocities,
    ScalarT* shapeTensorInverse,
    ScalarT* velocityGradient,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints,
    double dt
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* velGrad = velocityGradient;
  ScalarT velGradTr;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  ScalarT velStateX, velStateY, velStateZ;
  double neighborVolume, omega, temp;

  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> velGradFirstTermVector(9);
  ScalarT* velGradFirstTerm = &velGradFirstTermVector[0];

  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3,
      vel+=3, shapeTensorInv+=9, velGrad+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Zero out data
      *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
      *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
      *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
      *(velGradFirstTerm)   = 0.0 ; *(velGradFirstTerm+1) = 0.0 ; *(velGradFirstTerm+2) = 0.0 ;
      *(velGradFirstTerm+3) = 0.0 ; *(velGradFirstTerm+4) = 0.0 ; *(velGradFirstTerm+5) = 0.0 ;
      *(velGradFirstTerm+6) = 0.0 ; *(velGradFirstTerm+7) = 0.0 ; *(velGradFirstTerm+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminantN[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // The velState is the relative difference in velocities of the nodes at
        // each end of a bond. i.e., v_j - v_i
        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        temp = (1.0 - *bondDamagePtr) * omega * neighborVolume / (deformedBondLength * deformedBondLength);

        *(shapeTensor)   += temp * deformedBondX * deformedBondX;
        *(shapeTensor+1) += temp * deformedBondX * deformedBondY;
        *(shapeTensor+2) += temp * deformedBondX * deformedBondZ;
        *(shapeTensor+3) += temp * deformedBondY * deformedBondX;
        *(shapeTensor+4) += temp * deformedBondY * deformedBondY;
        *(shapeTensor+5) += temp * deformedBondY * deformedBondZ;
        *(shapeTensor+6) += temp * deformedBondZ * deformedBondX;
        *(shapeTensor+7) += temp * deformedBondZ * deformedBondY;
        *(shapeTensor+8) += temp * deformedBondZ * deformedBondZ;

        *(velGradFirstTerm)   += temp * velStateX * deformedBondX;
        *(velGradFirstTerm+1) += temp * velStateX * deformedBondY;
        *(velGradFirstTerm+2) += temp * velStateX * deformedBondZ;
        *(velGradFirstTerm+3) += temp * velStateY * deformedBondX;
        *(velGradFirstTerm+4) += temp * velStateY * deformedBondY;
        *(velGradFirstTerm+5) += temp * velStateY * deformedBondZ;
        *(velGradFirstTerm+6) += temp * velStateZ * deformedBondX;
        *(velGradFirstTerm+7) += temp * velStateZ * deformedBondY;
        *(velGradFirstTerm+8) += temp * velStateZ * deformedBondZ;
      }
      
      inversionReturnCode = Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
      if(inversionReturnCode > 0){
        returnCode = inversionReturnCode;
        std::cout << inversionErrorMessage;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // Matrix multiply the first term and the shape tensor inverse to compute
      // the velocity gradient
      MatrixMultiply(false, false, 1.0, velGradFirstTerm, shapeTensorInv, velGrad);

      // update volume
      // J dot = J . tr(L)
      velGradTr = *velGrad + *(velGrad+4) + *(velGrad+8);
      jacobianDeterminantNP1[iID] = jacobianDeterminantN[iID] * (1.0 + velGradTr*dt);
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors; 
    }
  }

  return returnCode;
}

//This function computes the node-level velocity gradient
template<typename ScalarT>
int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient
(
    const double* volume,
    const ScalarT* jacobianDeterminantN,
    ScalarT* jacobianDeterminantNP1,
    const double* horizon,
    const ScalarT* coordinates,
    const ScalarT* velocities,
    ScalarT* shapeTensorInverse,
    ScalarT* velocityGradient,
    ScalarT* velocityGradientX,
    ScalarT* velocityGradientY,
    ScalarT* velocityGradientZ,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints,
    double dt
)
{
  int returnCode = 0;

  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  ScalarT* shapeTensorInv = shapeTensorInverse;
  ScalarT* velGrad = velocityGradient;
  ScalarT* velGradX = velocityGradientX;
  ScalarT* velGradY = velocityGradientY;
  ScalarT* velGradZ = velocityGradientZ;
  ScalarT velGradTr;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  ScalarT velStateX, velStateY, velStateZ;
  double neighborVolume, omega, temp;

  std::vector<ScalarT> shapeTensorVector(9);
  ScalarT* shapeTensor = &shapeTensorVector[0];
  ScalarT shapeTensorDeterminant;

  std::vector<ScalarT> velGradFirstTermVector(9);
  ScalarT* velGradFirstTerm = &velGradFirstTermVector[0];

  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, vel+=3, shapeTensorInv+=9, 
      velGrad+=9, velGradX+=3, velGradY+=3, velGradZ+=3, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Zero out data
      *(shapeTensor)   = 0.0 ; *(shapeTensor+1) = 0.0 ; *(shapeTensor+2) = 0.0 ;
      *(shapeTensor+3) = 0.0 ; *(shapeTensor+4) = 0.0 ; *(shapeTensor+5) = 0.0 ;
      *(shapeTensor+6) = 0.0 ; *(shapeTensor+7) = 0.0 ; *(shapeTensor+8) = 0.0 ;
      *(velGradFirstTerm)   = 0.0 ; *(velGradFirstTerm+1) = 0.0 ; *(velGradFirstTerm+2) = 0.0 ;
      *(velGradFirstTerm+3) = 0.0 ; *(velGradFirstTerm+4) = 0.0 ; *(velGradFirstTerm+5) = 0.0 ;
      *(velGradFirstTerm+6) = 0.0 ; *(velGradFirstTerm+7) = 0.0 ; *(velGradFirstTerm+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminantN[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // The velState is the relative difference in velocities of the nodes at
        // each end of a bond. i.e., v_j - v_i
        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        omega = MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        temp = (1.0 - *bondDamagePtr) * omega * neighborVolume / (deformedBondLength * deformedBondLength);

        *(shapeTensor)   += temp * deformedBondX * deformedBondX;
        *(shapeTensor+1) += temp * deformedBondX * deformedBondY;
        *(shapeTensor+2) += temp * deformedBondX * deformedBondZ;
        *(shapeTensor+3) += temp * deformedBondY * deformedBondX;
        *(shapeTensor+4) += temp * deformedBondY * deformedBondY;
        *(shapeTensor+5) += temp * deformedBondY * deformedBondZ;
        *(shapeTensor+6) += temp * deformedBondZ * deformedBondX;
        *(shapeTensor+7) += temp * deformedBondZ * deformedBondY;
        *(shapeTensor+8) += temp * deformedBondZ * deformedBondZ;

        *(velGradFirstTerm)   += temp * velStateX * deformedBondX;
        *(velGradFirstTerm+1) += temp * velStateX * deformedBondY;
        *(velGradFirstTerm+2) += temp * velStateX * deformedBondZ;
        *(velGradFirstTerm+3) += temp * velStateY * deformedBondX;
        *(velGradFirstTerm+4) += temp * velStateY * deformedBondY;
        *(velGradFirstTerm+5) += temp * velStateY * deformedBondZ;
        *(velGradFirstTerm+6) += temp * velStateZ * deformedBondX;
        *(velGradFirstTerm+7) += temp * velStateZ * deformedBondY;
        *(velGradFirstTerm+8) += temp * velStateZ * deformedBondZ;
      }
      
      inversionReturnCode = Invert3by3Matrix(shapeTensor, shapeTensorDeterminant, shapeTensorInv);
      if(inversionReturnCode > 0){
        returnCode = inversionReturnCode;
        std::cout << inversionErrorMessage;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // Matrix multiply the first term and the shape tensor inverse to compute
      // the velocity gradient
      MatrixMultiply(false, false, 1.0, velGradFirstTerm, shapeTensorInv, velGrad);

      *(velGradX+0) = *(velGrad+0); *(velGradX+1) = *(velGrad+1); *(velGradX+2) = *(velGrad+2); 
      *(velGradY+0) = *(velGrad+3); *(velGradY+1) = *(velGrad+4); *(velGradY+2) = *(velGrad+5); 
      *(velGradZ+0) = *(velGrad+6); *(velGradZ+1) = *(velGrad+7); *(velGradZ+2) = *(velGrad+8); 

      // update volume
      // J dot = J . tr(L)
      velGradTr = *velGrad + *(velGrad+4) + *(velGrad+8);
      jacobianDeterminantNP1[iID] = jacobianDeterminantN[iID] * (1.0 + velGradTr*dt);
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors; 
    }
  }

  return returnCode;
}

//This function computes bond-level velocity gradient
template<typename ScalarT>
void computeBondLevelVelocityGradient
(
    const ScalarT* coordinates,
    const ScalarT* velocities,
    const ScalarT* velocityGradient,
    ScalarT* bondLevelVelocityGradientXX,
    ScalarT* bondLevelVelocityGradientXY,
    ScalarT* bondLevelVelocityGradientXZ,
    ScalarT* bondLevelVelocityGradientYX,
    ScalarT* bondLevelVelocityGradientYY,
    ScalarT* bondLevelVelocityGradientYZ,
    ScalarT* bondLevelVelocityGradientZX,
    ScalarT* bondLevelVelocityGradientZY,
    ScalarT* bondLevelVelocityGradientZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
)
{
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* velGrad = velocityGradient;
  const ScalarT* neighborVelGrad;

  ScalarT* bondLevelVelGradXX = bondLevelVelocityGradientXX;
  ScalarT* bondLevelVelGradXY = bondLevelVelocityGradientXY;
  ScalarT* bondLevelVelGradXZ = bondLevelVelocityGradientXZ;
  ScalarT* bondLevelVelGradYX = bondLevelVelocityGradientYX;
  ScalarT* bondLevelVelGradYY = bondLevelVelocityGradientYY;
  ScalarT* bondLevelVelGradYZ = bondLevelVelocityGradientYZ;
  ScalarT* bondLevelVelGradZX = bondLevelVelocityGradientZX;
  ScalarT* bondLevelVelGradZY = bondLevelVelocityGradientZY;
  ScalarT* bondLevelVelGradZZ = bondLevelVelocityGradientZZ;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT temp;

  std::vector<ScalarT> meanVelGradVector(9);
  ScalarT* meanVelGrad = &meanVelGradVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, coord+=3, vel+=3, velGrad+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          bondLevelVelGradXX++, bondLevelVelGradXY++, bondLevelVelGradXZ++, 
          bondLevelVelGradYX++, bondLevelVelGradYY++, bondLevelVelGradYZ++,
          bondLevelVelGradZX++, bondLevelVelGradZY++, bondLevelVelGradZZ++){

        neighborIndex = *neighborListPtr;
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;
        neighborVelGrad = velocityGradient + 9*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // The velState is the relative difference in velocities of the nodes at
        // each end of a bond. i.e., v_j - v_i
        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        // average of the two points 
        for(int i=0; i<9; i++)
          *(meanVelGrad+i) = 0.5 * (*(velGrad+i) + *(neighborVelGrad+i));

        deformedBondLengthSq = deformedBondLength * deformedBondLength;

        temp = *(meanVelGrad+0) * deformedBondX + *(meanVelGrad+1) * deformedBondY + *(meanVelGrad+2) * deformedBondZ;
        *bondLevelVelGradXX = *(meanVelGrad+0) + (velStateX - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradXY = *(meanVelGrad+1) + (velStateX - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradXZ = *(meanVelGrad+2) + (velStateX - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+3) * deformedBondX + *(meanVelGrad+4) * deformedBondY + *(meanVelGrad+5) * deformedBondZ;
        *bondLevelVelGradYX = *(meanVelGrad+3) + (velStateY - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradYY = *(meanVelGrad+4) + (velStateY - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradYZ = *(meanVelGrad+5) + (velStateY - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+6) * deformedBondX + *(meanVelGrad+7) * deformedBondY + *(meanVelGrad+8) * deformedBondZ;
        *bondLevelVelGradZX = *(meanVelGrad+6) + (velStateZ - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradZY = *(meanVelGrad+7) + (velStateZ - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradZZ = *(meanVelGrad+8) + (velStateZ - temp) * deformedBondZ/deformedBondLengthSq;
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondLevelVelGradXX += numNeighbors; bondLevelVelGradXY += numNeighbors; bondLevelVelGradXZ += numNeighbors; 
      bondLevelVelGradYX += numNeighbors; bondLevelVelGradYY += numNeighbors; bondLevelVelGradYZ += numNeighbors;
      bondLevelVelGradZX += numNeighbors; bondLevelVelGradZY += numNeighbors; bondLevelVelGradZZ += numNeighbors;
    }
  }
}

//This function computes bond-level velocity gradient
template<typename ScalarT>
void computeBondLevelVelocityGradient
(
    const ScalarT* coordinates,
    const ScalarT* velocities,
    const ScalarT* velocityGradientX,
    const ScalarT* velocityGradientY,
    const ScalarT* velocityGradientZ,
    ScalarT* bondLevelVelocityGradientXX,
    ScalarT* bondLevelVelocityGradientXY,
    ScalarT* bondLevelVelocityGradientXZ,
    ScalarT* bondLevelVelocityGradientYX,
    ScalarT* bondLevelVelocityGradientYY,
    ScalarT* bondLevelVelocityGradientYZ,
    ScalarT* bondLevelVelocityGradientZX,
    ScalarT* bondLevelVelocityGradientZY,
    ScalarT* bondLevelVelocityGradientZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
)
{
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const ScalarT* vel = velocities;
  const ScalarT* neighborVel;
  const ScalarT* velGradX = velocityGradientX;
  const ScalarT* velGradY = velocityGradientY;
  const ScalarT* velGradZ = velocityGradientZ;
  const ScalarT* neighborVelGradX;
  const ScalarT* neighborVelGradY;
  const ScalarT* neighborVelGradZ;

  ScalarT* bondLevelVelGradXX = bondLevelVelocityGradientXX;
  ScalarT* bondLevelVelGradXY = bondLevelVelocityGradientXY;
  ScalarT* bondLevelVelGradXZ = bondLevelVelocityGradientXZ;
  ScalarT* bondLevelVelGradYX = bondLevelVelocityGradientYX;
  ScalarT* bondLevelVelGradYY = bondLevelVelocityGradientYY;
  ScalarT* bondLevelVelGradYZ = bondLevelVelocityGradientYZ;
  ScalarT* bondLevelVelGradZX = bondLevelVelocityGradientZX;
  ScalarT* bondLevelVelGradZY = bondLevelVelocityGradientZY;
  ScalarT* bondLevelVelGradZZ = bondLevelVelocityGradientZZ;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  ScalarT velStateX, velStateY, velStateZ;
  ScalarT temp;

  std::vector<ScalarT> meanVelGradVector(9);
  ScalarT* meanVelGrad = &meanVelGradVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, coord+=3, vel+=3, 
      velGradX+=3,  velGradY+=3, velGradZ+=3, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          bondLevelVelGradXX++, bondLevelVelGradXY++, bondLevelVelGradXZ++, 
          bondLevelVelGradYX++, bondLevelVelGradYY++, bondLevelVelGradYZ++,
          bondLevelVelGradZX++, bondLevelVelGradZY++, bondLevelVelGradZZ++){

        neighborIndex = *neighborListPtr;
        neighborCoord = coordinates + 3*neighborIndex;
        neighborVel = velocities + 3*neighborIndex;
        neighborVelGradX = velocityGradientX + 3*neighborIndex;
        neighborVelGradY = velocityGradientY + 3*neighborIndex;
        neighborVelGradZ = velocityGradientZ + 3*neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);

        // The velState is the relative difference in velocities of the nodes at
        // each end of a bond. i.e., v_j - v_i
        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        // average of the two points 
        for(int i=0; i<3; i++){
          *(meanVelGrad+i) = 0.5 * (*(velGradX+i) + *(neighborVelGradX+i));
          *(meanVelGrad+i+3) = 0.5 * (*(velGradY+i) + *(neighborVelGradY+i));
          *(meanVelGrad+i+6) = 0.5 * (*(velGradZ+i) + *(neighborVelGradZ+i));
        }

        deformedBondLengthSq = deformedBondLength * deformedBondLength;

        velStateX = *(neighborVel)   - *(vel);
        velStateY = *(neighborVel+1) - *(vel+1);
        velStateZ = *(neighborVel+2) - *(vel+2);

        temp = *(meanVelGrad+0) * deformedBondX + *(meanVelGrad+1) * deformedBondY + *(meanVelGrad+2) * deformedBondZ;
        *bondLevelVelGradXX = *(meanVelGrad+0) + (velStateX - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradXY = *(meanVelGrad+1) + (velStateX - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradXZ = *(meanVelGrad+2) + (velStateX - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+3) * deformedBondX + *(meanVelGrad+4) * deformedBondY + *(meanVelGrad+5) * deformedBondZ;
        *bondLevelVelGradYX = *(meanVelGrad+3) + (velStateY - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradYY = *(meanVelGrad+4) + (velStateY - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradYZ = *(meanVelGrad+5) + (velStateY - temp) * deformedBondZ/deformedBondLengthSq;

        temp = *(meanVelGrad+6) * deformedBondX + *(meanVelGrad+7) * deformedBondY + *(meanVelGrad+8) * deformedBondZ;
        *bondLevelVelGradZX = *(meanVelGrad+6) + (velStateZ - temp) * deformedBondX/deformedBondLengthSq;
        *bondLevelVelGradZY = *(meanVelGrad+7) + (velStateZ - temp) * deformedBondY/deformedBondLengthSq;
        *bondLevelVelGradZZ = *(meanVelGrad+8) + (velStateZ - temp) * deformedBondZ/deformedBondLengthSq;
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondLevelVelGradXX += numNeighbors; bondLevelVelGradXY += numNeighbors; bondLevelVelGradXZ += numNeighbors; 
      bondLevelVelGradYX += numNeighbors; bondLevelVelGradYY += numNeighbors; bondLevelVelGradYZ += numNeighbors;
      bondLevelVelGradZX += numNeighbors; bondLevelVelGradZY += numNeighbors; bondLevelVelGradZZ += numNeighbors;
    }
  }
}

//This function updates the node-level deformation gradient based on velocity gradient
template<typename ScalarT>
void updateDeformationGradient
(
    const ScalarT* velocityGradient,
    const ScalarT* deformationGradientN,
    ScalarT* deformationGradientNP1,
    const double* flyingPointFlag,
    int numPoints,
    double dt
)
{
  const ScalarT* velGrad = velocityGradient;
  const ScalarT* defGradN = deformationGradientN;
  ScalarT* defGradNP1 = deformationGradientNP1;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> FdotVector(9);
  ScalarT* Fdot = &FdotVector[0];

  for(int iID=0 ; iID<numPoints ; ++iID, velGrad+=9, defGradN+=9, defGradNP1+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // L = Fdot . inv(F)
      // Fdot = L . F
      // F_NP1 = F_N + Fdot . dt
      for(int i=0; i<9; i++)
        *(defGradNP1+i) = *(defGradN+i);

      MatrixMultiply(false, false, 1.0, velGrad, defGradN, Fdot);

      for(int i=0; i<9; i++)
        *(defGradNP1+i) += *(Fdot+i) * dt;
    }
  }
}

template<typename ScalarT>
void computeGreenLagrangeStrain
(
    const ScalarT* deformationGradient,
    ScalarT* greenLagrangeStrain,
    const double* flyingPointFlag,
    int numPoints
)
{
  // Green-Lagrange Strain E = 0.5*(F^T F - I)

  const ScalarT* defGrad = deformationGradient;
  ScalarT* strain = greenLagrangeStrain;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> tempVector(9);
  ScalarT* temp = &tempVector[0];

  for(int iID=0 ; iID<numPoints ; ++iID, defGrad+=9, strain+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      MatrixMultiply(true, false, 1.0, defGrad, defGrad, temp);

      for(int i=0; i<9; i++)
        *(strain+i) = 0.5 * *(temp+i);

      *(strain+0) -= 0.5;
      *(strain+4) -= 0.5;
      *(strain+8) -= 0.5;
    }
  }
}

//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
//This function computes the node-level values
template<typename ScalarT>
int computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor(
    const ScalarT* velocityGradient,
    const ScalarT* leftStretchTensorN,
    const ScalarT* rotationTensorN,
    ScalarT* leftStretchTensorNP1,
    ScalarT* rotationTensorNP1,
    ScalarT* unrotatedRateOfDeformation,
    const double* flyingPointFlag,
    int numPoints,
    double dt
)
{
  int returnCode = 0;

  const ScalarT* eulerianVelGrad = velocityGradient;
  const ScalarT* leftStretchN = leftStretchTensorN;
  const ScalarT* rotTensorN = rotationTensorN;

  ScalarT* leftStretchNP1 = leftStretchTensorNP1;
  ScalarT* rotTensorNP1 = rotationTensorNP1;
  ScalarT* unrotRateOfDef = unrotatedRateOfDeformation;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> rateOfDefVector(9) ; ScalarT* rateOfDef = &rateOfDefVector[0];
  std::vector<ScalarT> spinVector(9) ; ScalarT* spin = &spinVector[0];
  std::vector<ScalarT> tempVector(9) ; ScalarT* temp = &tempVector[0];
  std::vector<ScalarT> tempInvVector(9) ; ScalarT* tempInv = &tempInvVector[0];
  std::vector<ScalarT> OmegaTensorVector(9) ; ScalarT* OmegaTensor = &OmegaTensorVector[0];
  std::vector<ScalarT> QMatrixVector(9) ; ScalarT* QMatrix = &QMatrixVector[0];
  std::vector<ScalarT> OmegaTensorSqVector(9) ; ScalarT* OmegaTensorSq = &OmegaTensorSqVector[0];
  std::vector<ScalarT> tempAVector(9) ; ScalarT* tempA = &tempAVector[0];
  std::vector<ScalarT> tempBVector(9) ; ScalarT* tempB = &tempBVector[0];
  std::vector<ScalarT> rateOfStretchVector(9) ; ScalarT* rateOfStretch = &rateOfStretchVector[0];

  ScalarT determinant;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  for(int iID=0 ; iID<numPoints ; ++iID, eulerianVelGrad+=9, rotTensorN+=9, 
      rotTensorNP1+=9, leftStretchNP1+=9, leftStretchN+=9, unrotRateOfDef+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
      *(rateOfDef)   = *(eulerianVelGrad);
      *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
      *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
      *(rateOfDef+3) = *(rateOfDef+1);
      *(rateOfDef+4) = *(eulerianVelGrad+4);
      *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
      *(rateOfDef+6) = *(rateOfDef+2);
      *(rateOfDef+7) = *(rateOfDef+5);
      *(rateOfDef+8) = *(eulerianVelGrad+8);

      // Compute spin tensor, W = 1/2 * (L - Lt)
      *(spin)   = 0.0;
      *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
      *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
      *(spin+3) = -1.0 * *(spin+1);
      *(spin+4) = 0.0;
      *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
      *(spin+6) = -1.0 * *(spin+2);
      *(spin+7) = -1.0 * *(spin+5);
      *(spin+8) = 0.0;
     
      //Following Flanagan & Taylor (T&F) 
      //
      //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
      //
      //where \epsilon_{ikj} is the alternator tensor.
      //
      //Components below copied from computer algebra solution to the expansion
      //above
      
      zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
             *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
             *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
      zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
             *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
             *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
      zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
             *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
             *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

      //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
      wX = 0.5 * ( *(spin+7) - *(spin+5) );
      wY = 0.5 * ( *(spin+2) - *(spin+6) );
      wZ = 0.5 * ( *(spin+3) - *(spin+1) );

      //Find trace(V)
      traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

      // Compute (trace(V) * I - V) store in temp
      *(temp)   = traceV - *(leftStretchN);
      *(temp+1) = - *(leftStretchN+1);
      *(temp+2) = - *(leftStretchN+2);
      *(temp+3) = - *(leftStretchN+3);
      *(temp+4) = traceV - *(leftStretchN+4);
      *(temp+5) = - *(leftStretchN+5);
      *(temp+6) = - *(leftStretchN+6);
      *(temp+7) = - *(leftStretchN+7);
      *(temp+8) = traceV - *(leftStretchN+8);

      // Compute the inverse of the temp matrix
      Invert3by3Matrix(temp, determinant, tempInv);
      if(inversionReturnCode > 0){
        returnCode = inversionReturnCode;
        std::cout << inversionErrorMessage;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
      omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
      omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
      omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

      //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
      *(OmegaTensor) = 0.0;
      *(OmegaTensor+1) = -omegaZ;
      *(OmegaTensor+2) = omegaY;
      *(OmegaTensor+3) = omegaZ;
      *(OmegaTensor+4) = 0.0;
      *(OmegaTensor+5) = -omegaX;
      *(OmegaTensor+6) = -omegaY;
      *(OmegaTensor+7) = omegaX;
      *(OmegaTensor+8) = 0.0;

      //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
      //is desirable for accuracy in implicit solves and has no effect on
      //explicit solves (other than a slight decrease in speed).
      //
      // Compute Q with (T&F Eq. 44)
      //
      // Omega^2 = w_i * w_i (T&F Eq. 42)
      OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
      // Omega = \sqrt{OmegaSq}
      Omega = sqrt(OmegaSq);

      // Avoid a potential divide-by-zero
      if(OmegaSq > 1.e-30){

        // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
        //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
        scaleFactor1 = sin(dt*Omega) / Omega;
        scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
        MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
        *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
        *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
        *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
        *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
        *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
        *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
        *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
        *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
        *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

      } else {
        *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
        *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
        *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
      };

      // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
      MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

      // Compute rate of stretch, Vdot = L*V - V*Omega
      // First tempA = L*V, 
      MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

      // tempB = V*Omega
      MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

      //Vdot = tempA - tempB
      for(int i=0 ; i<9 ; ++i)
        *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

      //V_STEP_NP1 = V_STEP_N + dt*Vdot
      for(int i=0 ; i<9 ; ++i)
        *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

      // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
      MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

      // d = Rt * temp
      MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);
    }
  }

  return returnCode;
}


//Performs kinematic computations following Flanagan and Taylor (1987), returns
//unrotated rate-of-deformation and rotation tensors
//This function computes the node-based values
template<typename ScalarT>
int computeBondLevelUnrotatedRateOfDeformationAndRotationTensor(
    const ScalarT* bondLevelVelocityGradientXX, 
    const ScalarT* bondLevelVelocityGradientXY, 
    const ScalarT* bondLevelVelocityGradientXZ,
    const ScalarT* bondLevelVelocityGradientYX, 
    const ScalarT* bondLevelVelocityGradientYY, 
    const ScalarT* bondLevelVelocityGradientYZ, 
    const ScalarT* bondLevelVelocityGradientZX,
    const ScalarT* bondLevelVelocityGradientZY,
    const ScalarT* bondLevelVelocityGradientZZ,
    const ScalarT* bondLevelLeftStretchTensorXXN,
    const ScalarT* bondLevelLeftStretchTensorXYN,
    const ScalarT* bondLevelLeftStretchTensorXZN,
    const ScalarT* bondLevelLeftStretchTensorYXN,
    const ScalarT* bondLevelLeftStretchTensorYYN,
    const ScalarT* bondLevelLeftStretchTensorYZN,
    const ScalarT* bondLevelLeftStretchTensorZXN,
    const ScalarT* bondLevelLeftStretchTensorZYN,
    const ScalarT* bondLevelLeftStretchTensorZZN,
    const ScalarT* bondLevelRotationTensorXXN, 
    const ScalarT* bondLevelRotationTensorXYN, 
    const ScalarT* bondLevelRotationTensorXZN, 
    const ScalarT* bondLevelRotationTensorYXN, 
    const ScalarT* bondLevelRotationTensorYYN, 
    const ScalarT* bondLevelRotationTensorYZN, 
    const ScalarT* bondLevelRotationTensorZXN, 
    const ScalarT* bondLevelRotationTensorZYN, 
    const ScalarT* bondLevelRotationTensorZZN, 
    ScalarT* bondLevelLeftStretchTensorXXNP1,
    ScalarT* bondLevelLeftStretchTensorXYNP1,
    ScalarT* bondLevelLeftStretchTensorXZNP1,
    ScalarT* bondLevelLeftStretchTensorYXNP1,
    ScalarT* bondLevelLeftStretchTensorYYNP1,
    ScalarT* bondLevelLeftStretchTensorYZNP1,
    ScalarT* bondLevelLeftStretchTensorZXNP1,
    ScalarT* bondLevelLeftStretchTensorZYNP1,
    ScalarT* bondLevelLeftStretchTensorZZNP1,
    ScalarT* bondLevelRotationTensorXXNP1,
    ScalarT* bondLevelRotationTensorXYNP1,
    ScalarT* bondLevelRotationTensorXZNP1,
    ScalarT* bondLevelRotationTensorYXNP1,
    ScalarT* bondLevelRotationTensorYYNP1,
    ScalarT* bondLevelRotationTensorYZNP1,
    ScalarT* bondLevelRotationTensorZXNP1,
    ScalarT* bondLevelRotationTensorZYNP1,
    ScalarT* bondLevelRotationTensorZZNP1,
    ScalarT* bondLevelUnrotatedRateOfDeformationXX,
    ScalarT* bondLevelUnrotatedRateOfDeformationXY,
    ScalarT* bondLevelUnrotatedRateOfDeformationXZ,
    ScalarT* bondLevelUnrotatedRateOfDeformationYX,
    ScalarT* bondLevelUnrotatedRateOfDeformationYY,
    ScalarT* bondLevelUnrotatedRateOfDeformationYZ,
    ScalarT* bondLevelUnrotatedRateOfDeformationZX,
    ScalarT* bondLevelUnrotatedRateOfDeformationZY,
    ScalarT* bondLevelUnrotatedRateOfDeformationZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints,
    double dt
)
{
  int returnCode = 0;

  const ScalarT* velGradXX = bondLevelVelocityGradientXX;
  const ScalarT* velGradXY = bondLevelVelocityGradientXY;
  const ScalarT* velGradXZ = bondLevelVelocityGradientXZ;
  const ScalarT* velGradYX = bondLevelVelocityGradientYX;
  const ScalarT* velGradYY = bondLevelVelocityGradientYY;
  const ScalarT* velGradYZ = bondLevelVelocityGradientYZ;
  const ScalarT* velGradZX = bondLevelVelocityGradientZX;
  const ScalarT* velGradZY = bondLevelVelocityGradientZY;
  const ScalarT* velGradZZ = bondLevelVelocityGradientZZ;
  const ScalarT* leftStretchXXN = bondLevelLeftStretchTensorXXN;
  const ScalarT* leftStretchXYN = bondLevelLeftStretchTensorXYN;
  const ScalarT* leftStretchXZN = bondLevelLeftStretchTensorXZN;
  const ScalarT* leftStretchYXN = bondLevelLeftStretchTensorYXN;
  const ScalarT* leftStretchYYN = bondLevelLeftStretchTensorYYN;
  const ScalarT* leftStretchYZN = bondLevelLeftStretchTensorYZN;
  const ScalarT* leftStretchZXN = bondLevelLeftStretchTensorZXN;
  const ScalarT* leftStretchZYN = bondLevelLeftStretchTensorZYN;
  const ScalarT* leftStretchZZN = bondLevelLeftStretchTensorZZN;
  const ScalarT* rotTensorXXN = bondLevelRotationTensorXXN;
  const ScalarT* rotTensorXYN = bondLevelRotationTensorXYN;
  const ScalarT* rotTensorXZN = bondLevelRotationTensorXZN;
  const ScalarT* rotTensorYXN = bondLevelRotationTensorYXN;
  const ScalarT* rotTensorYYN = bondLevelRotationTensorYYN;
  const ScalarT* rotTensorYZN = bondLevelRotationTensorYZN;
  const ScalarT* rotTensorZXN = bondLevelRotationTensorZXN;
  const ScalarT* rotTensorZYN = bondLevelRotationTensorZYN;
  const ScalarT* rotTensorZZN = bondLevelRotationTensorZZN;

  ScalarT* leftStretchXXNP1 = bondLevelLeftStretchTensorXXNP1;
  ScalarT* leftStretchXYNP1 = bondLevelLeftStretchTensorXYNP1;
  ScalarT* leftStretchXZNP1 = bondLevelLeftStretchTensorXZNP1;
  ScalarT* leftStretchYXNP1 = bondLevelLeftStretchTensorYXNP1;
  ScalarT* leftStretchYYNP1 = bondLevelLeftStretchTensorYYNP1;
  ScalarT* leftStretchYZNP1 = bondLevelLeftStretchTensorYZNP1;
  ScalarT* leftStretchZXNP1 = bondLevelLeftStretchTensorZXNP1;
  ScalarT* leftStretchZYNP1 = bondLevelLeftStretchTensorZYNP1;
  ScalarT* leftStretchZZNP1 = bondLevelLeftStretchTensorZZNP1;
  ScalarT* rotTensorXXNP1 = bondLevelRotationTensorXXNP1;
  ScalarT* rotTensorXYNP1 = bondLevelRotationTensorXYNP1;
  ScalarT* rotTensorXZNP1 = bondLevelRotationTensorXZNP1;
  ScalarT* rotTensorYXNP1 = bondLevelRotationTensorYXNP1;
  ScalarT* rotTensorYYNP1 = bondLevelRotationTensorYYNP1;
  ScalarT* rotTensorYZNP1 = bondLevelRotationTensorYZNP1;
  ScalarT* rotTensorZXNP1 = bondLevelRotationTensorZXNP1;
  ScalarT* rotTensorZYNP1 = bondLevelRotationTensorZYNP1;
  ScalarT* rotTensorZZNP1 = bondLevelRotationTensorZZNP1;
  ScalarT* unrotRateOfDefXX = bondLevelUnrotatedRateOfDeformationXX;
  ScalarT* unrotRateOfDefXY = bondLevelUnrotatedRateOfDeformationXY;
  ScalarT* unrotRateOfDefXZ = bondLevelUnrotatedRateOfDeformationXZ;
  ScalarT* unrotRateOfDefYX = bondLevelUnrotatedRateOfDeformationYX;
  ScalarT* unrotRateOfDefYY = bondLevelUnrotatedRateOfDeformationYY;
  ScalarT* unrotRateOfDefYZ = bondLevelUnrotatedRateOfDeformationYZ;
  ScalarT* unrotRateOfDefZX = bondLevelUnrotatedRateOfDeformationZX;
  ScalarT* unrotRateOfDefZY = bondLevelUnrotatedRateOfDeformationZY;
  ScalarT* unrotRateOfDefZZ = bondLevelUnrotatedRateOfDeformationZZ;

  const double* flyingPointFlg = flyingPointFlag;

  std::vector<ScalarT> eulerianVelGradVector(9) ; ScalarT* eulerianVelGrad = &eulerianVelGradVector[0];
  std::vector<ScalarT> leftStretchNVector(9) ; ScalarT* leftStretchN = &leftStretchNVector[0];
  std::vector<ScalarT> leftStretchNP1Vector(9) ; ScalarT* leftStretchNP1 = &leftStretchNP1Vector[0];
  std::vector<ScalarT> rotTensorNVector(9) ; ScalarT* rotTensorN = &rotTensorNVector[0];
  std::vector<ScalarT> rotTensorNP1Vector(9) ; ScalarT* rotTensorNP1 = &rotTensorNP1Vector[0];
  std::vector<ScalarT> unrotRateOfDefVector(9) ; ScalarT* unrotRateOfDef = &unrotRateOfDefVector[0];

  std::vector<ScalarT> rateOfDefVector(9) ; ScalarT* rateOfDef = &rateOfDefVector[0];
  std::vector<ScalarT> spinVector(9) ; ScalarT* spin = &spinVector[0];
  std::vector<ScalarT> tempVector(9) ; ScalarT* temp = &tempVector[0];
  std::vector<ScalarT> tempInvVector(9) ; ScalarT* tempInv = &tempInvVector[0];
  std::vector<ScalarT> OmegaTensorVector(9) ; ScalarT* OmegaTensor = &OmegaTensorVector[0];
  std::vector<ScalarT> QMatrixVector(9) ; ScalarT* QMatrix = &QMatrixVector[0];
  std::vector<ScalarT> OmegaTensorSqVector(9) ; ScalarT* OmegaTensorSq = &OmegaTensorSqVector[0];
  std::vector<ScalarT> tempAVector(9) ; ScalarT* tempA = &tempAVector[0];
  std::vector<ScalarT> tempBVector(9) ; ScalarT* tempB = &tempBVector[0];
  std::vector<ScalarT> rateOfStretchVector(9) ; ScalarT* rateOfStretch = &rateOfStretchVector[0];

  ScalarT determinant;
  ScalarT omegaX, omegaY, omegaZ;
  ScalarT zX, zY, zZ;
  ScalarT wX, wY, wZ;
  ScalarT traceV, Omega, OmegaSq, scaleFactor1, scaleFactor2;
  int inversionReturnCode(0);
  std::string inversionErrorMessage = 
    "**** Error:  computeShapeTensorInverseAndApproximateVelocityGradient: Non-invertible matrix\n";

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // All is bond level.
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          velGradXX++, velGradXY++, velGradXZ++, 
          velGradYX++, velGradYY++, velGradYZ++, 
          velGradZX++, velGradZY++, velGradZZ++,
          leftStretchXXN++, leftStretchXYN++, leftStretchXZN++, 
          leftStretchYXN++, leftStretchYYN++, leftStretchYZN++, 
          leftStretchZXN++, leftStretchZYN++, leftStretchZZN++,
          rotTensorXXN++, rotTensorXYN++, rotTensorXZN++,
          rotTensorYXN++, rotTensorYYN++, rotTensorYZN++,
          rotTensorZXN++, rotTensorZYN++, rotTensorZZN++,
          leftStretchXXNP1++, leftStretchXYNP1++, leftStretchXZNP1++, 
          leftStretchYXNP1++, leftStretchYYNP1++, leftStretchYZNP1++, 
          leftStretchZXNP1++, leftStretchZYNP1++, leftStretchZZNP1++,
          rotTensorXXNP1++, rotTensorXYNP1++, rotTensorXZNP1++,
          rotTensorYXNP1++, rotTensorYYNP1++, rotTensorYZNP1++,
          rotTensorZXNP1++, rotTensorZYNP1++, rotTensorZZNP1++,
          unrotRateOfDefXX++, unrotRateOfDefXY++, unrotRateOfDefXZ++,
          unrotRateOfDefYX++, unrotRateOfDefYY++, unrotRateOfDefYZ++,
          unrotRateOfDefZX++, unrotRateOfDefZY++, unrotRateOfDefZZ++){

        neighborIndex = *neighborListPtr;

        // Store in a tensor form 
        *(eulerianVelGrad+0) = *velGradXX; *(eulerianVelGrad+1) = *velGradXY; *(eulerianVelGrad+2) = *velGradXZ;
        *(eulerianVelGrad+3) = *velGradYX; *(eulerianVelGrad+4) = *velGradYY; *(eulerianVelGrad+5) = *velGradYZ;
        *(eulerianVelGrad+6) = *velGradZX; *(eulerianVelGrad+7) = *velGradZY; *(eulerianVelGrad+8) = *velGradZZ;
        *(leftStretchN+0) = *leftStretchXXN; *(leftStretchN+1) = *leftStretchXYN; *(leftStretchN+2) = *leftStretchXZN;
        *(leftStretchN+3) = *leftStretchYXN; *(leftStretchN+4) = *leftStretchYYN; *(leftStretchN+5) = *leftStretchYZN;
        *(leftStretchN+6) = *leftStretchZXN; *(leftStretchN+7) = *leftStretchZYN; *(leftStretchN+8) = *leftStretchZZN;
        *(rotTensorN+0) = *rotTensorXXN; *(rotTensorN+1) = *rotTensorXYN; *(rotTensorN+2) = *rotTensorXZN;
        *(rotTensorN+3) = *rotTensorYXN; *(rotTensorN+4) = *rotTensorYYN; *(rotTensorN+5) = *rotTensorYZN;
        *(rotTensorN+6) = *rotTensorZXN; *(rotTensorN+7) = *rotTensorZYN; *(rotTensorN+8) = *rotTensorZZN;

        // Compute rate-of-deformation tensor, D = 1/2 * (L + Lt)
        *(rateOfDef)   = *(eulerianVelGrad);
        *(rateOfDef+1) = 0.5 * ( *(eulerianVelGrad+1) + *(eulerianVelGrad+3) );
        *(rateOfDef+2) = 0.5 * ( *(eulerianVelGrad+2) + *(eulerianVelGrad+6) );
        *(rateOfDef+3) = *(rateOfDef+1);
        *(rateOfDef+4) = *(eulerianVelGrad+4);
        *(rateOfDef+5) = 0.5 * ( *(eulerianVelGrad+5) + *(eulerianVelGrad+7) );
        *(rateOfDef+6) = *(rateOfDef+2);
        *(rateOfDef+7) = *(rateOfDef+5);
        *(rateOfDef+8) = *(eulerianVelGrad+8);

        // Compute spin tensor, W = 1/2 * (L - Lt)
        *(spin)   = 0.0;
        *(spin+1) = 0.5 * ( *(eulerianVelGrad+1) - *(eulerianVelGrad+3) );
        *(spin+2) = 0.5 * ( *(eulerianVelGrad+2) - *(eulerianVelGrad+6) );
        *(spin+3) = -1.0 * *(spin+1);
        *(spin+4) = 0.0;
        *(spin+5) = 0.5 * ( *(eulerianVelGrad+5) - *(eulerianVelGrad+7) );
        *(spin+6) = -1.0 * *(spin+2);
        *(spin+7) = -1.0 * *(spin+5);
        *(spin+8) = 0.0;
       
        //Following Flanagan & Taylor (T&F) 
        //
        //Find the vector z_i = \epsilon_{ikj} * D_{jm} * V_{mk} (T&F Eq. 13)
        //
        //where \epsilon_{ikj} is the alternator tensor.
        //
        //Components below copied from computer algebra solution to the expansion
        //above
        zX = - *(leftStretchN+2) *  *(rateOfDef+3) -  *(leftStretchN+5) *  *(rateOfDef+4) - 
               *(leftStretchN+8) *  *(rateOfDef+5) +  *(leftStretchN+1) *  *(rateOfDef+6) + 
               *(leftStretchN+4) *  *(rateOfDef+7) +  *(leftStretchN+7) *  *(rateOfDef+8);
        zY =   *(leftStretchN+2) *  *(rateOfDef)   +  *(leftStretchN+5) *  *(rateOfDef+1) + 
               *(leftStretchN+8) *  *(rateOfDef+2) -  *(leftStretchN)   *  *(rateOfDef+6) - 
               *(leftStretchN+3) *  *(rateOfDef+7) -  *(leftStretchN+6) *  *(rateOfDef+8);
        zZ = - *(leftStretchN+1) *  *(rateOfDef)   -  *(leftStretchN+4) *  *(rateOfDef+1) - 
               *(leftStretchN+7) *  *(rateOfDef+2) +  *(leftStretchN)   *  *(rateOfDef+3) + 
               *(leftStretchN+3) *  *(rateOfDef+4) +  *(leftStretchN+6) *  *(rateOfDef+5);

        //Find the vector w_i = -1/2 * \epsilon_{ijk} * W_{jk} (T&F Eq. 11)
        wX = 0.5 * ( *(spin+7) - *(spin+5) );
        wY = 0.5 * ( *(spin+2) - *(spin+6) );
        wZ = 0.5 * ( *(spin+3) - *(spin+1) );

        //Find trace(V)
        traceV = *(leftStretchN) + *(leftStretchN+4) + *(leftStretchN+8);

        // Compute (trace(V) * I - V) store in temp
        *(temp)   = traceV - *(leftStretchN);
        *(temp+1) = - *(leftStretchN+1);
        *(temp+2) = - *(leftStretchN+2);
        *(temp+3) = - *(leftStretchN+3);
        *(temp+4) = traceV - *(leftStretchN+4);
        *(temp+5) = - *(leftStretchN+5);
        *(temp+6) = - *(leftStretchN+6);
        *(temp+7) = - *(leftStretchN+7);
        *(temp+8) = traceV - *(leftStretchN+8);

        // Compute the inverse of the temp matrix
        Invert3by3Matrix(temp, determinant, tempInv);
        if(inversionReturnCode > 0){
          returnCode = inversionReturnCode;
          std::cout << inversionErrorMessage;
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

        //Find omega vector, i.e. \omega = w +  (trace(V) I - V)^(-1) * z (T&F Eq. 12)
        omegaX =  wX + *(tempInv)   * zX + *(tempInv+1) * zY + *(tempInv+2) * zZ;
        omegaY =  wY + *(tempInv+3) * zX + *(tempInv+4) * zY + *(tempInv+5) * zZ;
        omegaZ =  wZ + *(tempInv+6) * zX + *(tempInv+7) * zY + *(tempInv+8) * zZ;

        //Find the tensor \Omega_{ij} = \epsilon_{ikj} * w_k (T&F Eq. 10)
        *(OmegaTensor) = 0.0;
        *(OmegaTensor+1) = -omegaZ;
        *(OmegaTensor+2) = omegaY;
        *(OmegaTensor+3) = omegaZ;
        *(OmegaTensor+4) = 0.0;
        *(OmegaTensor+5) = -omegaX;
        *(OmegaTensor+6) = -omegaY;
        *(OmegaTensor+7) = omegaX;
        *(OmegaTensor+8) = 0.0;

        //Increment R with (T&F Eq. 36 and 44) as opposed to solving (T&F 39) this
        //is desirable for accuracy in implicit solves and has no effect on
        //explicit solves (other than a slight decrease in speed).
        //
        // Compute Q with (T&F Eq. 44)
        //
        // Omega^2 = w_i * w_i (T&F Eq. 42)
        OmegaSq = omegaX*omegaX + omegaY*omegaY + omegaZ*omegaZ;
        // Omega = \sqrt{OmegaSq}
        Omega = sqrt(OmegaSq);

        // Avoid a potential divide-by-zero
        if(OmegaSq > 1.e-30){

          // Compute Q = I + sin( dt * Omega ) * OmegaTensor / Omega - (1. - cos(dt * Omega)) * omegaTensor^2 / OmegaSq
          //           = I + scaleFactor1 * OmegaTensor + scaleFactor2 * OmegaTensorSq
          scaleFactor1 = sin(dt*Omega) / Omega;
          scaleFactor2 = -(1.0 - cos(dt*Omega)) / OmegaSq;
          MatrixMultiply(false, false, 1.0, OmegaTensor, OmegaTensor, OmegaTensorSq);
          *(QMatrix)   = 1.0 + scaleFactor1 * *(OmegaTensor)   + scaleFactor2 * *(OmegaTensorSq)   ;
          *(QMatrix+1) =       scaleFactor1 * *(OmegaTensor+1) + scaleFactor2 * *(OmegaTensorSq+1) ;
          *(QMatrix+2) =       scaleFactor1 * *(OmegaTensor+2) + scaleFactor2 * *(OmegaTensorSq+2) ;
          *(QMatrix+3) =       scaleFactor1 * *(OmegaTensor+3) + scaleFactor2 * *(OmegaTensorSq+3) ;
          *(QMatrix+4) = 1.0 + scaleFactor1 * *(OmegaTensor+4) + scaleFactor2 * *(OmegaTensorSq+4) ;
          *(QMatrix+5) =       scaleFactor1 * *(OmegaTensor+5) + scaleFactor2 * *(OmegaTensorSq+5) ;
          *(QMatrix+6) =       scaleFactor1 * *(OmegaTensor+6) + scaleFactor2 * *(OmegaTensorSq+6) ;
          *(QMatrix+7) =       scaleFactor1 * *(OmegaTensor+7) + scaleFactor2 * *(OmegaTensorSq+7) ;
          *(QMatrix+8) = 1.0 + scaleFactor1 * *(OmegaTensor+8) + scaleFactor2 * *(OmegaTensorSq+8) ;

        } else {
          *(QMatrix)   = 1.0 ; *(QMatrix+1) = 0.0 ; *(QMatrix+2) = 0.0 ;
          *(QMatrix+3) = 0.0 ; *(QMatrix+4) = 1.0 ; *(QMatrix+5) = 0.0 ;
          *(QMatrix+6) = 0.0 ; *(QMatrix+7) = 0.0 ; *(QMatrix+8) = 1.0 ;
        };

        // Compute R_STEP_NP1 = QMatrix * R_STEP_N (T&F Eq. 36)
        MatrixMultiply(false, false, 1.0, QMatrix, rotTensorN, rotTensorNP1);

        // Compute rate of stretch, Vdot = L*V - V*Omega
        // First tempA = L*V, 
        MatrixMultiply(false, false, 1.0, eulerianVelGrad, leftStretchN, tempA);

        // tempB = V*Omega
        MatrixMultiply(false, false, 1.0, leftStretchN, OmegaTensor, tempB);

        //Vdot = tempA - tempB
        for(int i=0 ; i<9 ; ++i)
          *(rateOfStretch+i) = *(tempA+i) - *(tempB+i);

        //V_STEP_NP1 = V_STEP_N + dt*Vdot
        for(int i=0 ; i<9 ; ++i)
          *(leftStretchNP1+i) = *(leftStretchN+i) + dt * *(rateOfStretch+i);

        // Compute the unrotated rate-of-deformation, d, i.e., temp = D * R
        MatrixMultiply(false, false, 1.0, rateOfDef, rotTensorNP1, temp);

        // d = Rt * temp
        MatrixMultiply(true, false, 1.0, rotTensorNP1, temp, unrotRateOfDef);

        // Store back in element-wise format
        *leftStretchXXNP1 = *(leftStretchNP1+0); *leftStretchXYNP1 = *(leftStretchNP1+1); *leftStretchXZNP1 = *(leftStretchNP1+2);
        *leftStretchYXNP1 = *(leftStretchNP1+3); *leftStretchYYNP1 = *(leftStretchNP1+4); *leftStretchYZNP1 = *(leftStretchNP1+5);
        *leftStretchZXNP1 = *(leftStretchNP1+6); *leftStretchZYNP1 = *(leftStretchNP1+7); *leftStretchZZNP1 = *(leftStretchNP1+8);
        *rotTensorXXNP1 = *(rotTensorNP1+0); *rotTensorXYNP1 = *(rotTensorNP1+1); *rotTensorXZNP1 = *(rotTensorNP1+2);
        *rotTensorYXNP1 = *(rotTensorNP1+3); *rotTensorYYNP1 = *(rotTensorNP1+4); *rotTensorYZNP1 = *(rotTensorNP1+5);
        *rotTensorZXNP1 = *(rotTensorNP1+6); *rotTensorZYNP1 = *(rotTensorNP1+7); *rotTensorZZNP1 = *(rotTensorNP1+8);
        *unrotRateOfDefXX = *(unrotRateOfDef+0); *unrotRateOfDefXY = *(unrotRateOfDef+1); *unrotRateOfDefXZ = *(unrotRateOfDef+2);
        *unrotRateOfDefYX = *(unrotRateOfDef+3); *unrotRateOfDefYY = *(unrotRateOfDef+4); *unrotRateOfDefYZ = *(unrotRateOfDef+5);
        *unrotRateOfDefZX = *(unrotRateOfDef+6); *unrotRateOfDefZY = *(unrotRateOfDef+7); *unrotRateOfDefZZ = *(unrotRateOfDef+8);
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      velGradXX += numNeighbors; velGradXY += numNeighbors; velGradXZ += numNeighbors; 
      velGradYX += numNeighbors; velGradYY += numNeighbors; velGradYZ += numNeighbors; 
      velGradZX += numNeighbors; velGradZY += numNeighbors; velGradZZ += numNeighbors;
      leftStretchXXN += numNeighbors; leftStretchXYN += numNeighbors; leftStretchXZN += numNeighbors; 
      leftStretchYXN += numNeighbors; leftStretchYYN += numNeighbors; leftStretchYZN += numNeighbors; 
      leftStretchZXN += numNeighbors; leftStretchZYN += numNeighbors; leftStretchZZN += numNeighbors;
      rotTensorXXN += numNeighbors; rotTensorXYN += numNeighbors; rotTensorXZN += numNeighbors;
      rotTensorYXN += numNeighbors; rotTensorYYN += numNeighbors; rotTensorYZN += numNeighbors;
      rotTensorZXN += numNeighbors; rotTensorZYN += numNeighbors; rotTensorZZN += numNeighbors;
      leftStretchXXNP1 += numNeighbors; leftStretchXYNP1 += numNeighbors; leftStretchXZNP1 += numNeighbors; 
      leftStretchYXNP1 += numNeighbors; leftStretchYYNP1 += numNeighbors; leftStretchYZNP1 += numNeighbors; 
      leftStretchZXNP1 += numNeighbors; leftStretchZYNP1 += numNeighbors; leftStretchZZNP1 += numNeighbors;
      rotTensorXXNP1 += numNeighbors; rotTensorXYNP1 += numNeighbors; rotTensorXZNP1 += numNeighbors;
      rotTensorYXNP1 += numNeighbors; rotTensorYYNP1 += numNeighbors; rotTensorYZNP1 += numNeighbors;
      rotTensorZXNP1 += numNeighbors; rotTensorZYNP1 += numNeighbors; rotTensorZZNP1 += numNeighbors;
      unrotRateOfDefXX += numNeighbors; unrotRateOfDefXY += numNeighbors; unrotRateOfDefXZ += numNeighbors;
      unrotRateOfDefYX += numNeighbors; unrotRateOfDefYY += numNeighbors; unrotRateOfDefYZ += numNeighbors;
      unrotRateOfDefZX += numNeighbors; unrotRateOfDefZY += numNeighbors; unrotRateOfDefZZ += numNeighbors;
    }
  }

  return returnCode;
}

template<typename ScalarT>
void rotateCauchyStress
(
    const ScalarT* rotationTensor,
    const ScalarT* unrotatedCauchyStress,
    ScalarT* rotatedCauchyStress,
    const double* flyingPointFlag,
    int numPoints
)
{
  const ScalarT* rotTensor = rotationTensor;
  const ScalarT* unrotatedStress = unrotatedCauchyStress;
  ScalarT* rotatedStress = rotatedCauchyStress;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT temp[9];

  for(int iID=0 ; iID<numPoints ; ++iID, 
        rotTensor+=9, unrotatedStress+=9, rotatedStress+=9, flyingPointFlg++){ 

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // temp = \sigma_unrot * Rt
      CORRESPONDENCE::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
      // \sigma_rot = R * temp
      CORRESPONDENCE::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);
    }
  }
}


template<typename ScalarT>
void rotateBondLevelCauchyStress(
    const ScalarT* bondLevelRotationTensorXX,
    const ScalarT* bondLevelRotationTensorXY,
    const ScalarT* bondLevelRotationTensorXZ,
    const ScalarT* bondLevelRotationTensorYX,
    const ScalarT* bondLevelRotationTensorYY,
    const ScalarT* bondLevelRotationTensorYZ,
    const ScalarT* bondLevelRotationTensorZX,
    const ScalarT* bondLevelRotationTensorZY,
    const ScalarT* bondLevelRotationTensorZZ,
    const ScalarT* bondLevelUnrotatedCauchyStressXX,
    const ScalarT* bondLevelUnrotatedCauchyStressXY,
    const ScalarT* bondLevelUnrotatedCauchyStressXZ,
    const ScalarT* bondLevelUnrotatedCauchyStressYX,
    const ScalarT* bondLevelUnrotatedCauchyStressYY,
    const ScalarT* bondLevelUnrotatedCauchyStressYZ,
    const ScalarT* bondLevelUnrotatedCauchyStressZX,
    const ScalarT* bondLevelUnrotatedCauchyStressZY,
    const ScalarT* bondLevelUnrotatedCauchyStressZZ,
    ScalarT* bondLevelRotatedCauchyStressXX,
    ScalarT* bondLevelRotatedCauchyStressXY,
    ScalarT* bondLevelRotatedCauchyStressXZ,
    ScalarT* bondLevelRotatedCauchyStressYX,
    ScalarT* bondLevelRotatedCauchyStressYY,
    ScalarT* bondLevelRotatedCauchyStressYZ,
    ScalarT* bondLevelRotatedCauchyStressZX,
    ScalarT* bondLevelRotatedCauchyStressZY,
    ScalarT* bondLevelRotatedCauchyStressZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
)
{
  const ScalarT* rotTensorXX = bondLevelRotationTensorXX;
  const ScalarT* rotTensorXY = bondLevelRotationTensorXY;
  const ScalarT* rotTensorXZ = bondLevelRotationTensorXZ;
  const ScalarT* rotTensorYX = bondLevelRotationTensorYX;
  const ScalarT* rotTensorYY = bondLevelRotationTensorYY;
  const ScalarT* rotTensorYZ = bondLevelRotationTensorYZ;
  const ScalarT* rotTensorZX = bondLevelRotationTensorZX;
  const ScalarT* rotTensorZY = bondLevelRotationTensorZY;
  const ScalarT* rotTensorZZ = bondLevelRotationTensorZZ;
  const ScalarT* unrotatedStressXX = bondLevelUnrotatedCauchyStressXX;
  const ScalarT* unrotatedStressXY = bondLevelUnrotatedCauchyStressXY;
  const ScalarT* unrotatedStressXZ = bondLevelUnrotatedCauchyStressXZ;
  const ScalarT* unrotatedStressYX = bondLevelUnrotatedCauchyStressYX;
  const ScalarT* unrotatedStressYY = bondLevelUnrotatedCauchyStressYY;
  const ScalarT* unrotatedStressYZ = bondLevelUnrotatedCauchyStressYZ;
  const ScalarT* unrotatedStressZX = bondLevelUnrotatedCauchyStressZX;
  const ScalarT* unrotatedStressZY = bondLevelUnrotatedCauchyStressZY;
  const ScalarT* unrotatedStressZZ = bondLevelUnrotatedCauchyStressZZ;
  ScalarT* rotatedStressXX = bondLevelRotatedCauchyStressXX;
  ScalarT* rotatedStressXY = bondLevelRotatedCauchyStressXY;
  ScalarT* rotatedStressXZ = bondLevelRotatedCauchyStressXZ;
  ScalarT* rotatedStressYX = bondLevelRotatedCauchyStressYX;
  ScalarT* rotatedStressYY = bondLevelRotatedCauchyStressYY;
  ScalarT* rotatedStressYZ = bondLevelRotatedCauchyStressYZ;
  ScalarT* rotatedStressZX = bondLevelRotatedCauchyStressZX;
  ScalarT* rotatedStressZY = bondLevelRotatedCauchyStressZY;
  ScalarT* rotatedStressZZ = bondLevelRotatedCauchyStressZZ;

  const double* flyingPointFlg = flyingPointFlag;

  ScalarT unrotatedStress[9], rotTensor[9], rotatedStress[9], temp[9];

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // All is bond level.
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++,
          rotTensorXX++, rotTensorXY++, rotTensorXZ++, 
          rotTensorYX++, rotTensorYY++, rotTensorYZ++, 
          rotTensorZX++, rotTensorZY++, rotTensorZZ++, 
          unrotatedStressXX++, unrotatedStressXY++, unrotatedStressXZ++, 
          unrotatedStressYX++, unrotatedStressYY++, unrotatedStressYZ++, 
          unrotatedStressZX++, unrotatedStressZY++, unrotatedStressZZ++, 
          rotatedStressXX++, rotatedStressXY++, rotatedStressXZ++, 
          rotatedStressYX++, rotatedStressYY++, rotatedStressYZ++, 
          rotatedStressZX++, rotatedStressZY++, rotatedStressZZ++){

        // write in matrix form 
        rotTensor[0] = *rotTensorXX; rotTensor[1] = *rotTensorXY; rotTensor[2] = *rotTensorXZ;
        rotTensor[3] = *rotTensorYX; rotTensor[4] = *rotTensorYY; rotTensor[5] = *rotTensorYZ;
        rotTensor[6] = *rotTensorZX; rotTensor[7] = *rotTensorZY; rotTensor[8] = *rotTensorZZ;
        unrotatedStress[0] = *unrotatedStressXX; unrotatedStress[1] = *unrotatedStressXY; unrotatedStress[2] = *unrotatedStressXZ;
        unrotatedStress[3] = *unrotatedStressYX; unrotatedStress[4] = *unrotatedStressYY; unrotatedStress[5] = *unrotatedStressYZ;
        unrotatedStress[6] = *unrotatedStressZX; unrotatedStress[7] = *unrotatedStressZY; unrotatedStress[8] = *unrotatedStressZZ;

        // temp = \sigma_unrot * Rt
        CORRESPONDENCE::MatrixMultiply(false, true, 1.0, unrotatedStress, rotTensor, temp);
        // \sigma_rot = R * temp
        CORRESPONDENCE::MatrixMultiply(false, false, 1.0, rotTensor, temp, rotatedStress);

        // update bond-level data field
        *rotatedStressXX = rotatedStress[0]; *rotatedStressXY = rotatedStress[1]; *rotatedStressXZ = rotatedStress[2]; 
        *rotatedStressYX = rotatedStress[3]; *rotatedStressYY = rotatedStress[4]; *rotatedStressYZ = rotatedStress[5]; 
        *rotatedStressZX = rotatedStress[6]; *rotatedStressZY = rotatedStress[7]; *rotatedStressZZ = rotatedStress[8]; 
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      rotTensorXX += numNeighbors; rotTensorXY += numNeighbors; rotTensorXZ += numNeighbors; 
      rotTensorYX += numNeighbors; rotTensorYY += numNeighbors; rotTensorYZ += numNeighbors; 
      rotTensorZX += numNeighbors; rotTensorZY += numNeighbors; rotTensorZZ += numNeighbors; 
      unrotatedStressXX += numNeighbors; unrotatedStressXY += numNeighbors; unrotatedStressXZ += numNeighbors; 
      unrotatedStressYX += numNeighbors; unrotatedStressYY += numNeighbors; unrotatedStressYZ += numNeighbors; 
      unrotatedStressZX += numNeighbors; unrotatedStressZY += numNeighbors; unrotatedStressZZ += numNeighbors; 
      rotatedStressXX += numNeighbors; rotatedStressXY += numNeighbors; rotatedStressXZ += numNeighbors; 
      rotatedStressYX += numNeighbors; rotatedStressYY += numNeighbors; rotatedStressYZ += numNeighbors; 
      rotatedStressZX += numNeighbors; rotatedStressZY += numNeighbors; rotatedStressZZ += numNeighbors;
    }
  }
}

template<typename ScalarT>
void computeNonhomogeneityIntegral
(
    const double* volume,
    const double* weightedVolume,
    const ScalarT* jacobianDeterminant,
    const double* horizon,
    const ScalarT* coordinates,
    const ScalarT* bondLevelCauchyStressXX,
    const ScalarT* bondLevelCauchyStressXY,
    const ScalarT* bondLevelCauchyStressXZ,
    const ScalarT* bondLevelCauchyStressYX,
    const ScalarT* bondLevelCauchyStressYY,
    const ScalarT* bondLevelCauchyStressYZ,
    const ScalarT* bondLevelCauchyStressZX,
    const ScalarT* bondLevelCauchyStressZY,
    const ScalarT* bondLevelCauchyStressZZ,
    ScalarT* nonhomogeneousIntegral,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int  numPoints
)
{
  const double* delta = horizon;
  const ScalarT* coord = coordinates;
  const ScalarT* neighborCoord;
  const double* w0 = weightedVolume;
  const double* neighborW0;
  const ScalarT* stressXX = bondLevelCauchyStressXX;
  const ScalarT* stressXY = bondLevelCauchyStressXY;
  const ScalarT* stressXZ = bondLevelCauchyStressXZ;
  const ScalarT* stressYX = bondLevelCauchyStressYX;
  const ScalarT* stressYY = bondLevelCauchyStressYY;
  const ScalarT* stressYZ = bondLevelCauchyStressYZ;
  const ScalarT* stressZX = bondLevelCauchyStressZX;
  const ScalarT* stressZY = bondLevelCauchyStressZY;
  const ScalarT* stressZZ = bondLevelCauchyStressZZ;
  ScalarT* integral = nonhomogeneousIntegral;

  const double* flyingPointFlg = flyingPointFlag;
  const double *bondDamagePtr = bondDamage;

  ScalarT deformedBondX, deformedBondY, deformedBondZ, deformedBondLength, deformedBondLengthSq;
  double neighborVolume, omega, scalarTemp;

  std::vector<ScalarT> tempVector(9);
  ScalarT* temp = &tempVector[0];

  std::vector<ScalarT> stressVector(9);
  ScalarT* stress = &stressVector[0];

  std::vector<ScalarT> integrandVector(9);
  ScalarT* integrand = &integrandVector[0];

  int neighborIndex, numNeighbors;
  const int *neighborListPtr = neighborhoodList; 
  for(int iID=0 ; iID<numPoints ; ++iID, delta++, coord+=3, w0++, integral+=9, flyingPointFlg++){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // Zero out data
      *(integral)   = 0.0 ; *(integral+1) = 0.0 ; *(integral+2) = 0.0 ;
      *(integral+3) = 0.0 ; *(integral+4) = 0.0 ; *(integral+5) = 0.0 ;
      *(integral+6) = 0.0 ; *(integral+7) = 0.0 ; *(integral+8) = 0.0 ;

      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, bondDamagePtr++,
          stressXX++, stressXY++, stressXZ++, 
          stressYX++, stressYY++, stressYZ++, 
          stressZX++, stressZY++, stressZZ++){

        neighborIndex = *neighborListPtr;
        neighborVolume = jacobianDeterminant[neighborIndex] * volume[neighborIndex];
        neighborCoord = coordinates + 3*neighborIndex;
        neighborW0 = weightedVolume + neighborIndex;

        deformedBondX = *(neighborCoord)   - *(coord);
        deformedBondY = *(neighborCoord+1) - *(coord+1);
        deformedBondZ = *(neighborCoord+2) - *(coord+2);
        deformedBondLength = sqrt(deformedBondX*deformedBondX +
                                  deformedBondY*deformedBondY +
                                  deformedBondZ*deformedBondZ);
        deformedBondLengthSq = deformedBondX*deformedBondX +
                               deformedBondY*deformedBondY +
                               deformedBondZ*deformedBondZ;

        // write the stress in matrix form 
        stress[0] = *stressXX; stress[1] = *stressXY; stress[2] = *stressXZ; 
        stress[3] = *stressYX; stress[4] = *stressYY; stress[5] = *stressYZ; 
        stress[6] = *stressZX; stress[7] = *stressZY; stress[8] = *stressZZ; 

        // delta_jp - (y_j y_p)/|y|^2
        *(temp+0) = 1.0 - deformedBondX * deformedBondX / deformedBondLengthSq;
        *(temp+1) = - deformedBondX * deformedBondY / deformedBondLengthSq;
        *(temp+2) = - deformedBondX * deformedBondZ / deformedBondLengthSq;
        *(temp+3) = *(temp+1);
        *(temp+4) = 1.0 - deformedBondY * deformedBondY / deformedBondLengthSq;
        *(temp+5) = - deformedBondY * deformedBondZ / deformedBondLengthSq;
        *(temp+6) = *(temp+2);
        *(temp+7) = *(temp+5);
        *(temp+8) = 1.0 - deformedBondZ * deformedBondZ / deformedBondLengthSq;

        // Matrix multiply the stress and the second term to compute the integrand
        MatrixMultiply(false, false, 1.0, stress, temp, integrand);

        omega = (1.0 - *bondDamagePtr) * MATERIAL_EVALUATION::scalarInfluenceFunction(deformedBondLength, *delta);

        if(omega > 0.0){
          scalarTemp = omega * (0.5 / *w0 + 0.5 / *neighborW0) * neighborVolume;

          for(int i=0; i<9; i++)
            *(integral+i) += scalarTemp * *(integrand+i);
        }
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;

      bondDamagePtr += numNeighbors; 
      stressXX += numNeighbors; stressXY += numNeighbors; stressXZ += numNeighbors; 
      stressYX += numNeighbors; stressYY += numNeighbors; stressYZ += numNeighbors; 
      stressZX += numNeighbors; stressZY += numNeighbors; stressZZ += numNeighbors;
    }
  }
}


/** Explicit template instantiation for double. */

template void TransposeMatrix<double>
(
    const double* matrix,
    double* transpose
);

template void MatrixMultiply<double>
(
    bool transA,
    bool transB,
    double alpha,
    const double* a,
    const double* b,
    double* result
);

template void rotateCauchyStress<double>
(
    const double* rotationTensor,
    const double* unrotatedCauchyStress,
    double* rotatedCauchyStress,
    int numPoints
 );

template int Invert3by3Matrix<double>
(
    const double* matrix,
    double& determinant,
     double* inverse
);

template int computeShapeTensorInverseAndApproximateDeformationGradient<double>
(
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const double* coordinates,
    double* shapeTensorInverse,
    double* deformationGradient,
    const int* neighborhoodList,
    int numPoints
);

template int computeUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* volume,
    const double* horizon,
    const double* modelCoordinates,
    const double* velocities,
    const double* deformationGradient,
    const double* shapeTensorInverse,
    const double* leftStretchTensorN,
    const double* rotationTensorN,
    double* leftStretchTensorNP1,
    double* rotationTensorNP1,
    double* unrotatedRateOfDeformation,
    const int* neighborhoodList,
    int numPoints,
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
    const double* horizon,
    const double* modelCoordinates,
    const double* coordinates,
    const double* deformationGradient,
    double* hourglassForceDensity,
    const int* neighborhoodList,
    int numPoints,
    double bulkModulus,
    double hourglassCoefficient
);

template void setOnesOnDiagonalFullTensor<double>
(
    double* tensor,
    int numPoints
);

template void computeUndamagedWeightedVolume<double>
(
    const double* volume,
    double* weightedVolume,
    const double* jacobianDeterminant,
    const double* horizon,
    const double* coordinates,
    const int* neighborhoodList,
    int numPoints
);

template void computeWeightedVolume<double>
(
    const double* volume,
    double* weightedVolume,
    const double* jacobianDeterminant,
    const double* horizon,
    const double* coordinates,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
);

template int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient<double>
(
    const double* volume,
    const double* jacobianDeterminantN,
    double* jacobianDeterminantNP1,
    const double* horizon,
    const double* coordinates,
    const double* velocities,
    double* shapeTensorInverse,
    double* velocityGradient,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints,
    double dt
);

template int computeShapeTensorInverseAndApproximateNodeLevelVelocityGradient<double>
(
    const double* volume,
    const double* jacobianDeterminantN,
    double* jacobianDeterminantNP1,
    const double* horizon,
    const double* coordinates,
    const double* velocities,
    double* shapeTensorInverse,
    double* velocityGradient,
    double* velocityGradientX,
    double* velocityGradientY,
    double* velocityGradientZ,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints,
    double dt
);

template void computeBondLevelVelocityGradient<double>
(
    const double* coordinates,
    const double* velocities,
    const double* velocityGradient,
    double* bondLevelVelocityGradientXX,
    double* bondLevelVelocityGradientXY,
    double* bondLevelVelocityGradientXZ,
    double* bondLevelVelocityGradientYX,
    double* bondLevelVelocityGradientYY,
    double* bondLevelVelocityGradientYZ,
    double* bondLevelVelocityGradientZX,
    double* bondLevelVelocityGradientZY,
    double* bondLevelVelocityGradientZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
);

template void computeBondLevelVelocityGradient<double>
(
    const double* coordinates,
    const double* velocities,
    const double* velocityGradientX,
    const double* velocityGradientY,
    const double* velocityGradientZ,
    double* bondLevelVelocityGradientXX,
    double* bondLevelVelocityGradientXY,
    double* bondLevelVelocityGradientXZ,
    double* bondLevelVelocityGradientYX,
    double* bondLevelVelocityGradientYY,
    double* bondLevelVelocityGradientYZ,
    double* bondLevelVelocityGradientZX,
    double* bondLevelVelocityGradientZY,
    double* bondLevelVelocityGradientZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
);

template int computeNodeLevelUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* velocityGradient,
    const double* leftStretchTensorN,
    const double* rotationTensorN,
    double* leftStretchTensorNP1,
    double* rotationTensorNP1,
    double* unrotatedRateOfDeformation,
    const double* flyingPointFlag,
    int numPoints,
    double dt
);

template void updateDeformationGradient<double>
(
    const double* velocityGradient,
    const double* deformationGradientN,
    double* deformationGradientNP1,
    const double* flyingPointFlag,
    int numPoints,
    double dt
);

template void computeGreenLagrangeStrain<double>
(
    const double* deformationGradient,
    double* greenLagrangeStrain,
    const double* flyingPointFlag,
    int numPoints
);

template int computeBondLevelUnrotatedRateOfDeformationAndRotationTensor<double>
(
    const double* bondLevelVelocityGradientXX, 
    const double* bondLevelVelocityGradientXY, 
    const double* bondLevelVelocityGradientXZ,
    const double* bondLevelVelocityGradientYX, 
    const double* bondLevelVelocityGradientYY, 
    const double* bondLevelVelocityGradientYZ, 
    const double* bondLevelVelocityGradientZX,
    const double* bondLevelVelocityGradientZY,
    const double* bondLevelVelocityGradientZZ,
    const double* bondLevelLeftStretchTensorXXN,
    const double* bondLevelLeftStretchTensorXYN,
    const double* bondLevelLeftStretchTensorXZN,
    const double* bondLevelLeftStretchTensorYXN,
    const double* bondLevelLeftStretchTensorYYN,
    const double* bondLevelLeftStretchTensorYZN,
    const double* bondLevelLeftStretchTensorZXN,
    const double* bondLevelLeftStretchTensorZYN,
    const double* bondLevelLeftStretchTensorZZN,
    const double* bondLevelRotationTensorXXN, 
    const double* bondLevelRotationTensorXYN, 
    const double* bondLevelRotationTensorXZN, 
    const double* bondLevelRotationTensorYXN, 
    const double* bondLevelRotationTensorYYN, 
    const double* bondLevelRotationTensorYZN, 
    const double* bondLevelRotationTensorZXN, 
    const double* bondLevelRotationTensorZYN, 
    const double* bondLevelRotationTensorZZN, 
    double* bondLevelLeftStretchTensorXXNP1,
    double* bondLevelLeftStretchTensorXYNP1,
    double* bondLevelLeftStretchTensorXZNP1,
    double* bondLevelLeftStretchTensorYXNP1,
    double* bondLevelLeftStretchTensorYYNP1,
    double* bondLevelLeftStretchTensorYZNP1,
    double* bondLevelLeftStretchTensorZXNP1,
    double* bondLevelLeftStretchTensorZYNP1,
    double* bondLevelLeftStretchTensorZZNP1,
    double* bondLevelRotationTensorXXNP1,
    double* bondLevelRotationTensorXYNP1,
    double* bondLevelRotationTensorXZNP1,
    double* bondLevelRotationTensorYXNP1,
    double* bondLevelRotationTensorYYNP1,
    double* bondLevelRotationTensorYZNP1,
    double* bondLevelRotationTensorZXNP1,
    double* bondLevelRotationTensorZYNP1,
    double* bondLevelRotationTensorZZNP1,
    double* bondLevelUnrotatedRateOfDeformationXX,
    double* bondLevelUnrotatedRateOfDeformationXY,
    double* bondLevelUnrotatedRateOfDeformationXZ,
    double* bondLevelUnrotatedRateOfDeformationYX,
    double* bondLevelUnrotatedRateOfDeformationYY,
    double* bondLevelUnrotatedRateOfDeformationYZ,
    double* bondLevelUnrotatedRateOfDeformationZX,
    double* bondLevelUnrotatedRateOfDeformationZY,
    double* bondLevelUnrotatedRateOfDeformationZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints,
    double dt
);

template void rotateCauchyStress<double>
(
    const double* rotationTensor,
    const double* unrotatedCauchyStress,
    double* rotatedCauchyStress,
    const double* flyingPointFlag,
    int numPoints
);

template void rotateBondLevelCauchyStress(
    const double* bondLevelRotationTensorXX,
    const double* bondLevelRotationTensorXY,
    const double* bondLevelRotationTensorXZ,
    const double* bondLevelRotationTensorYX,
    const double* bondLevelRotationTensorYY,
    const double* bondLevelRotationTensorYZ,
    const double* bondLevelRotationTensorZX,
    const double* bondLevelRotationTensorZY,
    const double* bondLevelRotationTensorZZ,
    const double* bondLevelUnrotatedCauchyStressXX,
    const double* bondLevelUnrotatedCauchyStressXY,
    const double* bondLevelUnrotatedCauchyStressXZ,
    const double* bondLevelUnrotatedCauchyStressYX,
    const double* bondLevelUnrotatedCauchyStressYY,
    const double* bondLevelUnrotatedCauchyStressYZ,
    const double* bondLevelUnrotatedCauchyStressZX,
    const double* bondLevelUnrotatedCauchyStressZY,
    const double* bondLevelUnrotatedCauchyStressZZ,
    double* bondLevelRotatedCauchyStressXX,
    double* bondLevelRotatedCauchyStressXY,
    double* bondLevelRotatedCauchyStressXZ,
    double* bondLevelRotatedCauchyStressYX,
    double* bondLevelRotatedCauchyStressYY,
    double* bondLevelRotatedCauchyStressYZ,
    double* bondLevelRotatedCauchyStressZX,
    double* bondLevelRotatedCauchyStressZY,
    double* bondLevelRotatedCauchyStressZZ,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    int numPoints
);

template void computeNonhomogeneityIntegral<double>
(
    const double* volume,
    const double* weightedVolume,
    const double* jacobianDeterminant,
    const double* horizon,
    const double* coordinates,
    const double* bondLevelCauchyStressXX,
    const double* bondLevelCauchyStressXY,
    const double* bondLevelCauchyStressXZ,
    const double* bondLevelCauchyStressYX,
    const double* bondLevelCauchyStressYY,
    const double* bondLevelCauchyStressYZ,
    const double* bondLevelCauchyStressZX,
    const double* bondLevelCauchyStressZY,
    const double* bondLevelCauchyStressZZ,
    double* nonhomogeneousIntegral,
    const double* flyingPointFlag,
    const double* bondDamage,
    const int* neighborhoodList,
    int numPoints
);


/** Explicit template instantiation for Sacado::Fad::DFad<double>. */

template void TransposeMatrix<Sacado::Fad::DFad<double> >
(
    const Sacado::Fad::DFad<double>* matrix,
    Sacado::Fad::DFad<double>* transpose
);

template void MatrixMultiply<Sacado::Fad::DFad<double> >
(
    bool transA,
    bool transB,
    Sacado::Fad::DFad<double> alpha,
    const Sacado::Fad::DFad<double>* a,
    const Sacado::Fad::DFad<double>* b,
    Sacado::Fad::DFad<double>* result
);

template int Invert3by3Matrix<Sacado::Fad::DFad<double> >
(
    const Sacado::Fad::DFad<double>* matrix,
    Sacado::Fad::DFad<double>& determinant,
    Sacado::Fad::DFad<double>* inverse
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

}