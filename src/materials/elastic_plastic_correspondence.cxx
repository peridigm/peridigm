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

#include "elastic_plastic_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h>

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticPerfectlyPlasticCauchyStress
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* cauchyStressN, 
    ScalarT* cauchyStressNP1, 
    ScalarT* vonMisesStress,
    const ScalarT* equivalentPlasticStrainN,
    ScalarT* equivalentPlasticStrainNP1,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
)
{
  
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* stressN = cauchyStressN;
  ScalarT* stressNP1 = cauchyStressNP1;

  ScalarT* vmStress = vonMisesStress;

  const ScalarT* eqpsN = equivalentPlasticStrainN;
  ScalarT* eqpsNP1 = equivalentPlasticStrainNP1;

  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];
  ScalarT deviatoricStressN[9];
  ScalarT deviatoricStressMagnitudeN;

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunction;

  for(int iID=0 ; iID<numPoints ; ++iID, rateOfDef+=9, stressN+=9,
      stressNP1+=9, ++vmStress,++eqpsN,++eqpsNP1){

      //strainInc = dt * rateOfDef
      for(int i = 0; i < 9; i++){
          strainInc[i] = *(rateOfDef+i)*dt;
          deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //Compute an elastic ``trail stress''
      for(int i = 0; i < 9; i++){
          *(stressNP1+i) = *(stressN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(stressNP1) += bulkMod*dilatationInc;
      *(stressNP1+4) += bulkMod*dilatationInc;
      *(stressNP1+8) += bulkMod*dilatationInc;

      sphericalStressNP1 = (*(stressNP1) + *(stressNP1+4) + *(stressNP1+8))/3.0;

      // Compute the ``trial'' von Mises stress
      for(int i = 0; i < 9; i++){
          deviatoricStressNP1[i] = *(stressNP1+i);
      }
      deviatoricStressNP1[0] -= sphericalStressNP1;
      deviatoricStressNP1[4] -= sphericalStressNP1;
      deviatoricStressNP1[8] -= sphericalStressNP1;

      // Compute \sigma_ij * \sigma_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
          for (int i = 0; i < 3; i++) {
              tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
      }

      *vmStress = sqrt(3.0/2.0*tempScalar);

      yieldFunction = yieldStress;

      //If true, the step is plastic and we need to return to the yield
      //surface.  
      if(*vmStress > yieldFunction){

        // Avoid divide-by-zero
        deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));
        //For perfectly plastic deformations we just have a constant factor to 
        //multiply the trial deviatoric stress
        tempScalar = sqrt(2.0/3.0)*yieldFunction/deviatoricStressMagnitudeNP1;

        // Return the deviatoric stress to the yield surface
        for (int i = 0; i < 9; i++) {
          deviatoricStressNP1[i] *= tempScalar; 
            *(stressNP1+i) = deviatoricStressNP1[i];
        }

        // Update the Cauchy Stress
        *stressNP1 += sphericalStressNP1;
        *(stressNP1+4) += sphericalStressNP1;
        *(stressNP1+8) += sphericalStressNP1;

        // Update the von Mises stress now that the state of stress is on the
        // yield surface
        tempScalar = 0.0;
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
        }

        *vmStress = sqrt(3.0/2.0*tempScalar);

        /////////////////////////////////////////////////////////
        //
        // Update the equivalent plastic strain
        //
        // The algorithm below is generic and should not need to be modified for
        // any J2 plasticity yield surface.  It uses the difference in the yield
        // surface location at the NP1 and N steps to increment eqps regardless
        // of how the plastic multiplier was found in the yield surface
        // evaluation.
        //
        // First go back to step N and compute deviatoric stress and its
        // magnitude.  We didn't do this earlier because it wouldn't be necassary
        // if the step is elastic.
        sphericalStressN = (*stressN + *(stressN+4) + *(stressN+8))/3.0;

        for (int i = 0; i < 9; i++) {
          deviatoricStressN[i] = *(stressN+i);
        }
        deviatoricStressN[0] -= sphericalStressN;
        deviatoricStressN[4] -= sphericalStressN;
        deviatoricStressN[8] -= sphericalStressN;

        tempScalar = 0.0;
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
          }
        }

        //Ensure that this is at least a very small number to avoid a divide by
        //zero
        deviatoricStressMagnitudeN = std::max(1.0e-20,sqrt(tempScalar));

        for (int i = 0; i < 9; i++) {
          //tempA -- The plastic deviatoric strain increment tensor \Delta e_{plastic} = \Delta e_{total} - \Delta e_{elastic}
          tempA[i] = deviatoricStrainInc[i] - (deviatoricStressNP1[i] - deviatoricStressN[i]) / 2.0 / shearMod;
          //tempB -- Deviatoric stress increment.  This is effectively an average of the deviatoric stress 
          //direction unit tensors at the half-step between steps NP1 and N
          tempB[i] = (deviatoricStressNP1[i]/deviatoricStressMagnitudeNP1 + 
                      deviatoricStressN[i]/deviatoricStressMagnitudeN)/2.0;
        }

        // Contract the two tensors. This represents a projection of the plastic
        // strain increment tensor onto the "direction" of deviatoric stress
        // increment
        tempScalar = 0.0;
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            tempScalar += tempA[i+3*j] * tempB[i+3*j];
          }
        }

        // Increment the plastic strain
        *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);

      } else {
        // The step is elastic
        *eqpsNP1 = *eqpsN;
      };

  }
}

template<typename ScalarT>
void updateElasticPerfectlyPlasticCauchyStress
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* cauchyStressN, 
    ScalarT* cauchyStressNP1, 
    ScalarT* vonMisesStress,
    const ScalarT* equivalentPlasticStrainN,
    ScalarT* equivalentPlasticStrainNP1,
    ScalarT* stressTriaxiality,
    const double* flyingPointFlag,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
)
{
  
  const ScalarT* rateOfDef = unrotatedRateOfDeformation;
  const ScalarT* stressN = cauchyStressN;
  ScalarT* stressNP1 = cauchyStressNP1;

  ScalarT* vmStress = vonMisesStress;

  const ScalarT* eqpsN = equivalentPlasticStrainN;
  ScalarT* eqpsNP1 = equivalentPlasticStrainNP1;

  ScalarT* triaxiality = stressTriaxiality;
  const double* flyingPointFlg = flyingPointFlag;

  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];
  ScalarT deviatoricStressN[9];
  ScalarT deviatoricStressMagnitudeN;

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunction;

  for(int iID=0 ; iID<numPoints ; ++iID, rateOfDef+=9, stressN+=9, stressNP1+=9, 
      ++vmStress, ++eqpsN, ++eqpsNP1, ++triaxiality, ++flyingPointFlg){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      //strainInc = dt * rateOfDef
      for(int i = 0; i < 9; i++){
        strainInc[i] = *(rateOfDef+i)*dt;
        deviatoricStrainInc[i] = strainInc[i];
      }

      //dilatation
      dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

      //deviatoric strain
      deviatoricStrainInc[0] -= dilatationInc/3.0;
      deviatoricStrainInc[4] -= dilatationInc/3.0;
      deviatoricStrainInc[8] -= dilatationInc/3.0;

      //Compute an elastic ``trail stress''
      for(int i = 0; i < 9; i++){
        *(stressNP1+i) = *(stressN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(stressNP1) += bulkMod*dilatationInc;
      *(stressNP1+4) += bulkMod*dilatationInc;
      *(stressNP1+8) += bulkMod*dilatationInc;

      sphericalStressNP1 = (*(stressNP1) + *(stressNP1+4) + *(stressNP1+8))/3.0;

      // Compute the ``trial'' von Mises stress
      for(int i = 0; i < 9; i++){
        deviatoricStressNP1[i] = *(stressNP1+i);
      }
      deviatoricStressNP1[0] -= sphericalStressNP1;
      deviatoricStressNP1[4] -= sphericalStressNP1;
      deviatoricStressNP1[8] -= sphericalStressNP1;

      // Compute \sigma_ij * \sigma_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
        for(int i = 0; i < 3; i++){
          tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
        }
      }

      *vmStress = sqrt(3.0/2.0*tempScalar);

      yieldFunction = yieldStress;

      //If true, the step is plastic and we need to return to the yield
      //surface.  
      if(*vmStress > yieldFunction) {

        // Avoid divide-by-zero
        deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));
        //For perfectly plastic deformations we just have a constant factor to 
        //multiply the trial deviatoric stress
        tempScalar = sqrt(2.0/3.0)*yieldFunction/deviatoricStressMagnitudeNP1;

        // Return the deviatoric stress to the yield surface
        for (int i = 0; i < 9; i++) {
          deviatoricStressNP1[i] *= tempScalar; 
          *(stressNP1+i) = deviatoricStressNP1[i];
        }

        // Update the Cauchy Stress
        *stressNP1 += sphericalStressNP1;
        *(stressNP1+4) += sphericalStressNP1;
        *(stressNP1+8) += sphericalStressNP1;

        // Update the von Mises stress now that the state of stress is on the
        // yield surface
        tempScalar = 0.0;
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
        }

        *vmStress = sqrt(3.0/2.0*tempScalar);

        /////////////////////////////////////////////////////////
        //
        // Update the equivalent plastic strain          //
        // The algorithm below is generic and should not need to be modified for
        // any J2 plasticity yield surface.  It uses the difference in the yield
        // surface location at the NP1 and N steps to increment eqps regardless
        // of how the plastic multiplier was found in the yield surface
        // evaluation.
        //
        // First go back to step N and compute deviatoric stress and its
        // magnitude.  We didn't do this earlier because it wouldn't be necassary
        // if the step is elastic.
        sphericalStressN = (*stressN + *(stressN+4) + *(stressN+8))/3.0;

        for(int i = 0; i < 9; i++){
          deviatoricStressN[i] = *(stressN+i);
        }
        deviatoricStressN[0] -= sphericalStressN;
        deviatoricStressN[4] -= sphericalStressN;
        deviatoricStressN[8] -= sphericalStressN;

        tempScalar = 0.0;
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
          }
        }

        //Ensure that this is at least a very small number to avoid a divide by
        //zero
        deviatoricStressMagnitudeN = std::max(1.0e-20,sqrt(tempScalar));

        for(int i = 0; i < 9; i++){
          //tempA -- The plastic deviatoric strain increment tensor \Delta e_{plastic} = \Delta e_{total} - \Delta e_{elastic}
          tempA[i] = deviatoricStrainInc[i] - (deviatoricStressNP1[i] - deviatoricStressN[i]) / 2.0 / shearMod;
          //tempB -- Deviatoric stress increment.  This is effectively an average of the deviatoric stress 
          //direction unit tensors at the half-step between steps NP1 and N
          tempB[i] = (deviatoricStressNP1[i]/deviatoricStressMagnitudeNP1 + 
                      deviatoricStressN[i]/deviatoricStressMagnitudeN)/2.0;
        }
          
        // Contract the two tensors. This represents a projection of the plastic
        // strain increment tensor onto the "direction" of deviatoric stress
        // increment
        tempScalar = 0.0;
        for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
            tempScalar += tempA[i+3*j] * tempB[i+3*j];
          }
        }

        // Increment the plastic strain
        *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);

      } else {
        // The step is elastic
        *eqpsNP1 = *eqpsN;
      };

      // compute stress triaxiality
      *triaxiality = (*(stressNP1+0) + *(stressNP1+4) + *(stressNP1+8)) / (3.0 * std::max(1.0e-20,*vmStress));
    }
  }
}

template<typename ScalarT>
void updateBondLevelElasticPerfectlyPlasticCauchyStress
(
    const ScalarT* bondLevelUnrotatedRateOfDeformationXX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationXY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationXZ, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationYZ, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZX, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZY, 
    const ScalarT* bondLevelUnrotatedRateOfDeformationZZ, 
    const ScalarT* bondLevelUnrotatedCauchyStressXXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressXYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressXZN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressYZN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZXN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZYN, 
    const ScalarT* bondLevelUnrotatedCauchyStressZZN, 
    ScalarT* bondLevelUnrotatedCauchyStressXXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressXYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressXZNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressYZNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZXNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZYNP1, 
    ScalarT* bondLevelUnrotatedCauchyStressZZNP1, 
    ScalarT* bondLevelVonMisesStress,
    const ScalarT* bondLevelEquivalentPlasticStrainN,
    ScalarT* bondLevelEquivalentPlasticStrainNP1,
    ScalarT* bondLevelStressTriaxiality,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
)
{

  const ScalarT* rateOfDefXX = bondLevelUnrotatedRateOfDeformationXX;
  const ScalarT* rateOfDefXY = bondLevelUnrotatedRateOfDeformationXY;
  const ScalarT* rateOfDefXZ = bondLevelUnrotatedRateOfDeformationXZ;
  const ScalarT* rateOfDefYX = bondLevelUnrotatedRateOfDeformationYX;
  const ScalarT* rateOfDefYY = bondLevelUnrotatedRateOfDeformationYY;
  const ScalarT* rateOfDefYZ = bondLevelUnrotatedRateOfDeformationYZ;
  const ScalarT* rateOfDefZX = bondLevelUnrotatedRateOfDeformationZX;
  const ScalarT* rateOfDefZY = bondLevelUnrotatedRateOfDeformationZY;
  const ScalarT* rateOfDefZZ = bondLevelUnrotatedRateOfDeformationZZ;
  const ScalarT* stressXXN = bondLevelUnrotatedCauchyStressXXN;
  const ScalarT* stressXYN = bondLevelUnrotatedCauchyStressXYN;
  const ScalarT* stressXZN = bondLevelUnrotatedCauchyStressXZN;
  const ScalarT* stressYXN = bondLevelUnrotatedCauchyStressYXN;
  const ScalarT* stressYYN = bondLevelUnrotatedCauchyStressYYN;
  const ScalarT* stressYZN = bondLevelUnrotatedCauchyStressYZN;
  const ScalarT* stressZXN = bondLevelUnrotatedCauchyStressZXN;
  const ScalarT* stressZYN = bondLevelUnrotatedCauchyStressZYN;
  const ScalarT* stressZZN = bondLevelUnrotatedCauchyStressZZN;
  ScalarT* stressXXNP1 = bondLevelUnrotatedCauchyStressXXNP1;
  ScalarT* stressXYNP1 = bondLevelUnrotatedCauchyStressXYNP1;
  ScalarT* stressXZNP1 = bondLevelUnrotatedCauchyStressXZNP1;
  ScalarT* stressYXNP1 = bondLevelUnrotatedCauchyStressYXNP1;
  ScalarT* stressYYNP1 = bondLevelUnrotatedCauchyStressYYNP1;
  ScalarT* stressYZNP1 = bondLevelUnrotatedCauchyStressYZNP1;
  ScalarT* stressZXNP1 = bondLevelUnrotatedCauchyStressZXNP1;
  ScalarT* stressZYNP1 = bondLevelUnrotatedCauchyStressZYNP1;
  ScalarT* stressZZNP1 = bondLevelUnrotatedCauchyStressZZNP1;
  
  ScalarT* vmStress = bondLevelVonMisesStress;

  const ScalarT* eqpsN = bondLevelEquivalentPlasticStrainN;
  ScalarT* eqpsNP1 = bondLevelEquivalentPlasticStrainNP1;

  ScalarT* triaxiality = bondLevelStressTriaxiality;
  const double* flyingPointFlg = flyingPointFlag;

  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];
  ScalarT deviatoricStressN[9];
  ScalarT deviatoricStressMagnitudeN;

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunction;

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID, ++flyingPointFlg){

    // if the node is not flying, update the values. Otherwise, just skip
    if(*flyingPointFlg < 0.0){

      // All is bond level.
      numNeighbors = *neighborListPtr; neighborListPtr++;
      for(int n=0; n<numNeighbors; n++, neighborListPtr++, 
          rateOfDefXX++, rateOfDefXY++, rateOfDefXZ++, 
          rateOfDefYX++, rateOfDefYY++, rateOfDefYZ++, 
          rateOfDefZX++, rateOfDefZY++, rateOfDefZZ++,
          stressXXN++, stressXYN++, stressXZN++, 
          stressYXN++, stressYYN++, stressYZN++, 
          stressZXN++, stressZYN++, stressZZN++,
          stressXXNP1++, stressXYNP1++, stressXZNP1++, 
          stressYXNP1++, stressYYNP1++, stressYZNP1++, 
          stressZXNP1++, stressZYNP1++, stressZZNP1++,
          vmStress++, eqpsN++, eqpsNP1++, triaxiality++){

        //strainInc = dt * rateOfDef
        strainInc[0] = *rateOfDefXX*dt; strainInc[1] = *rateOfDefXY*dt; strainInc[2] = *rateOfDefXZ*dt;
        strainInc[3] = *rateOfDefYX*dt; strainInc[4] = *rateOfDefYY*dt; strainInc[5] = *rateOfDefYZ*dt;
        strainInc[6] = *rateOfDefZX*dt; strainInc[7] = *rateOfDefZY*dt; strainInc[8] = *rateOfDefZZ*dt;

        //dilatation
        dilatationInc = strainInc[0] + strainInc[4] + strainInc[8];

        //deviatoric strain
        for(int i = 0; i < 9; i++){
          deviatoricStrainInc[i] = strainInc[i];
        }
        deviatoricStrainInc[0] -= dilatationInc/3.0;
        deviatoricStrainInc[4] -= dilatationInc/3.0;
        deviatoricStrainInc[8] -= dilatationInc/3.0;

        //Compute an elastic ``trail stress''
        *stressXXNP1 = *stressXXN + deviatoricStrainInc[0]*2.0*shearMod;
        *stressXYNP1 = *stressXYN + deviatoricStrainInc[1]*2.0*shearMod;
        *stressXZNP1 = *stressXZN + deviatoricStrainInc[2]*2.0*shearMod;
        *stressYXNP1 = *stressYXN + deviatoricStrainInc[3]*2.0*shearMod;
        *stressYYNP1 = *stressYYN + deviatoricStrainInc[4]*2.0*shearMod;
        *stressYZNP1 = *stressYZN + deviatoricStrainInc[5]*2.0*shearMod;
        *stressZXNP1 = *stressZXN + deviatoricStrainInc[6]*2.0*shearMod;
        *stressZYNP1 = *stressZYN + deviatoricStrainInc[7]*2.0*shearMod;
        *stressZZNP1 = *stressZZN + deviatoricStrainInc[8]*2.0*shearMod;

        *stressXXNP1 += bulkMod*dilatationInc;
        *stressYYNP1 += bulkMod*dilatationInc;
        *stressZZNP1 += bulkMod*dilatationInc;

        sphericalStressNP1 = (*stressXXNP1 + *stressYYNP1 + *stressZZNP1)/3.0;

        // Compute the ``trial'' von Mises stress
        deviatoricStressNP1[0] = *stressXXNP1; deviatoricStressNP1[1] = *stressXYNP1; deviatoricStressNP1[2] = *stressXZNP1;
        deviatoricStressNP1[3] = *stressYXNP1; deviatoricStressNP1[4] = *stressYYNP1; deviatoricStressNP1[5] = *stressYZNP1;
        deviatoricStressNP1[6] = *stressZXNP1; deviatoricStressNP1[7] = *stressZYNP1; deviatoricStressNP1[8] = *stressZZNP1;

        deviatoricStressNP1[0] -= sphericalStressNP1;
        deviatoricStressNP1[4] -= sphericalStressNP1;
        deviatoricStressNP1[8] -= sphericalStressNP1;

        // Compute \stress_ij * \stress_ij
        tempScalar = 0.0;
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
        }

        *vmStress = sqrt(3.0/2.0*tempScalar);

        yieldFunction = yieldStress;

        //If true, the step is plastic and we need to return to the yield
        //surface.  
        if(*vmStress > yieldFunction){

          // Avoid divide-by-zero
          deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));
          //For perfectly plastic deformations we just have a constant factor to 
          //multiply the trial deviatoric stress
          tempScalar = sqrt(2.0/3.0)*yieldFunction/deviatoricStressMagnitudeNP1;

          // Return the deviatoric stress to the yield surface
          for(int i = 0; i < 9; i++){
            deviatoricStressNP1[i] *= tempScalar; 
          }
          *stressXXNP1 = deviatoricStressNP1[0];
          *stressXYNP1 = deviatoricStressNP1[1];
          *stressXZNP1 = deviatoricStressNP1[2];
          *stressYXNP1 = deviatoricStressNP1[3];
          *stressYYNP1 = deviatoricStressNP1[4];
          *stressYZNP1 = deviatoricStressNP1[5];
          *stressZXNP1 = deviatoricStressNP1[6];
          *stressZYNP1 = deviatoricStressNP1[7];
          *stressZZNP1 = deviatoricStressNP1[8];

          // Update the Cauchy Stress
          *stressXXNP1 += sphericalStressNP1;
          *stressYYNP1 += sphericalStressNP1;
          *stressZZNP1 += sphericalStressNP1;

          // Update the von Mises stress now that the state of stress is on the
          // yield surface
          tempScalar = 0.0;
          for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
              tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
            }
          }

          *vmStress = sqrt(3.0/2.0*tempScalar);

          /////////////////////////////////////////////////////////
          //
          // Update the equivalent plastic strain
          //
          // The algorithm below is generic and should not need to be modified for
          // any J2 plasticity yield surface.  It uses the difference in the yield
          // surface location at the NP1 and N steps to increment eqps regardless
          // of how the plastic multiplier was found in the yield surface
          // evaluation.
          //
          // First go back to step N and compute deviatoric stress and its
          // magnitude.  We didn't do this earlier because it wouldn't be necessary
          // if the step is elastic.
          sphericalStressN = (*stressXXN + *stressYYN + *stressZZN)/3.0;

          deviatoricStressN[0] = *stressXXN; deviatoricStressN[1] = *stressXYN; deviatoricStressN[2] = *stressXZN;
          deviatoricStressN[3] = *stressYXN; deviatoricStressN[4] = *stressYYN; deviatoricStressN[5] = *stressYZN;
          deviatoricStressN[6] = *stressZXN; deviatoricStressN[7] = *stressZYN; deviatoricStressN[8] = *stressZZN;

          deviatoricStressN[0] -= sphericalStressN;
          deviatoricStressN[4] -= sphericalStressN;
          deviatoricStressN[8] -= sphericalStressN;

          tempScalar = 0.0;
          for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
              tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
            }
          }

          //Ensure that this is at least a very small number to avoid a divide by
          //zero
          deviatoricStressMagnitudeN = std::max(1.0e-20,sqrt(tempScalar));

          for(int i = 0; i < 9; i++){
            //tempA -- The plastic deviatoric strain increment tensor \Delta e_{plastic} = \Delta e_{total} - \Delta e_{elastic}
            tempA[i] = deviatoricStrainInc[i] - (deviatoricStressNP1[i] - deviatoricStressN[i]) / 2.0 / shearMod;
            //tempB -- Deviatoric stress increment.  This is effectively an average of the deviatoric stress 
            //direction unit tensors at the half-step between steps NP1 and N
            tempB[i] = (deviatoricStressNP1[i]/deviatoricStressMagnitudeNP1 + 
                        deviatoricStressN[i]/deviatoricStressMagnitudeN)/2.0;
          }

          // Contract the two tensors. This represents a projection of the plastic
          // strain increment tensor onto the "direction" of deviatoric stress
          // increment
          tempScalar = 0.0;
          for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
              tempScalar += tempA[i+3*j] * tempB[i+3*j];
            }
          }

          // Increment the plastic strain
          *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);

        } else {
          // The step is elastic
          *eqpsNP1 = *eqpsN;
        };

        // compute stress triaxiality
        *triaxiality = (*stressXXNP1 + *stressYYNP1 + *stressZZNP1) / (3.0 * std::max(1.0e-20,*vmStress));
      }
    }
    else{
      // adjust the neighborhood pointer for the next point
      numNeighbors = *neighborListPtr; 
      neighborListPtr += numNeighbors+1;
      rateOfDefXX += numNeighbors; rateOfDefXY += numNeighbors; rateOfDefXZ += numNeighbors; 
      rateOfDefYX += numNeighbors; rateOfDefYY += numNeighbors; rateOfDefYZ += numNeighbors; 
      rateOfDefZX += numNeighbors; rateOfDefZY += numNeighbors; rateOfDefZZ += numNeighbors;
      stressXXN += numNeighbors; stressXYN += numNeighbors; stressXZN += numNeighbors; 
      stressYXN += numNeighbors; stressYYN += numNeighbors; stressYZN += numNeighbors; 
      stressZXN += numNeighbors; stressZYN += numNeighbors; stressZZN += numNeighbors;
      stressXXNP1 += numNeighbors; stressXYNP1 += numNeighbors; stressXZNP1 += numNeighbors; 
      stressYXNP1 += numNeighbors; stressYYNP1 += numNeighbors; stressYZNP1 += numNeighbors; 
      stressZXNP1 += numNeighbors; stressZYNP1 += numNeighbors; stressZZNP1 += numNeighbors;
      vmStress += numNeighbors; eqpsN += numNeighbors; eqpsNP1 += numNeighbors; triaxiality += numNeighbors;
    }
  }
}

// Explicit template instantiation for double
template void updateElasticPerfectlyPlasticCauchyStress<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
);

template void updateElasticPerfectlyPlasticCauchyStress<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    double* stressTriaxiality,
    const double* flyingPointFlag,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
);

template void updateBondLevelElasticPerfectlyPlasticCauchyStress<double>
(
    const double* bondLevelUnrotatedRateOfDeformationXX, 
    const double* bondLevelUnrotatedRateOfDeformationXY, 
    const double* bondLevelUnrotatedRateOfDeformationXZ, 
    const double* bondLevelUnrotatedRateOfDeformationYX, 
    const double* bondLevelUnrotatedRateOfDeformationYY, 
    const double* bondLevelUnrotatedRateOfDeformationYZ, 
    const double* bondLevelUnrotatedRateOfDeformationZX, 
    const double* bondLevelUnrotatedRateOfDeformationZY, 
    const double* bondLevelUnrotatedRateOfDeformationZZ, 
    const double* bondLevelUnrotatedCauchyStressXXN, 
    const double* bondLevelUnrotatedCauchyStressXYN, 
    const double* bondLevelUnrotatedCauchyStressXZN, 
    const double* bondLevelUnrotatedCauchyStressYXN, 
    const double* bondLevelUnrotatedCauchyStressYYN, 
    const double* bondLevelUnrotatedCauchyStressYZN, 
    const double* bondLevelUnrotatedCauchyStressZXN, 
    const double* bondLevelUnrotatedCauchyStressZYN, 
    const double* bondLevelUnrotatedCauchyStressZZN, 
    double* bondLevelUnrotatedCauchyStressXXNP1, 
    double* bondLevelUnrotatedCauchyStressXYNP1, 
    double* bondLevelUnrotatedCauchyStressXZNP1, 
    double* bondLevelUnrotatedCauchyStressYXNP1, 
    double* bondLevelUnrotatedCauchyStressYYNP1, 
    double* bondLevelUnrotatedCauchyStressYZNP1, 
    double* bondLevelUnrotatedCauchyStressZXNP1, 
    double* bondLevelUnrotatedCauchyStressZYNP1, 
    double* bondLevelUnrotatedCauchyStressZZNP1, 
    double* bondLevelVonMisesStress,
    const double* bondLevelEquivalentPlasticStrainN,
    double* bondLevelEquivalentPlasticStrainNP1,
    double* bondLevelStressTriaxiality,
    const double* flyingPointFlag,
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void updateElasticPerfectlyPlasticCauchyStress<Sacado::Fad::DFad<double> >
(
    const Sacado::Fad::DFad<double>* unrotatedRateOfDeformation, 
    const Sacado::Fad::DFad<double>* cauchyStressN, 
    Sacado::Fad::DFad<double>* cauchyStressNP1, 
    Sacado::Fad::DFad<double>* vonMisesStress,
    const Sacado::Fad::DFad<double>* equivalentPlasticStrainN,
    Sacado::Fad::DFad<double>* equivalentPlasticStrainNP1,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double dt
);

}