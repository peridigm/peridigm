//! \file isotropic_hardening_correspondence.cxx

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

#include "isotropic_hardening_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h>

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticIsotropicHardeningPlasticCauchyStress
(
    const double* modelCoordinates,
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
    const double hardMod,
    const bool isFlaw,
    const double flawLocationX,
    const double flawLocationY,
    const double flawLocationZ,
    const double flawSize,
    const double flawMagnitude,
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

  ScalarT scalarDeviatoricStrainInc;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunctionVal;


  ScalarT deltaLambda;

  double reducedYieldStress;
  const double* modelCoord = modelCoordinates;


  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, 
      rateOfDef+=9, stressN+=9, stressNP1+=9,
      ++vmStress,++eqpsN,++eqpsNP1){

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

    // Compute \S_ij * \S_ij
    tempScalar = 0.0;
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
      }
    }

    // Avoid divide-by-zero
    deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));

    *vmStress = sqrt(3.0/2.0*tempScalar);

    // Reduce yield stress if flaws are present
    if(isFlaw){
      //Increment the pointer
      reducedYieldStress = yieldStress * (1.0 - flawMagnitude 
                            * exp( (
                            (- ( *(modelCoord) - flawLocationX)) * (*(modelCoord) - flawLocationX) -
                            (( *(modelCoord+1) - flawLocationY)) * (*(modelCoord+1) - flawLocationY) -
                            (( *(modelCoord+2) - flawLocationZ)) * (*(modelCoord+2) - flawLocationZ)
                            ) / flawSize / flawSize
                            ));
    } else {
      //Without flaws the reduced yield stress is the yield stress.
      reducedYieldStress = yieldStress;
    }

    // Elastic or plastic?
    if (*vmStress < reducedYieldStress){
      // The step is definitely elastic, so skip the yield surface
      // evaluation.
      *eqpsNP1 = *eqpsN;

    } else {
      // The step could be plastic, we have to solve for the current value
      // of the yield function to find out.  This is because the yield
      // function is dependent on eqps and can change over a load step.


      // Compute \S_ij * \epsilon_inc_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
        for(int i = 0; i < 3; i++){
          tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStrainInc[i+3*j];
        }
      }

      scalarDeviatoricStrainInc = tempScalar / deviatoricStressMagnitudeNP1;

      // First go back to step N and compute deviatoric stress and its
      // magnitude.
      sphericalStressN = (*(stressN) + *(stressN+4) + *(stressN+8))/3.0;

      for(int i = 0; i < 9; i++){
        deviatoricStressN[i] = *(stressN+i);
      }
      deviatoricStressN[0] -= sphericalStressN;
      deviatoricStressN[4] -= sphericalStressN;
      deviatoricStressN[8] -= sphericalStressN;

      // Compute \S_ij * \S_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
        for(int i = 0; i < 3; i++){
          tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
        }
      }

      deviatoricStressMagnitudeN = sqrt(tempScalar);

      //Solve for deltaLambda
      deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod + 3.0 * 
                    deviatoricStressMagnitudeN  - sqrt(6.0) * reducedYieldStress) /
                    (2.0 * (hardMod + 3.0 * shearMod));


      //Increment the plastic strain for the purposes of evaluating the
      //yield surface
      *eqpsNP1 = *eqpsN + sqrt(2.0 / 3.0) * deltaLambda;
      //Evaluate the extent of the yield surface
      yieldFunctionVal = reducedYieldStress + hardMod * (*eqpsNP1);

      //If true, the step is plastic and we need to return to the yield
      //surface.  
      if(*vmStress > yieldFunctionVal){

        //multiply the trial deviatoric stress
        tempScalar = sqrt(2.0/3.0)*yieldFunctionVal/deviatoricStressMagnitudeNP1;

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
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += tempA[i+3*j] * tempB[i+3*j];
          }
        }

        // Increment the plastic strain
        *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);

        } else {
          // The step is elastic
          *eqpsNP1 = *eqpsN;
        }
      }
  }
}

template<typename ScalarT>
void updateElasticIsotropicPowerlawHardeningPlasticCauchyStress
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
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
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

  ScalarT scalarDeviatoricStrainInc;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunctionVal;

  ScalarT deltaLambda;

  double hardMod, hardModN, hardModNP1;

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

      // Compute \S_ij * \S_ij
      tempScalar = 0.0;
      for(int j = 0; j < 3; j++){
        for(int i = 0; i < 3; i++){
          tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
        }
      }

      // Avoid divide-by-zero
      deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));

      *vmStress = sqrt(3.0/2.0*tempScalar);

      //Evaluate the yield surface
      yieldFunctionVal = yieldStress * pow(1.0 + *eqpsN / hardeningCharacteristicStrain, hardeningExponent);

      // Elastic or plastic?
      if(*vmStress < yieldFunctionVal){
        // elastic it is
        *eqpsNP1 = *eqpsN;

      } else {
        // The step is plastic and we need to update the yield surface
        // (hardening) and return to that.


        // Compute \S_ij * \epsilon_inc_ij
        tempScalar = 0.0;
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStrainInc[i+3*j];
          }
        }

        scalarDeviatoricStrainInc = tempScalar / deviatoricStressMagnitudeNP1;

        // First go back to step N and compute deviatoric stress and its
        // magnitude.
        sphericalStressN = (*(stressN) + *(stressN+4) + *(stressN+8))/3.0;

        for (int i = 0; i < 9; i++){
          deviatoricStressN[i] = *(stressN+i);
        }
        deviatoricStressN[0] -= sphericalStressN;
        deviatoricStressN[4] -= sphericalStressN;
        deviatoricStressN[8] -= sphericalStressN;

        // Compute \S_ij * \S_ij
        tempScalar = 0.0;
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
          }
        }

        deviatoricStressMagnitudeN = sqrt(tempScalar);

        //Solve for deltaLambda
        //deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - sqrt(6.0) * yieldFunctionVal 
            //+ 3.0 * deviatoricStressMagnitudeN) / (2.0 * (hardModN + 3.0 * shearMod));
        //deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - 3.0 * 
            //deviatoricStressMagnitudeN  - sqrt(6.0) * yieldStress) /
            //(2.0 * (hardMod + 3.0 * shearMod));

        // evaluate the hardening modulus at the previous step
        hardModN = yieldStress * hardeningExponent / hardeningCharacteristicStrain * pow(1.0 + *eqpsN / hardeningCharacteristicStrain, hardeningExponent-1.0);

        //Solve for a projected deltaLambda
        // dot(eps_d)_ij Q_ij - dot(s)/2mu - dot(lambda) = 0
        // Q_ij = s_ij / ||s||
        deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - sqrt(6.0) * yieldFunctionVal 
                      + 3.0 * deviatoricStressMagnitudeN) / (2.0 * (hardModN + 3.0 * shearMod));

        //Increment the plastic strain for the purposes of evaluating the
        //yield surface
        *eqpsNP1 = *eqpsN + sqrt(2.0 / 3.0) * deltaLambda;

        // evaluate the hardening modulus at the projected steps
        hardModNP1 = yieldStress * hardeningExponent / hardeningCharacteristicStrain * pow(1.0 + *eqpsNP1 / hardeningCharacteristicStrain, hardeningExponent-1.0);

        // take the average of the two moduli 
        hardMod = 0.5 * (hardModN + hardModNP1);

        //Solve for deltaLambda again
        deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - sqrt(6.0) * yieldFunctionVal 
                      + 3.0 * deviatoricStressMagnitudeN) / (2.0 * (hardMod + 3.0 * shearMod));
        //deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - 3.0 * 
            //deviatoricStressMagnitudeN  - sqrt(6.0) * yieldStress) /
            //(2.0 * (hardMod + 3.0 * shearMod));

        //Increment the plastic strain for the purposes of evaluating the
        //yield surface
        *eqpsNP1 = *eqpsN + sqrt(2.0 / 3.0) * deltaLambda;

        //Evaluate the extent of the yield surface
        yieldFunctionVal = yieldStress * pow(1.0 + *eqpsNP1 / hardeningCharacteristicStrain, hardeningExponent);

        //If true, the step is plastic and we need to return to the yield
        //surface.  
        if(*vmStress > yieldFunctionVal){

          //multiply the trial deviatoric stress
          tempScalar = sqrt(2.0/3.0)*yieldFunctionVal/deviatoricStressMagnitudeNP1;

          // Return the deviatoric stress to the yield surface
          for(int i = 0; i < 9; i++){
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
          for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
              tempScalar += tempA[i+3*j] * tempB[i+3*j];
            }
          }

          // Increment the plastic strain
          *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);

        } else {
          // The step is elastic
          *eqpsNP1 = *eqpsN;
        }

      // compute stress triaxiality
      *triaxiality = (*(stressNP1+0) + *(stressNP1+4) + *(stressNP1+8)) / (3.0 * std::max(1.0e-20,*vmStress));
      }
    }
  }
}

template<typename ScalarT>
void updateElasticIsotropicPowerlawHardeningPlasticCauchyStress
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* cauchyStressN, 
    ScalarT* cauchyStressNP1, 
    ScalarT* vonMisesStress,
    const ScalarT* equivalentPlasticStrainN,
    ScalarT* equivalentPlasticStrainNP1,
    ScalarT* stressTriaxiality,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
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

  for(int iID=0 ; iID<numPoints ; ++iID, rateOfDef+=9, stressN+=9, stressNP1+=9,
      ++vmStress, ++eqpsN, ++eqpsNP1, ++triaxiality){

    radialReturnPowerlaw(rateOfDef, 
                         stressN, 
                         stressNP1, 
                         vmStress,
                         eqpsN, 
                         eqpsNP1, 
                         bulkMod, 
                         shearMod,
                         yieldStress, 
                         hardeningCharacteristicStrain,
                         hardeningExponent,
                         dt);

    // compute stress triaxiality
    *triaxiality = (*(stressNP1+0) + *(stressNP1+4) + *(stressNP1+8)) / (3.0 * std::max(1.0e-20,*vmStress));
  }
}

template<typename ScalarT>
void updateElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* cauchyStressN, 
    ScalarT* cauchyStressNP1, 
    ScalarT* vonMisesStress,
    const ScalarT* equivalentPlasticStrainN,
    ScalarT* equivalentPlasticStrainNP1,
    ScalarT* stressTriaxiality,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
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

  for(int iID=0 ; iID<numPoints ; ++iID, rateOfDef+=9, stressN+=9, stressNP1+=9,
      ++vmStress, ++eqpsN, ++eqpsNP1, ++triaxiality){

    radialReturnSaturationExponential(rateOfDef, 
                                      stressN, 
                                      stressNP1, 
                                      vmStress,
                                      eqpsN, 
                                      eqpsNP1, 
                                      bulkMod, 
                                      shearMod,
                                      initialYieldStress, 
                                      saturatedYieldStress,
                                      exponentialConstant, 
                                      linearConstant, 
                                      dt);
  
    // compute stress triaxiality
    *triaxiality = (*(stressNP1+0) + *(stressNP1+4) + *(stressNP1+8)) / (3.0 * std::max(1.0e-20,*vmStress));
  }
}

template<typename ScalarT>
void updateBondLevelElasticIsotropicPowerlawHardeningPlasticCauchyStress
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
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
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

  const double* flyingPointFlg = flyingPointFlag;

  const ScalarT* eqpsN = bondLevelEquivalentPlasticStrainN;
  ScalarT* eqpsNP1 = bondLevelEquivalentPlasticStrainNP1;
  ScalarT* triaxiality = bondLevelStressTriaxiality;

  ScalarT strainInc[9];
  ScalarT deviatoricStrainInc[9];
  ScalarT deviatoricStressN[9];
  ScalarT deviatoricStressMagnitudeN;

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;

  ScalarT scalarDeviatoricStrainInc;

  ScalarT tempA[9];
  ScalarT tempB[9];

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunctionVal;

  ScalarT deltaLambda;

  double hardMod, hardModN, hardModNP1;

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

        // Compute \S_ij * \S_ij
        tempScalar = 0.0;
        for(int j = 0; j < 3; j++){
          for(int i = 0; i < 3; i++){
            tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStressNP1[i+3*j];
          }
        }

        // Avoid divide-by-zero
        deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));

        *vmStress = sqrt(3.0/2.0*tempScalar);

        //Evaluate the yield surface
        yieldFunctionVal = yieldStress * pow(1.0 + *eqpsN / hardeningCharacteristicStrain, hardeningExponent);

        // Elastic or plastic?
        if(*vmStress < yieldFunctionVal){
          // elastic it is
          *eqpsNP1 = *eqpsN;

        } else {
          // The step is plastic and we need to update the yield surface
          // (hardening) and return to that.

          // Compute \S_ij * \epsilon_inc_ij
          tempScalar = 0.0;
          for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
              tempScalar += deviatoricStressNP1[i+3*j] * deviatoricStrainInc[i+3*j];
            }
          }

          scalarDeviatoricStrainInc = tempScalar / deviatoricStressMagnitudeNP1;

          // First go back to step N and compute deviatoric stress and its
          // magnitude.
          sphericalStressN = (*stressXXN + *stressYYN + *stressZZN)/3.0;

          deviatoricStressN[0] = *stressXXN; deviatoricStressN[1] = *stressXYN; deviatoricStressN[2] = *stressXZN;
          deviatoricStressN[3] = *stressYXN; deviatoricStressN[4] = *stressYYN; deviatoricStressN[5] = *stressYZN;
          deviatoricStressN[6] = *stressZXN; deviatoricStressN[7] = *stressZYN; deviatoricStressN[8] = *stressZZN;
          deviatoricStressN[0] -= sphericalStressN;
          deviatoricStressN[4] -= sphericalStressN;
          deviatoricStressN[8] -= sphericalStressN;

          // Compute \S_ij * \S_ij
          tempScalar = 0.0;
          for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
              tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
            }
          }

          deviatoricStressMagnitudeN = sqrt(tempScalar);

          // evaluate the hardening modulus at the previous step
          hardModN = yieldStress * hardeningExponent / hardeningCharacteristicStrain * pow(1.0 + *eqpsN / hardeningCharacteristicStrain, hardeningExponent-1.0);

          //Solve for a projected deltaLambda
          // dot(eps_d)_ij Q_ij - dot(s)/2mu - dot(lambda) = 0
          // Q_ij = s_ij / ||s||
          deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - sqrt(6.0) * yieldFunctionVal 
              + 3.0 * deviatoricStressMagnitudeN) / (2.0 * (hardModN + 3.0 * shearMod));
          //deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - 3.0 * 
              //deviatoricStressMagnitudeN  - sqrt(6.0) * yieldStress) /
              //(2.0 * (hardMod + 3.0 * shearMod));

          //Increment the plastic strain for the purposes of evaluating the
          //yield surface
          *eqpsNP1 = *eqpsN + sqrt(2.0 / 3.0) * deltaLambda;

          // evaluate the hardening modulus at the projected steps
          hardModNP1 = yieldStress * hardeningExponent / hardeningCharacteristicStrain * pow(1.0 + *eqpsNP1 / hardeningCharacteristicStrain, hardeningExponent-1.0);

          // take the average of the two moduli 
          hardMod = 0.5 * (hardModN + hardModNP1);

          //Solve for deltaLambda again
          deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - sqrt(6.0) * yieldFunctionVal 
              + 3.0 * deviatoricStressMagnitudeN) / (2.0 * (hardMod + 3.0 * shearMod));
          //deltaLambda = (6.0 * scalarDeviatoricStrainInc * shearMod - 3.0 * 
              //deviatoricStressMagnitudeN  - sqrt(6.0) * yieldStress) /
              //(2.0 * (hardMod + 3.0 * shearMod));

          //Increment the plastic strain for the purposes of evaluating the
          //yield surface
          *eqpsNP1 = *eqpsN + sqrt(2.0 / 3.0) * deltaLambda;

          //Evaluate the extent of the yield surface
          yieldFunctionVal = yieldStress * pow(1.0 + *eqpsNP1 / hardeningCharacteristicStrain, hardeningExponent);

          // Now return the stress to the yield surface

          //multiply the trial deviatoric stress
          tempScalar = sqrt(2.0/3.0)*yieldFunctionVal/deviatoricStressMagnitudeNP1;

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
          // magnitude.  We didn't do this earlier because it wouldn't be necassary
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
          for(int j = 0; j < 3; j++){
            for(int i = 0; i < 3; i++){
              tempScalar += tempA[i+3*j] * tempB[i+3*j];
            }
          }

          // Increment the plastic strain
          *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);
        }

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

template<typename ScalarT>
void updateBondLevelElasticIsotropicPowerlawHardeningPlasticCauchyStress
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
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
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

  std::vector<ScalarT> rateOfDefVector(9);
  ScalarT* rateOfDef = &rateOfDefVector[0];

  std::vector<ScalarT> stressNVector(9);
  ScalarT* stressN = &stressNVector[0];

  std::vector<ScalarT> stressNP1Vector(9);
  ScalarT* stressNP1 = &stressNP1Vector[0];

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

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

      *(rateOfDef+0) = *rateOfDefXX; *(rateOfDef+1) = *rateOfDefXY; *(rateOfDef+2) = *rateOfDefXZ; 
      *(rateOfDef+3) = *rateOfDefYX; *(rateOfDef+4) = *rateOfDefYY; *(rateOfDef+5) = *rateOfDefYZ; 
      *(rateOfDef+6) = *rateOfDefZX; *(rateOfDef+7) = *rateOfDefZY; *(rateOfDef+8) = *rateOfDefZZ; 
      *(stressN+0) = *stressXXN; *(stressN+1) = *stressXYN; *(stressN+2) = *stressXZN; 
      *(stressN+3) = *stressYXN; *(stressN+4) = *stressYYN; *(stressN+5) = *stressYZN; 
      *(stressN+6) = *stressZXN; *(stressN+7) = *stressZYN; *(stressN+8) = *stressZZN; 

      radialReturnPowerlaw(rateOfDef, 
                           stressN, 
                           stressNP1, 
                           vmStress,
                           eqpsN, 
                           eqpsNP1, 
                           bulkMod, 
                           shearMod,
                           yieldStress, 
                           hardeningCharacteristicStrain,
                           hardeningExponent,
                           dt);

      *stressXXNP1 = *(stressNP1+0); *stressXYNP1 = *(stressNP1+1); *stressXZNP1 = *(stressNP1+2); 
      *stressYXNP1 = *(stressNP1+3); *stressYYNP1 = *(stressNP1+4); *stressYZNP1 = *(stressNP1+5); 
      *stressZXNP1 = *(stressNP1+6); *stressZYNP1 = *(stressNP1+7); *stressZZNP1 = *(stressNP1+8); 

      // compute stress triaxiality
      *triaxiality = (*stressXXNP1 + *stressYYNP1 + *stressZZNP1) / (3.0 * std::max(1.0e-20,*vmStress));
    }
  }
}

template<typename ScalarT>
void updateBondLevelElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress
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
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
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

  std::vector<ScalarT> rateOfDefVector(9);
  ScalarT* rateOfDef = &rateOfDefVector[0];

  std::vector<ScalarT> stressNVector(9);
  ScalarT* stressN = &stressNVector[0];

  std::vector<ScalarT> stressNP1Vector(9);
  ScalarT* stressNP1 = &stressNP1Vector[0];

  int numNeighbors;
  const int *neighborListPtr = neighborhoodList;
  for(int iID=0 ; iID<numPoints ; ++iID){

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

      *(rateOfDef+0) = *rateOfDefXX; *(rateOfDef+1) = *rateOfDefXY; *(rateOfDef+2) = *rateOfDefXZ; 
      *(rateOfDef+3) = *rateOfDefYX; *(rateOfDef+4) = *rateOfDefYY; *(rateOfDef+5) = *rateOfDefYZ; 
      *(rateOfDef+6) = *rateOfDefZX; *(rateOfDef+7) = *rateOfDefZY; *(rateOfDef+8) = *rateOfDefZZ; 
      *(stressN+0) = *stressXXN; *(stressN+1) = *stressXYN; *(stressN+2) = *stressXZN; 
      *(stressN+3) = *stressYXN; *(stressN+4) = *stressYYN; *(stressN+5) = *stressYZN; 
      *(stressN+6) = *stressZXN; *(stressN+7) = *stressZYN; *(stressN+8) = *stressZZN; 

      radialReturnSaturationExponential(rateOfDef, 
                                        stressN, 
                                        stressNP1, 
                                        vmStress,
                                        eqpsN, 
                                        eqpsNP1, 
                                        bulkMod, 
                                        shearMod,
                                        initialYieldStress, 
                                        saturatedYieldStress,
                                        exponentialConstant, 
                                        linearConstant, 
                                        dt);

      // compute stress triaxiality
      *triaxiality = (*stressXXNP1 + *stressYYNP1 + *stressZZNP1) / (3.0 * std::max(1.0e-20,*vmStress));

      *stressXXNP1 = *(stressNP1+0); *stressXYNP1 = *(stressNP1+1); *stressXZNP1 = *(stressNP1+2); 
      *stressYXNP1 = *(stressNP1+3); *stressYYNP1 = *(stressNP1+4); *stressYZNP1 = *(stressNP1+5); 
      *stressZXNP1 = *(stressNP1+6); *stressZYNP1 = *(stressNP1+7); *stressZZNP1 = *(stressNP1+8); 
    }
  }
}

template<typename ScalarT>
void radialReturnPowerlaw
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* cauchyStressN, 
    ScalarT* cauchyStressNP1, 
    ScalarT* vonMisesStress,
    const ScalarT* equivalentPlasticStrainN,
    ScalarT* equivalentPlasticStrainNP1,
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
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

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;
  ScalarT flowVector[9];
  ScalarT vmStress_trial;

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunctionVal;
  ScalarT hardenedYieldStress;
  ScalarT yieldFunctionDerivative;
  ScalarT deltaPlasticIncrement;
  ScalarT plasticIncrement;

  double tol = 1.0e-16 * bulkMod;
  int counter;
  int maxIterations = 100;

  double hardMod;

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
  *(stressNP1+0) += bulkMod*dilatationInc;
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

  // Compute \S_ij * \S_ij
  tempScalar = 0.0;
  for(int i = 0; i < 9; i++)
    tempScalar += deviatoricStressNP1[i] * deviatoricStressNP1[i];

  // Avoid divide-by-zero
  deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));

  vmStress_trial = sqrt(3.0/2.0*tempScalar);

  //Evaluate the current yield stress
  hardenedYieldStress = yieldStress * pow(1.0 + *eqpsN / hardeningCharacteristicStrain, hardeningExponent);

  //Evaluate the yield function
  yieldFunctionVal = vmStress_trial - hardenedYieldStress;

  // Elastic or plastic?
  if(yieldFunctionVal < 0.0){

    // elastic it is
    *eqpsNP1 = *eqpsN;

    *vmStress = vmStress_trial;

  } else {
    // The step is plastic and we need to update the yield surface
    // (hardening) and return to that.
    
    counter = 0;
    plasticIncrement = 0.0;
    *eqpsNP1 = *eqpsN;

    // Flow vector's direction does not change in radial return: this is N = df/dsigma
    for (int i = 0; i < 9; i++)
      flowVector[i] = deviatoricStressNP1[i] / deviatoricStressMagnitudeNP1;

    // Newton iteration to project stress to the yield surface
    while(fabs(yieldFunctionVal) > tol){

      counter++;

      if(counter > maxIterations){
        std::cout << "isotropic_hardening_correspondence.cxx::radialReturnPowerLaw: Too many iterations in the plastic radial return." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // evaluate the hardening modulus 
      hardMod = yieldStress * hardeningExponent / hardeningCharacteristicStrain * pow(1.0 + *eqpsNP1 / hardeningCharacteristicStrain, hardeningExponent-1.0);

      // --- (1) COMPUTE INCREMENT IN THE PLASTIC MULTIPLIER ---
      // Residual derivative
      yieldFunctionDerivative = 3.0*shearMod + hardMod;

      // Plastic increment
      deltaPlasticIncrement = (yieldFunctionVal/yieldFunctionDerivative);

      plasticIncrement += deltaPlasticIncrement;

      // --- (2) UPDATE --- 
      // Update eq. plastic strain
      *eqpsNP1 += deltaPlasticIncrement;

      // Update stress
      tempScalar = sqrt(6.0) * shearMod * deltaPlasticIncrement;
      for (int i = 0; i < 9; i++){
        *(stressNP1+i) -= tempScalar * flowVector[i];
      }

      //Evaluate the current yield stress
      hardenedYieldStress = yieldStress * pow(1.0 + *eqpsNP1 / hardeningCharacteristicStrain, hardeningExponent);

      //Evaluate the yield function
      yieldFunctionVal = vmStress_trial - 3.0 * shearMod * plasticIncrement - hardenedYieldStress;
    } 

    // Update von Mises stress
    *vmStress = vmStress_trial - 3.0 * shearMod * plasticIncrement;
  }
}

template<typename ScalarT>
void radialReturnSaturationExponential
(
    const ScalarT* unrotatedRateOfDeformation, 
    const ScalarT* cauchyStressN, 
    ScalarT* cauchyStressNP1, 
    ScalarT* vonMisesStress,
    const ScalarT* equivalentPlasticStrainN,
    ScalarT* equivalentPlasticStrainNP1,
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
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

  ScalarT deviatoricStressNP1[9];
  ScalarT deviatoricStressMagnitudeNP1;
  ScalarT flowVector[9];
  ScalarT vmStress_trial;

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT yieldFunctionVal;
  ScalarT hardenedYieldStress;
  ScalarT yieldFunctionDerivative;
  ScalarT deltaPlasticIncrement;
  ScalarT plasticIncrement;

  double tol = 1.0e-16 * bulkMod;
  int counter;
  int maxIterations = 100;

  double hardMod;

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
  *(stressNP1+0) += bulkMod*dilatationInc;
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

  // Compute \S_ij * \S_ij
  tempScalar = 0.0;
  for(int i = 0; i < 9; i++)
    tempScalar += deviatoricStressNP1[i] * deviatoricStressNP1[i];

  // Avoid divide-by-zero
  deviatoricStressMagnitudeNP1 = std::max(1.0e-20,sqrt(tempScalar));

  vmStress_trial = sqrt(3.0/2.0*tempScalar);

  //Evaluate the current yield stress
  hardenedYieldStress = initialYieldStress + (saturatedYieldStress - initialYieldStress) * ( 1.0 - exp(-exponentialConstant * *eqpsN)) + linearConstant * *eqpsN;

  //Evaluate the yield function
  yieldFunctionVal = vmStress_trial - hardenedYieldStress;

  // Elastic or plastic?
  if(yieldFunctionVal < 0.0){

    // elastic it is
    *eqpsNP1 = *eqpsN;

    *vmStress = vmStress_trial;

  } else {
    // The step is plastic and we need to update the yield surface
    // (hardening) and return to that.
    
    counter = 0;
    plasticIncrement = 0.0;
    *eqpsNP1 = *eqpsN;

    // Flow vector's direction does not change in radial return: this is N = df/dsigma
    for (int i = 0; i < 9; i++)
      flowVector[i] = deviatoricStressNP1[i] / deviatoricStressMagnitudeNP1;

    // Newton iteration to project stress to the yield surface
    while(fabs(yieldFunctionVal) > tol){

      counter++;

      if(counter > maxIterations){
        std::cout << "isotropic_hardening_correspondence.cxx::radialReturnSaturationExponential: Too many iterations in the plastic radial return." << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
      }

      // evaluate the hardening modulus 
      hardMod = exponentialConstant * (saturatedYieldStress - initialYieldStress) * exp(-exponentialConstant * *eqpsNP1) + linearConstant ;

      // --- (1) COMPUTE INCREMENT IN THE PLASTIC MULTIPLIER ---
      // Residual derivative
      yieldFunctionDerivative = 3.0*shearMod + hardMod;

      // Plastic increment
      deltaPlasticIncrement = (yieldFunctionVal/yieldFunctionDerivative);

      plasticIncrement += deltaPlasticIncrement;

      // --- (2) UPDATE --- 
      // Update eq. plastic strain
      *eqpsNP1 += deltaPlasticIncrement;

      // Update stress
      tempScalar = sqrt(6.0) * shearMod * deltaPlasticIncrement;
      for (int i = 0; i < 9; i++){
        *(stressNP1+i) -= tempScalar * flowVector[i];
      }

      //Evaluate the current yield stress
      hardenedYieldStress = initialYieldStress + (saturatedYieldStress - initialYieldStress) * ( 1.0 - exp(-exponentialConstant * *eqpsNP1)) + linearConstant * *eqpsNP1;

      //Evaluate the yield function
      yieldFunctionVal = vmStress_trial - 3.0 * shearMod * plasticIncrement - hardenedYieldStress;

    } 

    // Update von Mises stress
    *vmStress = vmStress_trial - 3.0 * shearMod * plasticIncrement;
  }
}

template void updateElasticIsotropicPowerlawHardeningPlasticCauchyStress<double>
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
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
    const double dt
);

template void updateElasticIsotropicPowerlawHardeningPlasticCauchyStress<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    double* stressTriaxiality,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
    const double dt
);

template void updateElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    double* stressTriaxiality,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
    const double dt
);

template void updateBondLevelElasticIsotropicPowerlawHardeningPlasticCauchyStress<double>
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
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
    const double dt
);

template void updateBondLevelElasticIsotropicPowerlawHardeningPlasticCauchyStress<double>
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
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
    const double dt
);

template void updateElasticIsotropicHardeningPlasticCauchyStress<double>
(
    const double* modelCoordinates,
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
    const double hardMod,
    const bool m_isFlaw,
    const double m_flawLocationX,
    const double m_flawLocationY,
    const double m_flawLocationZ,
    const double m_flawSize,
    const double m_flawMagnitude,
    const double dt
);

template void updateBondLevelElasticIsotropicSaturationExponentialHardeningPlasticCauchyStress<double>
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
    const int* neighborhoodList,
    const int numPoints, 
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
    const double dt
);

template void radialReturnPowerlaw<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    const double bulkMod,
    const double shearMod,
    const double yieldStress,
    const double hardeningCharacteristicStrain,
    const double hardeningExponent,
    const double dt
);

template void radialReturnSaturationExponential<double>
(
    const double* unrotatedRateOfDeformation, 
    const double* cauchyStressN, 
    double* cauchyStressNP1, 
    double* vonMisesStress,
    const double* equivalentPlasticStrainN,
    double* equivalentPlasticStrainNP1,
    const double bulkMod,
    const double shearMod,
    const double initialYieldStress,
    const double saturatedYieldStress,
    const double exponentialConstant,
    const double linearConstant,
    const double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void updateElasticIsotropicHardeningPlasticCauchyStress<Sacado::Fad::DFad<double> >
(
    const double* modelCoordinates,
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
    const double hardMod,
    const bool m_isFlaw,
    const double m_flawLocationX,
    const double m_flawLocationY,
    const double m_flawLocationZ,
    const double m_flawSize,
    const double m_flawMagnitude,
    const double dt
);

}
