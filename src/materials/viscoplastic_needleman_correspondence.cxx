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

#include "viscoplastic_needleman_correspondence.h"
#include "correspondence.h"
#include "material_utilities.h"
#include <Sacado.hpp>
#include <math.h>

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticViscoplasticCauchyStress
(
const double* modelCoordinates,
const ScalarT* unrotatedRateOfDeformationXX, 
const ScalarT* unrotatedRateOfDeformationXY, 
const ScalarT* unrotatedRateOfDeformationXZ, 
const ScalarT* unrotatedRateOfDeformationYX, 
const ScalarT* unrotatedRateOfDeformationYY, 
const ScalarT* unrotatedRateOfDeformationYZ, 
const ScalarT* unrotatedRateOfDeformationZX, 
const ScalarT* unrotatedRateOfDeformationZY, 
const ScalarT* unrotatedRateOfDeformationZZ, 
const ScalarT* cauchyStressNXX, 
const ScalarT* cauchyStressNXY, 
const ScalarT* cauchyStressNXZ, 
const ScalarT* cauchyStressNYX, 
const ScalarT* cauchyStressNYY, 
const ScalarT* cauchyStressNYZ, 
const ScalarT* cauchyStressNZX, 
const ScalarT* cauchyStressNZY, 
const ScalarT* cauchyStressNZZ, 
ScalarT* cauchyStressNP1XX, 
ScalarT* cauchyStressNP1XY, 
ScalarT* cauchyStressNP1XZ, 
ScalarT* cauchyStressNP1YX, 
ScalarT* cauchyStressNP1YY, 
ScalarT* cauchyStressNP1YZ, 
ScalarT* cauchyStressNP1ZX, 
ScalarT* cauchyStressNP1ZY, 
ScalarT* cauchyStressNP1ZZ, 
ScalarT* vonMisesStress,
const ScalarT* equivalentPlasticStrainN,
ScalarT* equivalentPlasticStrainNP1,
const int numPoints, 
const double bulkMod,
const double shearMod,
const double yieldStress,
const double strainHardExp,
const double rateHardExp, 
const double refStrainRate,
const double refStrain0,
const double refStrain1,
const bool isFlaw,
const double flawLocationX,
const double flawLocationY,
const double flawLocationZ,
const double flawSize,
const double flawMagnitude,
const double dt
)
{
  
  const ScalarT* rateOfDefXX = unrotatedRateOfDeformationXX;
  const ScalarT* rateOfDefXY = unrotatedRateOfDeformationXY;
  const ScalarT* rateOfDefXZ = unrotatedRateOfDeformationXZ;
  const ScalarT* rateOfDefYX = unrotatedRateOfDeformationYX;
  const ScalarT* rateOfDefYY = unrotatedRateOfDeformationYY;
  const ScalarT* rateOfDefYZ = unrotatedRateOfDeformationYZ;
  const ScalarT* rateOfDefZX = unrotatedRateOfDeformationZX;
  const ScalarT* rateOfDefZY = unrotatedRateOfDeformationZY;
  const ScalarT* rateOfDefZZ = unrotatedRateOfDeformationZZ;
  const ScalarT* stressNXX = cauchyStressNXX;
  const ScalarT* stressNXY = cauchyStressNXY;
  const ScalarT* stressNXZ = cauchyStressNXZ;
  const ScalarT* stressNYX = cauchyStressNYX;
  const ScalarT* stressNYY = cauchyStressNYY;
  const ScalarT* stressNYZ = cauchyStressNYZ;
  const ScalarT* stressNZX = cauchyStressNZX;
  const ScalarT* stressNZY = cauchyStressNZY;
  const ScalarT* stressNZZ = cauchyStressNZZ;
  ScalarT* stressNP1XX = cauchyStressNP1XX;
  ScalarT* stressNP1XY = cauchyStressNP1XY;
  ScalarT* stressNP1XZ = cauchyStressNP1XZ;
  ScalarT* stressNP1YX = cauchyStressNP1YX;
  ScalarT* stressNP1YY = cauchyStressNP1YY;
  ScalarT* stressNP1YZ = cauchyStressNP1YZ;
  ScalarT* stressNP1ZX = cauchyStressNP1ZX;
  ScalarT* stressNP1ZY = cauchyStressNP1ZY;
  ScalarT* stressNP1ZZ = cauchyStressNP1ZZ;

  ScalarT* vmStress = vonMisesStress;

  const ScalarT* eqpsN = equivalentPlasticStrainN;
  ScalarT* eqpsNP1 = equivalentPlasticStrainNP1;

  ScalarT strainIncXX, strainIncXY, strainIncXZ;
  ScalarT strainIncYX, strainIncYY, strainIncYZ;
  ScalarT strainIncZX, strainIncZY, strainIncZZ;

  ScalarT deviatoricStrainIncXX, deviatoricStrainIncXY, deviatoricStrainIncXZ;
  ScalarT deviatoricStrainIncYX, deviatoricStrainIncYY, deviatoricStrainIncYZ;
  ScalarT deviatoricStrainIncZX, deviatoricStrainIncZY, deviatoricStrainIncZZ;
  
  ScalarT deviatoricStressNXX, deviatoricStressNXY, deviatoricStressNXZ;
  ScalarT deviatoricStressNYX, deviatoricStressNYY, deviatoricStressNYZ;
  ScalarT deviatoricStressNZX, deviatoricStressNZY, deviatoricStressNZZ;
  ScalarT deviatoricStressMagnitudeN;

  ScalarT deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ;
  ScalarT deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ;
  ScalarT deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ;
  ScalarT deviatoricStressMagnitudeNP1;

  ScalarT scalarDeviatoricStrainInc;

  ScalarT tempAXX, tempAXY, tempAXZ;
  ScalarT tempAYX, tempAYY, tempAYZ;
  ScalarT tempAZX, tempAZY, tempAZZ;

  ScalarT tempBXX, tempBXY, tempBXZ;
  ScalarT tempBYX, tempBYY, tempBYZ;
  ScalarT tempBZX, tempBZY, tempBZZ;

  ScalarT dilatationInc;
  ScalarT sphericalStressN;
  ScalarT sphericalStressNP1;
  ScalarT tempScalar;
  ScalarT scaleFactor1;
  ScalarT scaleFactor2;
  ScalarT yieldFunctionVal;

  ScalarT rateTerm;
  ScalarT hardTerm;

  ScalarT deltaLambda;
  ScalarT deltaLambdaOld = 0.0;
  ScalarT a;
  ScalarT b;
  ScalarT c;
  ScalarT fb;
  ScalarT fc;
  ScalarT yfb;
  ScalarT yfc;

  double reducedYieldStress;
  const double* modelCoord = modelCoordinates;


  for(int iID=0 ; iID<numPoints ; ++iID, modelCoord+=3, 
        ++rateOfDefXX, ++rateOfDefXY, ++rateOfDefXZ,
        ++rateOfDefYX, ++rateOfDefYY, ++rateOfDefYZ,
        ++rateOfDefZX, ++rateOfDefZY, ++rateOfDefZZ,
        ++stressNXX, ++stressNXY, ++stressNXZ,
        ++stressNYX, ++stressNYY, ++stressNYZ,
        ++stressNZX, ++stressNZY, ++stressNZZ,
        ++stressNP1XX, ++stressNP1XY, ++stressNP1XZ,
        ++stressNP1YX, ++stressNP1YY, ++stressNP1YZ,
        ++stressNP1ZX, ++stressNP1ZY, ++stressNP1ZZ,
        ++vmStress,++eqpsN,++eqpsNP1){

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

      // Compute the elastic ``trial'' stress
      *stressNP1XX = (*stressNXX) + 2.0*shearMod*deviatoricStrainIncXX + bulkMod*dilatationInc;
      *stressNP1XY = (*stressNXY) + 2.0*shearMod*deviatoricStrainIncXY;
      *stressNP1XZ = (*stressNXZ) + 2.0*shearMod*deviatoricStrainIncXZ;
      *stressNP1YX = (*stressNYX) + 2.0*shearMod*deviatoricStrainIncYX;
      *stressNP1YY = (*stressNYY) + 2.0*shearMod*deviatoricStrainIncYY + bulkMod*dilatationInc;
      *stressNP1YZ = (*stressNYZ) + 2.0*shearMod*deviatoricStrainIncYZ;
      *stressNP1ZX = (*stressNZX) + 2.0*shearMod*deviatoricStrainIncZX;
      *stressNP1ZY = (*stressNZY) + 2.0*shearMod*deviatoricStrainIncZY;
      *stressNP1ZZ = (*stressNZZ) + 2.0*shearMod*deviatoricStrainIncZZ + bulkMod*dilatationInc;


      sphericalStressNP1 = ((*stressNP1XX) + (*stressNP1YY) + (*stressNP1ZZ))/3.0;

      // Compute the ``trial'' von Mises stress
      deviatoricStressNP1XX = *stressNP1XX - sphericalStressNP1;
      deviatoricStressNP1XY = *stressNP1XY;
      deviatoricStressNP1XZ = *stressNP1XZ;
      deviatoricStressNP1YX = *stressNP1YX;
      deviatoricStressNP1YY = *stressNP1YY - sphericalStressNP1;
      deviatoricStressNP1YZ = *stressNP1YZ;
      deviatoricStressNP1ZX = *stressNP1ZX;
      deviatoricStressNP1ZY = *stressNP1ZY;
      deviatoricStressNP1ZZ = *stressNP1ZZ - sphericalStressNP1;

      tempScalar = TensorContraction(deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ, 
                                     deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ, 
                                     deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ, 
                                     deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ, 
                                     deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ, 
                                     deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ); 

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
          // function is rate dependent and can change over a load step.

          tempScalar = TensorContraction(deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ, 
                                         deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ, 
                                         deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ, 
                                         deviatoricStrainIncXX, deviatoricStrainIncXY, deviatoricStrainIncXZ, 
                                         deviatoricStrainIncYX, deviatoricStrainIncYY, deviatoricStrainIncYZ, 
                                         deviatoricStrainIncZX, deviatoricStrainIncZY, deviatoricStrainIncZZ); 

          scalarDeviatoricStrainInc = tempScalar / deviatoricStressMagnitudeNP1;

          // First go back to step N and compute deviatoric stress and its
          // magnitude.
          sphericalStressN = ((*stressNXX) + (*stressNYY) + (*stressNZZ))/3.0;

          deviatoricStressNXX = *stressNXX - sphericalStressN;
          deviatoricStressNXY = *stressNXY;
          deviatoricStressNXZ = *stressNXZ;
          deviatoricStressNYX = *stressNYX;
          deviatoricStressNYY = *stressNYY - sphericalStressN;
          deviatoricStressNYZ = *stressNYZ;
          deviatoricStressNZX = *stressNZX;
          deviatoricStressNZY = *stressNZY;
          deviatoricStressNZZ = *stressNZZ - sphericalStressN;

          tempScalar = TensorContraction(deviatoricStressNXX, deviatoricStressNXY, deviatoricStressNXZ, 
                                         deviatoricStressNYX, deviatoricStressNYY, deviatoricStressNYZ, 
                                         deviatoricStressNZX, deviatoricStressNZY, deviatoricStressNZZ, 
                                         deviatoricStressNXX, deviatoricStressNXY, deviatoricStressNXZ, 
                                         deviatoricStressNYX, deviatoricStressNYY, deviatoricStressNYZ, 
                                         deviatoricStressNZX, deviatoricStressNZY, deviatoricStressNZZ); 

          deviatoricStressMagnitudeN = sqrt(tempScalar);

          //Verbose bi-section method to solve for deltaLambda
          //TODO: Implement a ``smarter'' faster root solve here.
          deltaLambda = 1.0;
          deltaLambdaOld = 0.0;
          a = 0.0;
          b = 1.0;
          c = (a+b)/2.0;

          //Bisection loop
          for(int iter = 0; iter < 100000; iter++){

              *eqpsNP1 = *eqpsN + b;
              yfb = ViscoplasticNeedlemanYieldFunction(b, *eqpsNP1, reducedYieldStress, strainHardExp, rateHardExp, refStrainRate,refStrain0, refStrain1, dt);
              *eqpsNP1 = *eqpsN + c;
              yfc = ViscoplasticNeedlemanYieldFunction(c, *eqpsNP1, reducedYieldStress, strainHardExp, rateHardExp, refStrainRate,refStrain0, refStrain1, dt);
       
              fb = scalarDeviatoricStrainInc - b - 1.0 / 2.0 / shearMod * (sqrt(2.0/3.0) * yfb - deviatoricStressMagnitudeN);
              fc = scalarDeviatoricStrainInc - c - 1.0 / 2.0 / shearMod * (sqrt(2.0/3.0) * yfc - deviatoricStressMagnitudeN);
              
              if(fb > 0.0 && fc > 0.0) {
                  b = c;
              } 
              else if(fb < 0.0 && fc < 0.0){
                  b = c;
              }
              else{
                  a = c;
              }

              deltaLambdaOld = c;
              c = (a+b)/2.0;
              deltaLambda = c;

              if(fabs(deltaLambda - deltaLambdaOld)/fabs(deltaLambda) < 1.0e-6){
                  //We're converged, stop
                  break;
              }
              else if (iter == 99999){
                  //Error message here!
                  std::cout << "Bisection method failed to converge in 1e6 iterations" << std::endl;
                  break;
              }
          }

          //Increment the plastic strain for the purposes of evaluating the
          //yield surface
          *eqpsNP1 = *eqpsN + deltaLambda;
          //Evaluate the extent of the yield surface with the result of
          //ViscoplasticNeedlemanFindRoot
          yieldFunctionVal = ViscoplasticNeedlemanYieldFunction(deltaLambda, 
                                                                *eqpsNP1, 
                                                                reducedYieldStress, 
                                                                strainHardExp, 
                                                                rateHardExp, 
                                                                refStrainRate, 
                                                                refStrain0, 
                                                                refStrain1, 
                                                                dt);

          //If true, the step is plastic and we need to return to the yield
          //surface.  
          if(*vmStress > yieldFunctionVal){

              //For perfectly plastic deformations we just have a constant factor to 
              //multiply the trial deviatoric stress
              tempScalar = sqrt(2.0/3.0)*yieldFunctionVal/deviatoricStressMagnitudeNP1;

              // Return the deviatoric stress to the yield surface
              deviatoricStressNP1XX *= tempScalar; 
              deviatoricStressNP1XY *= tempScalar;
              deviatoricStressNP1XZ *= tempScalar;
              deviatoricStressNP1YX *= tempScalar;
              deviatoricStressNP1YY *= tempScalar;
              deviatoricStressNP1YZ *= tempScalar;
              deviatoricStressNP1ZX *= tempScalar;
              deviatoricStressNP1ZY *= tempScalar;
              deviatoricStressNP1ZZ *= tempScalar;

              // Update the Cauchy Stress
              *stressNP1XX = deviatoricStressNP1XX + sphericalStressNP1;
              *stressNP1XY = deviatoricStressNP1XY;
              *stressNP1XZ = deviatoricStressNP1XZ;
              *stressNP1YX = deviatoricStressNP1YX;
              *stressNP1YY = deviatoricStressNP1YY + sphericalStressNP1;
              *stressNP1YZ = deviatoricStressNP1YZ;
              *stressNP1ZX = deviatoricStressNP1ZX;
              *stressNP1ZY = deviatoricStressNP1ZY;
              *stressNP1ZZ = deviatoricStressNP1ZZ + sphericalStressNP1;

              // Update the von Mises stress now that the state of stress is on the
              // yield surface
              tempScalar = TensorContraction(deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ, 
                                             deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ, 
                                             deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ, 
                                             deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ, 
                                             deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ, 
                                             deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ); 

              *vmStress = sqrt(3.0/2.0*tempScalar);

              /////////////////////////////////////////////////////////////////////
              //
              // Update the equivalent plastic strain
              //
              // The algorithm below is generic and should not need to be modified for
              // any J2 plasticity yield surface.  It uses the difference in the yield
              // surface location at the NP1 and N steps to increment eqps regardless
              // of how the plastic multiplier was found in the yield surface
              // evaluation.
              //
              /////////////////////////////////////////////////////////////////////
              //
              //Ensure that this is at least a very small number to avoid a divide by
              //zero
              deviatoricStressMagnitudeN = std::max(1.0e-20,deviatoricStressMagnitudeN);

              //The plastic deviatoric strain increment tensor \Delta e_{plastic} = \Delta e_{total} - \Delta e_{elastic}
              tempAXX = deviatoricStrainIncXX - (deviatoricStressNP1XX - deviatoricStressNXX) / 2.0 / shearMod; 
              tempAXY = deviatoricStrainIncXY - (deviatoricStressNP1XY - deviatoricStressNXY) / 2.0 / shearMod; 
              tempAXZ = deviatoricStrainIncXZ - (deviatoricStressNP1XZ - deviatoricStressNXZ) / 2.0 / shearMod; 
              tempAYX = deviatoricStrainIncYX - (deviatoricStressNP1YX - deviatoricStressNYX) / 2.0 / shearMod; 
              tempAYY = deviatoricStrainIncYY - (deviatoricStressNP1YY - deviatoricStressNYY) / 2.0 / shearMod; 
              tempAYZ = deviatoricStrainIncYZ - (deviatoricStressNP1YZ - deviatoricStressNYZ) / 2.0 / shearMod; 
              tempAZX = deviatoricStrainIncZX - (deviatoricStressNP1ZX - deviatoricStressNZX) / 2.0 / shearMod; 
              tempAZY = deviatoricStrainIncZY - (deviatoricStressNP1ZY - deviatoricStressNZY) / 2.0 / shearMod; 
              tempAZZ = deviatoricStrainIncZZ - (deviatoricStressNP1ZZ - deviatoricStressNZZ) / 2.0 / shearMod; 

              //Deviatoric stress increment.  This is effectively an average of the deviatoric stress 
              //direction unit tensors at the half-step between steps NP1 and N
              scaleFactor1 = 1.0/deviatoricStressMagnitudeNP1/2.0;
              scaleFactor2 = 1.0/deviatoricStressMagnitudeN/2.0;
              MatrixUpdate(scaleFactor1, scaleFactor2,
                           deviatoricStressNP1XX, deviatoricStressNP1XY, deviatoricStressNP1XZ,
                           deviatoricStressNP1YX, deviatoricStressNP1YY, deviatoricStressNP1YZ,
                           deviatoricStressNP1ZX, deviatoricStressNP1ZY, deviatoricStressNP1ZZ,
                           deviatoricStressNXX, deviatoricStressNXY, deviatoricStressNXZ,
                           deviatoricStressNYX, deviatoricStressNYY, deviatoricStressNYZ,
                           deviatoricStressNZX, deviatoricStressNZY, deviatoricStressNZZ,
                           tempBXX, tempBXY, tempBXZ,
                           tempBYX, tempBYY, tempBYZ,
                           tempBZX, tempBZY, tempBZZ);

              // Contract the two tensors. This represents a projection of the plastic
              // strain increment tensor onto the "direction" of deviatoric stress
              // increment
              tempScalar = TensorContraction(tempAXX, tempAXY, tempAXZ, 
                                             tempAYX, tempAYY, tempAYZ, 
                                             tempAZX, tempAZY, tempAZZ, 
                                             tempBXX, tempBXY, tempBXZ, 
                                             tempBYX, tempBYY, tempBYZ, 
                                             tempBZX, tempBZY, tempBZZ); 

              // Increment the plastic strain
              *eqpsNP1 = *eqpsN + std::max(0.0, sqrt(2.0/3.0) * tempScalar);

          } else {
              // The step is elastic
              *eqpsNP1 = *eqpsN;
          }
      }
  }
}

// Needleman yield function
template <typename ScalarT>
ScalarT ViscoplasticNeedlemanYieldFunction
(
 const ScalarT deltaLambda,
 const ScalarT eqps,
 const double yieldStress,
 const double strainHardExp,
 const double rateHardExp, 
 const double refStrainRate,
 const double refStrain0,
 const double refStrain1,
 const double dt
)
{
  ScalarT hardTerm = yieldStress * ((1.0 + pow(eqps/refStrain0, strainHardExp)) / (1.0 + pow(eqps/refStrain1,2.0)));
  ScalarT rateTerm = pow(sqrt(2.0/3.0) * deltaLambda / dt / refStrainRate, rateHardExp);

  return (hardTerm * rateTerm);
}


template <typename ScalarT>
ScalarT ViscoplasticNeedlemanFindRoot
(
 const ScalarT eqps,
 const ScalarT scalarDeviatoricStrainInc,
 const ScalarT deviatoricStressMagnitude,
 const double yieldStress,
 const double shearMod,
 const double strainHardExp,
 const double rateHardExp, 
 const double refStrainRate,
 const double refStrain0,
 const double refStrain1,
 const double dt
)
//    FindRoot seeks the root of a function F(X) in an interval [A,B].
//
//  Discussion:
//
//    The interval [A,B] must be a change of sign interval for F.
//    That is, F(A) and F(B) must be of opposite signs.  Then
//    assuming that F is continuous implies the existence of at least
//    one value C between A and B for which F(C) = 0.
//
//    The location of the zero is determined to within an accuracy
//    of 6 * MACHEPS * abs ( C ) + 2 * T.
//
//  Translated from the FORTRAN code in:
//
//    Richard Brent,
//    Algorithms for Minimization Without Derivatives,
//    Dover, 2002,
//    ISBN: 0-486-41998-3,
//    LC: QA402.5.B74.
//
{
  ScalarT c;
  ScalarT d;
  ScalarT e;
  ScalarT fa;
  ScalarT fb;
  ScalarT yfa;
  ScalarT yfb;
  ScalarT fc;
  ScalarT m;
  ScalarT macheps;
  ScalarT p;
  ScalarT q;
  ScalarT r;
  ScalarT s;
  ScalarT sa;
  ScalarT sb;
  ScalarT tol;

  ScalarT A = 0.0;
  ScalarT B = 1.0;
  ScalarT tolerence = 1.0e-8;
//
//  Make local copies of A and B.
//
  sa = A;
  sb = B;
  yfa = ViscoplasticNeedlemanYieldFunction(sa, eqps, yieldStress, strainHardExp, rateHardExp, refStrainRate,refStrain0, refStrain1, dt);
  yfb = ViscoplasticNeedlemanYieldFunction(sb, eqps, yieldStress, strainHardExp, rateHardExp, refStrainRate,refStrain0, refStrain1, dt);
   
  fa = scalarDeviatoricStrainInc - sa - 1.0/2.0/shearMod*( sqrt(2.0/3.0)*yfa - deviatoricStressMagnitude);
  fb = scalarDeviatoricStrainInc - sb - 1.0/2.0/shearMod*( sqrt(2.0/3.0)*yfb - deviatoricStressMagnitude);

  c = sa;
  fc = fa;
  e = sb - sa;
  d = e;

  macheps = Teuchos::ScalarTraits<double>::eps();

  for ( ; ; )
  {
    if ( fabs ( fc ) < fabs ( fb ) )
    {
      sa = sb;
      sb = c;
      c = sa;
      fa = fb;
      fb = fc;
      fc = fa;
    }

    tol = 2.0 * macheps * fabs ( sb ) + tolerence;
    m = 0.5 * ( c - sb );

    if ( fabs ( m ) <= tol || fb == 0.0 )
    {
      break;
    }

    if ( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
    {
      e = m;
      d = e;
    }
    else
    {
      s = fb / fa;

      if ( sa == c )
      {
        p = 2.0 * m * s;
        q = 1.0 - s;
      }
      else
      {
        q = fa / fc;
        r = fb / fc;
        p = s * ( 2.0 * m * q * ( q - r ) - ( sb - sa ) * ( r - 1.0 ) );
        q = ( q - 1.0 ) * ( r - 1.0 ) * ( s - 1.0 );
      }

      if ( 0.0 < p )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      if ( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
        p < fabs ( 0.5 * s * q ) )
      {
        d = p / q;
      }
      else
      {
        e = m;
        d = e;
      }
    }
    sa = sb;
    fa = fb;

    if ( tol < fabs ( d ) )
    {
      sb = sb + d;
    }
    else if ( 0.0 < m )
    {
      sb = sb + tol;
    }
    else
    {
      sb = sb - tol;
    }

    yfb = ViscoplasticNeedlemanYieldFunction( sb, eqps, yieldStress, strainHardExp, rateHardExp, refStrainRate,refStrain0, refStrain1, dt);
    fb = scalarDeviatoricStrainInc - sb - 1.0/2.0/shearMod*( sqrt(2.0/3.0)*yfb - deviatoricStressMagnitude);

    if ( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
    {
      c = sa;
      fc = fa;
      e = sb - sa;
      d = e;
    }
  }
  return sb;
}


// Explicit template instantiation for double
template void updateElasticViscoplasticCauchyStress<double>
(
const double* modelCoordinates,
const double* unrotatedRateOfDeformationXX, 
const double* unrotatedRateOfDeformationXY, 
const double* unrotatedRateOfDeformationXZ, 
const double* unrotatedRateOfDeformationYX, 
const double* unrotatedRateOfDeformationYY, 
const double* unrotatedRateOfDeformationYZ, 
const double* unrotatedRateOfDeformationZX, 
const double* unrotatedRateOfDeformationZY, 
const double* unrotatedRateOfDeformationZZ, 
const double* cauchyStressNXX, 
const double* cauchyStressNXY, 
const double* cauchyStressNXZ, 
const double* cauchyStressNYX, 
const double* cauchyStressNYY, 
const double* cauchyStressNYZ, 
const double* cauchyStressNZX, 
const double* cauchyStressNZY, 
const double* cauchyStressNZZ, 
double* cauchyStressNP1XX, 
double* cauchyStressNP1XY, 
double* cauchyStressNP1XZ, 
double* cauchyStressNP1YX, 
double* cauchyStressNP1YY, 
double* cauchyStressNP1YZ, 
double* cauchyStressNP1ZX, 
double* cauchyStressNP1ZY, 
double* cauchyStressNP1ZZ, 
double* vonMisesStress,
const double* equivalentPlasticStrainN,
double* equivalentPlasticStrainNP1,
const int numPoints, 
const double bulkMod,
const double shearMod,
const double yieldStress,
const double strainHardExp,
const double rateHardExp, 
const double refStrainRate,
const double refStrain0,
const double refStrain1,
const bool m_isFlaw,
const double m_flawLocationX,
const double m_flawLocationY,
const double m_flawLocationZ,
const double m_flawSize,
const double m_flawMagnitude,
const double dt
);

template double ViscoplasticNeedlemanYieldFunction<double>
(
 const double deltaLambda,
 const double eqps,
 const double yieldStress,
 const double strainHardExp,
 const double rateHardExp, 
 const double refStrainRate,
 const double refStrain0,
 const double refStrain1,
 const double dt
);

template double ViscoplasticNeedlemanFindRoot<double>
(
 const double eqps,
 const double scalarDeviatoricStrainInc,
 const double deviatoricStressMagnitude,
 const double yieldStress,
 const double shearMod,
 const double strainHardExp,
 const double rateHardExp, 
 const double refStrainRate,
 const double refStrain0,
 const double refStrain1,
 const double dt
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template Sacado::Fad::DFad<double> ViscoplasticNeedlemanFindRoot<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double> eqps,
 const Sacado::Fad::DFad<double> scalarDeviatoricStrainInc,
 const Sacado::Fad::DFad<double> deviatoricStressMagnitude,
 const double yieldStress,
 const double shearMod,
 const double strainHardExp,
 const double rateHardExp, 
 const double refStrainRate,
 const double refStrain0,
 const double refStrain1,
 const double dt
);

template Sacado::Fad::DFad<double> ViscoplasticNeedlemanYieldFunction<Sacado::Fad::DFad<double> >
(
 const Sacado::Fad::DFad<double> deltaLambda,
 const Sacado::Fad::DFad<double> eqps,
 const double yieldStress,
 const double strainHardExp,
 const double rateHardExp, 
 const double refStrainRate,
 const double refStrain0,
 const double refStrain1,
 const double dt
);

template void updateElasticViscoplasticCauchyStress<Sacado::Fad::DFad<double> >
(
const double* modelCoordinates,
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXX, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXY, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationXZ, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYX, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYY, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationYZ, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZX, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZY, 
const Sacado::Fad::DFad<double>* unrotatedRateOfDeformationZZ, 
const Sacado::Fad::DFad<double>* cauchyStressNXX, 
const Sacado::Fad::DFad<double>* cauchyStressNXY, 
const Sacado::Fad::DFad<double>* cauchyStressNXZ, 
const Sacado::Fad::DFad<double>* cauchyStressNYX, 
const Sacado::Fad::DFad<double>* cauchyStressNYY, 
const Sacado::Fad::DFad<double>* cauchyStressNYZ, 
const Sacado::Fad::DFad<double>* cauchyStressNZX, 
const Sacado::Fad::DFad<double>* cauchyStressNZY, 
const Sacado::Fad::DFad<double>* cauchyStressNZZ, 
Sacado::Fad::DFad<double>* cauchyStressNP1XX, 
Sacado::Fad::DFad<double>* cauchyStressNP1XY, 
Sacado::Fad::DFad<double>* cauchyStressNP1XZ, 
Sacado::Fad::DFad<double>* cauchyStressNP1YX, 
Sacado::Fad::DFad<double>* cauchyStressNP1YY, 
Sacado::Fad::DFad<double>* cauchyStressNP1YZ, 
Sacado::Fad::DFad<double>* cauchyStressNP1ZX, 
Sacado::Fad::DFad<double>* cauchyStressNP1ZY, 
Sacado::Fad::DFad<double>* cauchyStressNP1ZZ, 
Sacado::Fad::DFad<double>* vonMisesStress,
const Sacado::Fad::DFad<double>* equivalentPlasticStrainN,
Sacado::Fad::DFad<double>* equivalentPlasticStrainNP1,
const int numPoints, 
const double bulkMod,
const double shearMod,
const double yieldStress,
const double strainHardExp,
const double rateHardExp, 
const double refStrainRate,
const double refStrain0,
const double refStrain1,
const bool m_isFlaw,
const double m_flawLocationX,
const double m_flawLocationY,
const double m_flawLocationZ,
const double m_flawSize,
const double m_flawMagnitude,
const double dt
);

}
