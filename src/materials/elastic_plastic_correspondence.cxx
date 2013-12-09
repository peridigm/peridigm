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
      for (int i = 0; i < 9; i++) {
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
      for (int i = 0; i < 9; i++) {
          *(stressNP1+i) = *(stressN+i) + deviatoricStrainInc[i]*2.0*shearMod;
      }
      *(stressNP1) += bulkMod*dilatationInc;
      *(stressNP1+4) += bulkMod*dilatationInc;
      *(stressNP1+8) += bulkMod*dilatationInc;

      sphericalStressNP1 = (*(stressNP1) + *(stressNP1+4) + *(stressNP1+8))/3.0;

      // Compute the ``trial'' von Mises stress
      for (int i = 0; i < 9; i++) {
          deviatoricStressNP1[i] = *(stressNP1+i);
      }
      deviatoricStressNP1[0] -= sphericalStressNP1;
      deviatoricStressNP1[4] -= sphericalStressNP1;
      deviatoricStressNP1[8] -= sphericalStressNP1;

      // Compute \sigma_ij * \sigma_ij
      tempScalar = 0.0;
      for (int j = 0; j < 3; j++) {
          for (int i = 0; i < 3; i++) {
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
