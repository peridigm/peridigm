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
    for(int j = 0; j < 3; j++){
      for(int i = 0; i < 3; i++){
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
      // function is rate dependent and can change over a load step.


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

      for (int i = 0; i < 9; i++) {
        deviatoricStressN[i] = *(stressN+i);
      }
      deviatoricStressN[0] -= sphericalStressN;
      deviatoricStressN[4] -= sphericalStressN;
      deviatoricStressN[8] -= sphericalStressN;

      // Compute \S_ij * \S_ij
      tempScalar = 0.0;
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          tempScalar += deviatoricStressN[i+3*j] * deviatoricStressN[i+3*j];
        }
      }

      deviatoricStressMagnitudeN = sqrt(tempScalar);

      //Verbose bi-section method to solve for deltaLambda
      //TODO: Implement a ``smarter'' faster root solve here, i.e. Brent's
      //method
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
      *eqpsNP1 = *eqpsN + sqrt(2.0/3.0) * deltaLambda;
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
  ScalarT hardTerm = yieldStress * pow(1.0 + eqps/refStrain0, strainHardExp) / (1.0 + pow(eqps/refStrain1,2.0));
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

  for( ; ; )
  {
    if( fabs ( fc ) < fabs ( fb ) )
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

    if( fabs ( m ) <= tol || fb == 0.0 )
    {
      break;
    }

    if( fabs ( e ) < tol || fabs ( fa ) <= fabs ( fb ) )
    {
      e = m;
      d = e;
    }
    else
    {
      s = fb / fa;

      if( sa == c )
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

      if( 0.0 < p )
      {
        q = - q;
      }
      else
      {
        p = - p;
      }

      s = e;
      e = d;

      if( 2.0 * p < 3.0 * m * q - fabs ( tol * q ) &&
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

    if( tol < fabs ( d ) )
    {
      sb = sb + d;
    }
    else if( 0.0 < m )
    {
      sb = sb + tol;
    }
    else
    {
      sb = sb - tol;
    }

    yfb = ViscoplasticNeedlemanYieldFunction( sb, eqps, yieldStress, strainHardExp, rateHardExp, refStrainRate,refStrain0, refStrain1, dt);
    fb = scalarDeviatoricStrainInc - sb - 1.0/2.0/shearMod*( sqrt(2.0/3.0)*yfb - deviatoricStressMagnitude);

    if( ( 0.0 < fb && 0.0 < fc ) || ( fb <= 0.0 && fc <= 0.0 ) )
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
