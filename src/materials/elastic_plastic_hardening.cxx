//! \file elastic_plastic_hardening.cxx

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
#include <cmath>
#include <Sacado.hpp>
#include <float.h>
#include "elastic_plastic.h"
#include "elastic_plastic_hardening.h"
#include <complex>
#include "Peridigm_Constants.hpp"

namespace MATERIAL_EVALUATION {

template<typename ScalarT>
ScalarT sign(ScalarT v) {
    return v > 0.0 ? 1.0 : (v < 0.0 ? -1.0 : 0);
}

template<typename ScalarT>
void computeInternalForceIsotropicHardeningPlastic
(
    const double* xOverlap,
    const ScalarT* yNP1Overlap,
    const double* mOwned,
    const double* volumeOverlap,
    const ScalarT* dilatationOwned,
    const double* bondDamage,
    const double* scfOwned,
    const double* deviatoricPlasticExtensionStateN,
    ScalarT* deviatoricPlasticExtensionStateNp1,
    const double* lambdaN,
    ScalarT* lambdaNP1,
    ScalarT* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double SHEAR_MODULUS,
    double HORIZON,
    double yieldStress,
    double HARD_MODULUS
)
{

  /*
   * Compute processor local contribution to internal force
   */
  double K = BULK_MODULUS;
  double MU = SHEAR_MODULUS;
  double OMEGA=1.0;
  double DELTA=HORIZON;
  double H=HARD_MODULUS;
  /*
   * 3d variety of yield value 
   */
  double yieldValue = 25.0 * yieldStress * yieldStress / 8 / PeridigmNS::value_of_pi() / pow(DELTA,5);
  /*
   * Planar variety of yield value
   */
//    double THICKNESS=1.0;
//    double yieldValue = 225.0 * yieldStress * yieldStress / 8 / PeridigmNS::value_of_pi() / THICKNESS / pow(DELTA,4);
//    double yieldValue = 0.5 * pow(15*yieldStress/weightedVol,2) * PeridigmNS::value_of_pi() * THICKNESS * pow(DELTA,4) / 16.0;


  const double *xOwned = xOverlap;
  const ScalarT *yOwned = yNP1Overlap;
  const double *m = mOwned;
  const double *v = volumeOverlap;
  const ScalarT *theta = dilatationOwned;
  ScalarT *fOwned = fInternalOverlap;

  const int *neighPtr = localNeighborList;
  double cellVolume, alpha, dx_X, dy_X, dz_X, zeta, edpN;
  ScalarT dx_Y, dy_Y, dz_Y, dY, ed, tdTrial, t, ti, td;
  for(int p=0;p<numOwnedPoints;p++, xOwned +=3, yOwned +=3, fOwned+=3, m++, theta++, lambdaN++, lambdaNP1++, scfOwned++){

    int numNeigh = *neighPtr; neighPtr++;
    const double *X = xOwned;
    const ScalarT *Y = yOwned;
    double weightedVol = *m;
    alpha = *scfOwned * 15.0*MU/weightedVol;
    double selfCellVolume = v[p];
    ScalarT c = 3 * K * (*theta) * OMEGA / weightedVol;
    ScalarT deltaLambda=0.0;

    /*
     * Compute norm of trial stress
     */
    ScalarT tdNorm = 0.0;
    tdNorm = computeDeviatoricForceStateNorm(numNeigh,*theta,neighPtr,bondDamage,deviatoricPlasticExtensionStateN,X,Y,xOverlap,yNP1Overlap,v,alpha,OMEGA);

    /*
     * Evaluate yield function
     */
    double pointWiseYieldValue = *scfOwned * (*scfOwned) * yieldValue;

    /*
     * Compute lambdaNP1 using a backward Euler implicit scheme
     */
    if(tdNorm * tdNorm / 2 - pointWiseYieldValue > 0){
      deltaLambda = updateDeltaLambda(tdNorm, *lambdaN, pointWiseYieldValue, alpha, H);
      if(deltaLambda < 0.0 ){
        deltaLambda = 0.0;
      }
    } else {
      deltaLambda = 0.0;
    }


    //ScalarT f = tdNorm * tdNorm / 2 - pointWiseYieldValue - HARD_MODULUS*(deltaLambda + (*lambdaN));
    //ScalarT f = tdNorm * tdNorm / 2 - pointWiseYieldValue;
    bool elastic = true;

//    std::cout << "Point id = " << p << std::endl;
//    std::cout << "\tyieldStress/m^(4/5) = " << yieldStress/pow(weightedVol,4/5) << std::endl;
//    std::cout << "\tYield Value = " << yieldValue << "; tdNorm * tdNorm / 2 = " << tdNorm * tdNorm / 2 << std::endl;
    if(deltaLambda>0){
      /*
       * This step is incrementally plastic
       */
      //      std::cout << "\t PLASTIC" << std::endl;
      elastic = false;
      //deltaLambda= (*lambdaNP1 - *lambdaN);
      //deltaLambda=( tdNorm / sqrt(2.0*pointWiseYieldValue) - 1.0 ) / alpha;
      *lambdaNP1 = *lambdaN + deltaLambda;
    } else {
//       std::cout << "\t ELASTIC" << std::endl;
      *lambdaNP1 = *lambdaN;
    }

    for(int n=0;n<numNeigh;n++,neighPtr++,bondDamage++, deviatoricPlasticExtensionStateN++, deviatoricPlasticExtensionStateNp1++){
      int localId = *neighPtr;
      cellVolume = v[localId];
      const double *XP = &xOverlap[3*localId];
      const ScalarT *YP = &yNP1Overlap[3*localId];
      dx_X = XP[0]-X[0];
      dy_X = XP[1]-X[1];
      dz_X = XP[2]-X[2];
      zeta = sqrt(dx_X*dx_X+dy_X*dy_X+dz_X*dz_X);
      dx_Y = YP[0]-Y[0];
      dy_Y = YP[1]-Y[1];
      dz_Y = YP[2]-Y[2];
      dY = sqrt(dx_Y*dx_Y+dy_Y*dy_Y+dz_Y*dz_Y);
      /*
       * Deviatoric extension state
       */
      ed = dY-zeta-*theta*zeta/3;

      /*
       * Deviatoric plastic extension state from last step
       */
      edpN = *deviatoricPlasticExtensionStateN;

      /*
       * Compute trial stress
       */
      tdTrial = alpha * OMEGA * (ed - edpN);

      /*
       * Evaluate yield function
       */
      if(elastic){
        /*
         * Elastic case
         */
        td = tdTrial;

        /*
         * Therefore edpNp1 = edpN
         */
        *deviatoricPlasticExtensionStateNp1 = *deviatoricPlasticExtensionStateN;

      } else {
        /*
         * Compute deviatoric force state
         */
        //td = sqrt(2.0*pointWiseYieldValue) * tdTrial / tdNorm;
                td = tdTrial / (1+alpha*deltaLambda);

        /*
         * Update deviatoric plastic deformation state
         */
        *deviatoricPlasticExtensionStateNp1 = edpN + td * deltaLambda;

//        std::cout << "Neighbor Id = " << localId << "; Updating deviatoricPlasticExtensionState = " << *deviatoricPlasticExtensionState << std::endl;
      }
//      std::cout << "\tNeighbor Id = " << localId << "\n\ttd = " << td;
      /*
       * Compute isotropic part of force state
       */
      ti = c * zeta;

      /*
       * Force state (with damage)
       */
      double d=(1.0-*bondDamage);
      t = d*(ti + d*td);

      /*
       * Assemble pair wise force function
       */
      ScalarT fx = t * dx_Y / dY;
      ScalarT fy = t * dy_Y / dY;
      ScalarT fz = t * dz_Y / dY;

      *(fOwned+0) += fx*cellVolume;
      *(fOwned+1) += fy*cellVolume;
      *(fOwned+2) += fz*cellVolume;
      fInternalOverlap[3*localId+0] -= fx*selfCellVolume;
      fInternalOverlap[3*localId+1] -= fy*selfCellVolume;
      fInternalOverlap[3*localId+2] -= fz*selfCellVolume;
    }

  }
}

template<typename ScalarT>
ScalarT updateDeltaLambda
(
    ScalarT tdNorm,
    const double lambdaN,
    double pointWiseYieldValue,
    double alpha,
    double HARD_MODULUS
)
{   
    double  a  = 2./alpha + lambdaN + pointWiseYieldValue/HARD_MODULUS;
    double  b = (HARD_MODULUS + 2.*alpha*HARD_MODULUS*lambdaN + 2.*alpha*pointWiseYieldValue)/(alpha*alpha*HARD_MODULUS);
    ScalarT c = (2.*HARD_MODULUS*lambdaN + 2.*pointWiseYieldValue - tdNorm*tdNorm)/(2.*alpha*alpha*HARD_MODULUS);
    
    double  Q = (a*a - 3.*b)/9.;
    ScalarT R = (2.*a*a*a - 9.*a*b + 27.*c)/54.;
    
    //std::cout << "a1 = " << 2./alpha << " a3 = " << pointWiseYieldValue/HARD_MODULUS << std::endl;

    if( R*R < Q*Q*Q )
    {
        ScalarT theta = acos( R/sqrt(Q*Q*Q));

        ScalarT root1 = -2.*sqrt(Q)*cos(theta/3.) - a/3.;
        ScalarT root2 = -2.*sqrt(Q)*cos(theta/3. + 2.*PeridigmNS::value_of_pi()/3.) - a/3.;
        ScalarT root3 = -2.*sqrt(Q)*cos(theta/3. - 2.*PeridigmNS::value_of_pi()/3.) - a/3.;
        
        //if (root1 > 0 || root2 > 0 || root3 > 0)
            //std::cout << "sol 1 " << root1  << " sol 2 " << root2  << " sol 3 " << root3 << std::endl;

        ScalarT min = DBL_MAX;

        if( root1 > 0.0 && std::abs(root1) < std::abs(min) )
            min = root1;
        if( root2 > 0.0 && std::abs(root2) < std::abs(min) )
            min = root2;
        if( root3 > 0.0 && std::abs(root3) < std::abs(min) )
            min = root3;
        if( min == DBL_MAX ){
            //std::cout << "MAX_DOUBLE being returned for deltaLambda" << std::endl;
            //exit(1);
            min = 0.0;
        }

        //if (root1 > 0 || root2 > 0 || root3 > 0)
            //std::cout << "ans " << min << std::endl << std::endl; 
        
        return min;
    
    }
    else
    {
        ScalarT A = -sign(R)*pow( std::abs(R) + sqrt(R*R - Q*Q*Q) ,1./3.);
        if( std::abs(A) < 1e-50)
            std::cout << "A is VERY small! A = " << A << std::endl;
        
        ScalarT B = Q/A;
        

        //std::cout << "ans " << (A + B) - a/3.<< std::endl << std::endl; 
        ScalarT val = ( A + B ) - a/3.;
        if(val < 0.)
            return 0.;
        else
            //std::cout << "ans2 " << val << std::endl;
            return val;
    }

}

/** Explicit template instantiation for double. */
template void computeInternalForceIsotropicHardeningPlastic<double>
(
    const double* xOverlap,
    const double* yNP1Overlap,
    const double* mOwned,
    const double* volumeOverlap,
    const double* dilatationOwned,
    const double* bondDamage,
    const double* scfOwned,
    const double* deviatoricPlasticExtensionStateN,
    double* deviatoricPlasticExtensionStateNp1,
    const double* lambdaN,
    double* lambdaNP1,
    double* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double SHEAR_MODULUS,
    double HORIZON,
    double yieldStress,
    double HARD_MODULUS
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template void computeInternalForceIsotropicHardeningPlastic<Sacado::Fad::DFad<double> >
(
    const double* xOverlap,
    const Sacado::Fad::DFad<double>* yNP1Overlap,
    const double* mOwned,
    const double* volumeOverlap,
    const Sacado::Fad::DFad<double>* dilatationOwned,
    const double* bondDamage,
    const double* scfOwned,
    const double* deviatoricPlasticExtensionStateN,
    Sacado::Fad::DFad<double>* deviatoricPlasticExtensionStateNp1,
    const double* lambdaN,
    Sacado::Fad::DFad<double>* lambdaNP1,
    Sacado::Fad::DFad<double>* fInternalOverlap,
    const int*  localNeighborList,
    int numOwnedPoints,
    double BULK_MODULUS,
    double SHEAR_MODULUS,
    double HORIZON,
    double yieldStress,
    double HARD_MODULUS
);

/** Explicit template instantiation for int. */
template double sign<double> 
(
    double v
);

/** Explicit template instantiation for Sacado::Fad::DFad<int>. */
template Sacado::Fad::DFad<double> sign<Sacado::Fad::DFad<double> >
(
    Sacado::Fad::DFad<double> v
);


template double updateDeltaLambda<double>
(
    double tdNorm,
    const double lambdaN,
    double pointWiseYieldValue,
    double alpha,
    double HARD_MODULUS
);

/** Explicit template instantiation for Sacado::Fad::DFad<double>. */
template Sacado::Fad::DFad<double> updateDeltaLambda<Sacado::Fad::DFad<double> >
(
    Sacado::Fad::DFad<double> tdNorm,
    const double lambdaN,
    double pointWiseYieldValue,
    double alpha,
    double HARD_MODULUS
);

}
