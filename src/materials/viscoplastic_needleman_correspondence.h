//! \file elastic_plastic_correspondence.h

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
#ifndef VISCO_PLASTIC_NEEDLEMAN_CORRESPONDENCE_H
#define VISCO_PLASTIC_NEEDLEMAN_CORRESPONDENCE_H

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
    const bool m_isFlaw,
    const double m_flawLocationX,
    const double m_flawLocationY,
    const double m_flawLocationZ,
    const double m_flawSize,
    const double m_flawMagnitude,
    const double dt
);

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
);

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
);

}

#endif // VISCO_PLASTIC_NEEDLEMAN_CORRESPONDENCE_H
