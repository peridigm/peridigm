//! \file elastic_correspondence.h

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
#ifndef ELASTIC_CORRESPONDENCE_H
#define ELASTIC_CORRESPONDENCE_H

namespace CORRESPONDENCE {

template<typename ScalarT>
void updateElasticCauchyStress
(
const ScalarT* unrotatedRateOfDeformationXX, 
const ScalarT* unrotatedRateOfDeformationXY, 
const ScalarT* unrotatedRateOfDeformationXZ, 
const ScalarT* unrotatedRateOfDeformationYX, 
const ScalarT* unrotatedRateOfDeformationYY, 
const ScalarT* unrotatedRateOfDeformationYZ, 
const ScalarT* unrotatedRateOfDeformationZX, 
const ScalarT* unrotatedRateOfDeformationZY, 
const ScalarT* unrotatedRateOfDeformationZZ, 
const ScalarT* unrotatedCauchyStressNXX, 
const ScalarT* unrotatedCauchyStressNXY, 
const ScalarT* unrotatedCauchyStressNXZ, 
const ScalarT* unrotatedCauchyStressNYX, 
const ScalarT* unrotatedCauchyStressNYY, 
const ScalarT* unrotatedCauchyStressNYZ, 
const ScalarT* unrotatedCauchyStressNZX, 
const ScalarT* unrotatedCauchyStressNZY, 
const ScalarT* unrotatedCauchyStressNZZ, 
ScalarT* unrotatedCauchyStressNP1XX, 
ScalarT* unrotatedCauchyStressNP1XY, 
ScalarT* unrotatedCauchyStressNP1XZ, 
ScalarT* unrotatedCauchyStressNP1YX, 
ScalarT* unrotatedCauchyStressNP1YY, 
ScalarT* unrotatedCauchyStressNP1YZ, 
ScalarT* unrotatedCauchyStressNP1ZX, 
ScalarT* unrotatedCauchyStressNP1ZY, 
ScalarT* unrotatedCauchyyStressNP1ZZ, 
const int numPoints, 
const double bulkMod,
const double shearMod,
const double dt
);

}

#endif // ELASTIC_CORRESPONDENCE_H
