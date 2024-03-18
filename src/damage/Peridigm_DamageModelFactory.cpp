/*! \file Peridigm_DamageModelFactory.hpp */

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

#include <Teuchos_Assert.hpp>
#include "Peridigm_DamageModelFactory.hpp"
#include "Peridigm_CriticalStretchDamageModel.hpp"
#include "Peridigm_UserDefinedTimeDependentCriticalStretchDamageModel.hpp"
#include "Peridigm_InterfaceAwareDamageModel.hpp"
#include "Peridigm_InitialDamageModel.hpp"
#include "Peridigm_JohnsonCookDamageModel.hpp"
#include "Peridigm_VonMisesStressDamageModel.hpp"

Teuchos::RCP<PeridigmNS::DamageModel>
PeridigmNS::DamageModelFactory::create(const Teuchos::ParameterList& damageModelParams)
{
  Teuchos::RCP<PeridigmNS::DamageModel> damageModel;

  const string& damageModelName = damageModelParams.get<string>("Damage Model");

  if(damageModelName == "Critical Stretch")
    damageModel = Teuchos::rcp( new CriticalStretchDamageModel(damageModelParams) );
  else if(damageModelName == "Interface Aware")
    damageModel = Teuchos::rcp( new InterfaceAwareDamageModel(damageModelParams) );
  else if(damageModelName == "Time Dependent Critical Stretch")
    damageModel = Teuchos::rcp( new UserDefinedTimeDependentCriticalStretchDamageModel(damageModelParams) );
  else if(damageModelName == "Initial Damage")
    damageModel = Teuchos::rcp( new InitialDamageModel(damageModelParams) );
  else if(damageModelName == "Johnson Cook")
    damageModel = Teuchos::rcp( new JohnsonCookDamageModel(damageModelParams) );  
  else if(damageModelName == "Von Mises Stress")
    damageModel = Teuchos::rcp( new VonMisesStressDamageModel(damageModelParams) );  
  else {
    string invalidDamageModel("\n**** Unrecognized damage model type: ");
    invalidDamageModel += damageModelName;
    invalidDamageModel += ", must be \"Critical Stretch\", \"Time Dependent Critical Stretch\", \"Interface Aware\", or \"Initial Damage\".\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, invalidDamageModel);
  }

  return damageModel;
}
