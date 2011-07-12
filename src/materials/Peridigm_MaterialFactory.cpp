/*! \file Peridigm_MaterialFactory.hpp */

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

#include <Teuchos_TestForException.hpp>
#include "Peridigm_MaterialFactory.hpp"
#include "Peridigm_LinearElasticIsotropicMaterial.hpp"
#include "Peridigm_IsotropicElasticPlasticMaterial.hpp"
#include "Peridigm_ViscoelasticStandardLinearSolid.hpp"

Teuchos::RCP<PeridigmNS::Material>
PeridigmNS::MaterialFactory::create(Teuchos::RCP<const Teuchos::ParameterList>& materialParams)
{
  Teuchos::RCP<PeridigmNS::Material> materialModel;

  for(Teuchos::ParameterList::ConstIterator it = materialParams->begin() ; it != materialParams->end() ; it++){
    const string& name = it->first;
    const Teuchos::ParameterList& params = materialParams->sublist(name);
    if (name == "Linear Elastic")
      materialModel = Teuchos::rcp( new LinearElasticIsotropicMaterial(params) );
    else if (name == "Elastic Plastic")
      materialModel = Teuchos::rcp( new IsotropicElasticPlasticMaterial(params) );
    else if (name == "Viscoelastic Standard Linear Solid")
      materialModel = Teuchos::rcp( new ViscoelasticStandardLinearSolid(params) );
    else {
      string invalidMaterial("\n**** Unrecognized material model: ");
      invalidMaterial += name;
      invalidMaterial += ", must be Linear Elastic or Elastic Plastic or Viscoelastic Standard Linear Solid.\n";
      TEST_FOR_EXCEPT_MSG(true, invalidMaterial);
    }
  }
  
  return materialModel;
}

//Teuchos::rcp_implicit_cast<Material>(material)
