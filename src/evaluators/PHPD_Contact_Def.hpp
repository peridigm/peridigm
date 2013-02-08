/*! \file PHPD_Contact_Def.hpp */

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
#include <Phalanx_DataLayout.hpp>

//**********************************************************************
template<typename EvalT, typename Traits>
Contact<EvalT, Traits>::Contact(Teuchos::ParameterList& p) :
  m_verbose(false),
  m_num_pt(0)
{
  if(p.isParameter("Verbose"))
	 m_verbose = p.get<bool>("Verbose");

  contact_field_tag = 
    Teuchos::rcp(new PHX::Tag<ScalarT>("Contact",p.get< Teuchos::RCP<PHX::DataLayout> >("Dummy Data Layout")));

  this->addEvaluatedField(*contact_field_tag);

  this->setName("Contact");
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Contact<EvalT, Traits>::postRegistrationSetup(
                      typename Traits::SetupData d,
                      PHX::FieldManager<Traits>& fm)
{
}

//**********************************************************************
template<typename EvalT, typename Traits>
void Contact<EvalT, Traits>::evaluateFields(typename Traits::EvalData cellData)
{
  const double dt = *cellData.timeStep;

  std::vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = cellData.blocks->begin() ; blockIt != cellData.blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> contactNeighborhoodData = blockIt->getContactNeighborhoodData();
    const int numOwnedPoints = contactNeighborhoodData->NumOwnedPoints();
    const int* ownedIDs = contactNeighborhoodData->OwnedIDs();
    const int* neighborhoodList = contactNeighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::ContactModel> contactModel = blockIt->getContactModel();

    if(!contactModel.is_null())
      contactModel->computeForce(dt, 
                                 numOwnedPoints,
                                 ownedIDs,
                                 neighborhoodList,
                                 *dataManager);
  }
}

