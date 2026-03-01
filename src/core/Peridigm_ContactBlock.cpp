/*! \file Peridigm_ContactBlock.cpp */

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

#include "Peridigm_ContactBlock.hpp"
#include "Peridigm_Field.hpp"
#include <vector>
#include <set>

using std::vector;


void PeridigmNS::ContactBlock::initialize(Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOverlapScalarPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOwnedVectorPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOverlapVectorPointMap,
                                   Teuchos::RCP<const Epetra_BlockMap> globalOwnedScalarBondMap,
                                   Teuchos::RCP<const Epetra_Vector>   globalBlockIds,
                                   Teuchos::RCP<const PeridigmNS::NeighborhoodData> globalNeighborhoodData)
{
  // Initialize the base class to create maps and neighborhood data
  BlockBase::initialize(globalOwnedScalarPointMap,
                    globalOverlapScalarPointMap,
                    globalOwnedVectorPointMap,
                    globalOverlapVectorPointMap,
                    globalOwnedScalarBondMap,
                    globalBlockIds,
                    globalNeighborhoodData);

  // Initialize the data manager
  TEUCHOS_TEST_FOR_TERMINATION(contactModel.is_null(),
                              "\n**** Contact model must be set via ContactBlock::setContactModel() prior to calling ContactBlock::initialize()\n");

  // Collect all the required field Ids
  vector<int> fieldIds;
  // Ids passed in via setAuxiliaryFieldIds(), if any
  fieldIds.insert(fieldIds.end(), auxiliaryFieldIds.begin(), auxiliaryFieldIds.end());
  // Contact model fieldIds
  vector<int> contactModelFieldIds = contactModel->FieldIds();
  fieldIds.insert(fieldIds.end(), contactModelFieldIds.begin(), contactModelFieldIds.end());

  BlockBase::initializeDataManager(fieldIds);
}

void PeridigmNS::ContactBlock::rebalance(Teuchos::RCP<const Epetra_BlockMap> rebalancedGlobalOwnedScalarPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> rebalancedGlobalOverlapScalarPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> rebalancedGlobalOwnedVectorPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> rebalancedGlobalOverlapVectorPointMap,
                                  Teuchos::RCP<const Epetra_BlockMap> rebalancedGlobalOwnedScalarBondMap,
                                  Teuchos::RCP<const Epetra_Vector> rebalancedGlobalBlockIds,
                                  Teuchos::RCP<const PeridigmNS::NeighborhoodData> rebalancedGlobalNeighborhoodData)
{
  createMapsFromGlobalMaps(rebalancedGlobalOwnedScalarPointMap,
                           rebalancedGlobalOverlapScalarPointMap,
                           rebalancedGlobalOwnedVectorPointMap,
                           rebalancedGlobalOverlapVectorPointMap,
                           rebalancedGlobalOwnedScalarBondMap,
                           rebalancedGlobalBlockIds,
                           rebalancedGlobalNeighborhoodData);

  neighborhoodData = createNeighborhoodDataFromGlobalNeighborhoodData(rebalancedGlobalOverlapScalarPointMap,
                                                                      rebalancedGlobalNeighborhoodData);

  dataManager->rebalance(ownedScalarPointMap,
                         overlapScalarPointMap,
                         ownedVectorPointMap,
                         overlapVectorPointMap,
                         ownedScalarBondMap);
}
