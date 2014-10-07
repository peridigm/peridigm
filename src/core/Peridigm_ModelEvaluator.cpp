/*! \file Peridigm_ModelEvaluator.hpp */

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

#include "Peridigm_ModelEvaluator.hpp"

using namespace std;

PeridigmNS::ModelEvaluator::ModelEvaluator(){
}

PeridigmNS::ModelEvaluator::~ModelEvaluator(){
}

void 
PeridigmNS::ModelEvaluator::evalModel(Teuchos::RCP<Workset> workset) const
{
  const double dt = workset->timeStep;
  std::vector<PeridigmNS::Block>::iterator blockIt;

  // ---- Evaluate Damage ---

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<const PeridigmNS::DamageModel> damageModel = blockIt->getDamageModel();
    if(!damageModel.is_null()){
      Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
      const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
      const int* ownedIDs = neighborhoodData->OwnedIDs();
      const int* neighborhoodList = neighborhoodData->NeighborhoodList();
      Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
      damageModel->computeDamage(dt, 
                                 numOwnedPoints,
                                 ownedIDs,
                                 neighborhoodList,
                                 *dataManager);
    }
  }

  // ---- Evaluate Internal Force ----

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

    materialModel->computeForce(dt, 
                                numOwnedPoints,
                                ownedIDs,
                                neighborhoodList,
                                *dataManager);
  }

  // ---- Evaluate Contact ----
  if(!workset->contactManager.is_null())
    workset->contactManager->evaluateContactForce(dt);
}

void 
PeridigmNS::ModelEvaluator::evalJacobian(Teuchos::RCP<Workset> workset) const
{
  const double dt = workset->timeStep;
  std::vector<PeridigmNS::Block>::iterator blockIt;
  PeridigmNS::Material::JacobianType jacobianType = *(workset->jacobianType);
  PeridigmNS::SerialMatrix& jacobian = *(workset->jacobian);

  // ---- Compute the Tangent Stiffness Matrix ----

  for(blockIt = workset->blocks->begin() ; blockIt != workset->blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* ownedIDs = neighborhoodData->OwnedIDs();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::Material> materialModel = blockIt->getMaterialModel();

    materialModel->computeJacobian(dt, 
                                   numOwnedPoints,
                                   ownedIDs,
                                   neighborhoodList,
                                   *dataManager,
                                   jacobian,
                                   jacobianType);
  }
}
