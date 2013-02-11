/*! \file Peridigm_Compute_Strain_Energy.cpp */

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

#include <vector>

#include "Peridigm_Compute_Strain_Energy.hpp"
#include "Peridigm_Field.hpp"

//! Standard constructor.
PeridigmNS::Compute_Strain_Energy::Compute_Strain_Energy(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                         Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(params, epetraComm_), m_strainEnergyFId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_strainEnergyFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Strain_Energy");
  m_fieldIds.push_back(m_strainEnergyFId);
}

//! Destructor.
PeridigmNS::Compute_Strain_Energy::~Compute_Strain_Energy(){}

//! Compute the strain energy
int PeridigmNS::Compute_Strain_Energy::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  int retval = 0;

  Teuchos::RCP<const PeridigmNS::Material> materialModel;
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData;
  Teuchos::RCP<PeridigmNS::DataManager> dataManager;
  std::vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    materialModel = blockIt->getMaterialModel();
    neighborhoodData = blockIt->getNeighborhoodData();
    dataManager = blockIt->getDataManager();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* const ownedIDs = neighborhoodData->OwnedIDs();
    int* const neighborhoodList = neighborhoodData->NeighborhoodList();
    materialModel->computeStrainEnergy(0.0, numOwnedPoints, ownedIDs, neighborhoodList, *dataManager);
  }

  return(retval);
}

