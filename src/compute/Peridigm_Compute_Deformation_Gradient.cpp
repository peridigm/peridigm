/*! \file Peridigm_Compute_Deformation_Gradient.cpp */

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

#include "Peridigm_Compute_Deformation_Gradient.hpp"
#include "Peridigm_Field.hpp"

// Initialize ID generator
int PeridigmNS::Compute_Deformation_Gradient::myIDGenerator = 0;

//! Standard constructor.
PeridigmNS::Compute_Deformation_Gradient::Compute_Deformation_Gradient(Teuchos::RCP<const Teuchos::ParameterList> params,
                                         Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(params, epetraComm_), m_deformationGradientXXFId(-1), m_deformationGradientXYFId(-1), m_deformationGradientXZFId(-1),
                                  m_deformationGradientYXFId(-1), m_deformationGradientYYFId(-1), m_deformationGradientYZFId(-1),
                                  m_deformationGradientZXFId(-1), m_deformationGradientZYFId(-1), m_deformationGradientZZFId(-1)

{
  // Initialize my unique ID
  myID = myIDGenerator++;

  FieldManager& fieldManager = FieldManager::self();

  m_deformationGradientXXFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXX");
  m_deformationGradientXYFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXY");
  m_deformationGradientXZFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientXZ");
  m_deformationGradientYXFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYX");
  m_deformationGradientYYFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYY");
  m_deformationGradientYZFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientYZ");
  m_deformationGradientZXFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZX");
  m_deformationGradientZYFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZY");
  m_deformationGradientZZFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Deformation_GradientZZ");

  m_fieldIds.push_back(m_deformationGradientXXFId);
  m_fieldIds.push_back(m_deformationGradientXYFId);
  m_fieldIds.push_back(m_deformationGradientXZFId);
  m_fieldIds.push_back(m_deformationGradientYXFId);
  m_fieldIds.push_back(m_deformationGradientYYFId);
  m_fieldIds.push_back(m_deformationGradientYZFId);
  m_fieldIds.push_back(m_deformationGradientZXFId);
  m_fieldIds.push_back(m_deformationGradientZYFId);
  m_fieldIds.push_back(m_deformationGradientZZFId);
}

//! Destructor.
PeridigmNS::Compute_Deformation_Gradient::~Compute_Deformation_Gradient(){}

//! Fill the deformation gradient tensor
int PeridigmNS::Compute_Deformation_Gradient::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  int retval = 0;

  // Let object with ID 0 do all the work
  if (!myID) return(retval);

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
// DJL: Note to Mike, computeApproximateDeformationGradient() returns null, so this causes a compile error on latrobe
//     retval = retval || materialModel->computeApproximateDeformationGradient(numOwnedPoints, ownedIDs, neighborhoodList, *dataManager);
    materialModel->computeApproximateDeformationGradient(numOwnedPoints, ownedIDs, neighborhoodList, *dataManager);
  }

  // Warn if retval not zero
  if (retval && epetraComm->MyPID() == 0)
    cout << "**** Warning:  computeApproximateDeformationGradient class returned warning. Some elements may have too few bonds to accurately compute deformation gradient." << endl;

  return(retval);

}

