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
#include "../materials/correspondence.h"

//! Standard constructor.
PeridigmNS::Compute_Deformation_Gradient::Compute_Deformation_Gradient(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                                       Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                                       Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_),
    m_volumeFId(-1), m_modelCoordinatesFId(-1), m_coordinatesFId(-1),
    m_shapeTensorInverseFId(-1),
    m_deformationGradientFId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_volumeFId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_horizonFId               = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  m_modelCoordinatesFId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_shapeTensorInverseFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Shape_Tensor_Inverse");
  m_deformationGradientFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::FULL_TENSOR, PeridigmField::CONSTANT, "Deformation_Gradient");

  m_fieldIds.push_back(m_volumeFId);
  m_fieldIds.push_back(m_horizonFId);
  m_fieldIds.push_back(m_modelCoordinatesFId);
  m_fieldIds.push_back(m_coordinatesFId);
  m_fieldIds.push_back(m_shapeTensorInverseFId);
  m_fieldIds.push_back(m_deformationGradientFId);
}

//! Destructor.
PeridigmNS::Compute_Deformation_Gradient::~Compute_Deformation_Gradient(){}

//! Fill the deformation gradient tensor
int PeridigmNS::Compute_Deformation_Gradient::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  int retval = 0;

  std::vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* const neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    
    double *volume, *horizon, *modelCoordinates, *coordinates, *shapeTensorInverse, *deformationGradient;
    dataManager->getData(m_volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    dataManager->getData(m_horizonFId, PeridigmField::STEP_NONE)->ExtractView(&horizon);
    dataManager->getData(m_modelCoordinatesFId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
    dataManager->getData(m_coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);
    dataManager->getData(m_shapeTensorInverseFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverse);
    dataManager->getData(m_deformationGradientFId, PeridigmField::STEP_NONE)->ExtractView(&deformationGradient);

    retval = retval || CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                                  horizon,
                                                                                                  modelCoordinates,
                                                                                                  coordinates,
                                                                                                  shapeTensorInverse,
                                                                                                  deformationGradient,
                                                                                                  neighborhoodList,
                                                                                                  numOwnedPoints);
  }

  // Warn if retval not zero
  if (retval && epetraComm->MyPID() == 0)
    std::cout << "**** Warning:  computeApproximateDeformationGradient class returned warning. Some elements may have too few bonds to accurately compute deformation gradient." << std::endl;

  return(retval);
}

