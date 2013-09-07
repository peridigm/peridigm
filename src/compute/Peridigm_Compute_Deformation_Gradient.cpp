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

// Initialize ID generator
int PeridigmNS::Compute_Deformation_Gradient::myIDGenerator = 0;

//! Standard constructor.
PeridigmNS::Compute_Deformation_Gradient::Compute_Deformation_Gradient(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                                       Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                                       Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_),
    m_volumeFId(-1), m_modelCoordinatesFId(-1), m_coordinatesFId(-1),
    m_shapeTensorInverseXXFId(-1), m_shapeTensorInverseXYFId(-1), m_shapeTensorInverseXZFId(-1),
    m_shapeTensorInverseYXFId(-1), m_shapeTensorInverseYYFId(-1), m_shapeTensorInverseYZFId(-1),
    m_shapeTensorInverseZXFId(-1), m_shapeTensorInverseZYFId(-1), m_shapeTensorInverseZZFId(-1),
    m_deformationGradientXXFId(-1), m_deformationGradientXYFId(-1), m_deformationGradientXZFId(-1),
    m_deformationGradientYXFId(-1), m_deformationGradientYYFId(-1), m_deformationGradientYZFId(-1),
    m_deformationGradientZXFId(-1), m_deformationGradientZYFId(-1), m_deformationGradientZZFId(-1)
{
  // Initialize my unique ID
  myID = myIDGenerator++;

  FieldManager& fieldManager = FieldManager::self();
  m_volumeFId                = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  m_modelCoordinatesFId      = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_coordinatesFId           = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  m_shapeTensorInverseXXFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXX");
  m_shapeTensorInverseXYFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXY");
  m_shapeTensorInverseXZFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseXZ");
  m_shapeTensorInverseYXFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYX");
  m_shapeTensorInverseYYFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYY");
  m_shapeTensorInverseYZFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseYZ");
  m_shapeTensorInverseZXFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZX");
  m_shapeTensorInverseZYFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZY");
  m_shapeTensorInverseZZFId  = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Shape_Tensor_InverseZZ");
  m_deformationGradientXXFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientXX");
  m_deformationGradientXYFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientXY");
  m_deformationGradientXZFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientXZ");
  m_deformationGradientYXFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientYX");
  m_deformationGradientYYFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientYY");
  m_deformationGradientYZFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientYZ");
  m_deformationGradientZXFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientZX");
  m_deformationGradientZYFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientZY");
  m_deformationGradientZZFId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Deformation_GradientZZ");

  m_fieldIds.push_back(m_volumeFId);
  m_fieldIds.push_back(m_modelCoordinatesFId);
  m_fieldIds.push_back(m_coordinatesFId);
  m_fieldIds.push_back(m_shapeTensorInverseXXFId);
  m_fieldIds.push_back(m_shapeTensorInverseXYFId);
  m_fieldIds.push_back(m_shapeTensorInverseXZFId);
  m_fieldIds.push_back(m_shapeTensorInverseYXFId);
  m_fieldIds.push_back(m_shapeTensorInverseYYFId);
  m_fieldIds.push_back(m_shapeTensorInverseYZFId);
  m_fieldIds.push_back(m_shapeTensorInverseZXFId);
  m_fieldIds.push_back(m_shapeTensorInverseZYFId);
  m_fieldIds.push_back(m_shapeTensorInverseZZFId);
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

  std::vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    int* const neighborhoodList = neighborhoodData->NeighborhoodList();
    double horizon = blockIt->getMaterialModel()->Horizon();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = blockIt->getDataManager();
    
    double *volume, *modelCoordinates, *coordinates;
    dataManager->getData(m_volumeFId, PeridigmField::STEP_NONE)->ExtractView(&volume);
    dataManager->getData(m_modelCoordinatesFId, PeridigmField::STEP_NONE)->ExtractView(&modelCoordinates);
    dataManager->getData(m_coordinatesFId, PeridigmField::STEP_NP1)->ExtractView(&coordinates);

    double *shapeTensorInverseXX, *shapeTensorInverseXY, *shapeTensorInverseXZ;
    double *shapeTensorInverseYX, *shapeTensorInverseYY, *shapeTensorInverseYZ;
    double *shapeTensorInverseZX, *shapeTensorInverseZY, *shapeTensorInverseZZ;
    dataManager->getData(m_shapeTensorInverseXXFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseXX);
    dataManager->getData(m_shapeTensorInverseXYFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseXY);
    dataManager->getData(m_shapeTensorInverseXZFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseXZ);
    dataManager->getData(m_shapeTensorInverseYXFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseYX);
    dataManager->getData(m_shapeTensorInverseYYFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseYY);
    dataManager->getData(m_shapeTensorInverseYZFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseYZ);
    dataManager->getData(m_shapeTensorInverseZXFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseZX);
    dataManager->getData(m_shapeTensorInverseZYFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseZY);
    dataManager->getData(m_shapeTensorInverseZZFId, PeridigmField::STEP_NONE)->ExtractView(&shapeTensorInverseZZ);

    double *deformationGradientXX, *deformationGradientXY, *deformationGradientXZ;
    double *deformationGradientYX, *deformationGradientYY, *deformationGradientYZ;
    double *deformationGradientZX, *deformationGradientZY, *deformationGradientZZ;
    dataManager->getData(m_deformationGradientXXFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientXX);
    dataManager->getData(m_deformationGradientXYFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientXY);
    dataManager->getData(m_deformationGradientXZFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientXZ);
    dataManager->getData(m_deformationGradientYXFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientYX);
    dataManager->getData(m_deformationGradientYYFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientYY);
    dataManager->getData(m_deformationGradientYZFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientYZ);
    dataManager->getData(m_deformationGradientZXFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientZX);
    dataManager->getData(m_deformationGradientZYFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientZY);
    dataManager->getData(m_deformationGradientZZFId, PeridigmField::STEP_NP1)->ExtractView(&deformationGradientZZ);
    
    retval = retval || CORRESPONDENCE::computeShapeTensorInverseAndApproximateDeformationGradient(volume,
                                                                                                  modelCoordinates,
                                                                                                  coordinates,
                                                                                                  shapeTensorInverseXX,
                                                                                                  shapeTensorInverseXY,
                                                                                                  shapeTensorInverseXZ,
                                                                                                  shapeTensorInverseYX,
                                                                                                  shapeTensorInverseYY,
                                                                                                  shapeTensorInverseYZ,
                                                                                                  shapeTensorInverseZX,
                                                                                                  shapeTensorInverseZY,
                                                                                                  shapeTensorInverseZZ,
                                                                                                  deformationGradientXX,
                                                                                                  deformationGradientXY,
                                                                                                  deformationGradientXZ,
                                                                                                  deformationGradientYX,
                                                                                                  deformationGradientYY,
                                                                                                  deformationGradientYZ,
                                                                                                  deformationGradientZX,
                                                                                                  deformationGradientZY,
                                                                                                  deformationGradientZZ,
                                                                                                  neighborhoodList,
                                                                                                  numOwnedPoints,
                                                                                                  horizon);
  }

  // Warn if retval not zero
  if (retval && epetraComm->MyPID() == 0)
    std::cout << "**** Warning:  computeApproximateDeformationGradient class returned warning. Some elements may have too few bonds to accurately compute deformation gradient." << std::endl;

  return(retval);
}

