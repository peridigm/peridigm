/*! \file Peridigm_Compute_Nearest_Point_Data.cpp */

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

#include "Peridigm_Compute_Nearest_Point_Data.hpp"
#include "Peridigm_Field.hpp"

using namespace std;

//! Standard constructor.
PeridigmNS::Compute_Nearest_Point_Data::Compute_Nearest_Point_Data(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                                   Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(params, epetraComm_), m_elementId(-1.0), m_blockId(-1), m_verbose(false), m_elementIdFieldId(-1),
    m_modelCoordinatesFieldId(-1), m_variableFieldId(-1), m_outputFieldId(-1), m_outputXFieldId(-1),
    m_outputYFieldId(-1), m_outputZFieldId(-1)
{
  m_positionX = params->get<double>("X");
  m_positionY = params->get<double>("Y");
  m_positionZ = params->get<double>("Z");
  m_variable = params->get<string>("Variable");
  m_outputLabel = params->get<string>("Output Label");

  if(params->isParameter("Verbose"))
    m_verbose = params->get<bool>("Verbose");

  FieldManager& fieldManager = FieldManager::self();
  m_elementIdFieldId = fieldManager.getFieldId("Element_Id");
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_variableFieldId = fieldManager.getFieldId(m_variable);
  m_fieldIds.push_back(m_elementIdFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_variableFieldId);

  PeridigmField::Length length = fieldManager.getFieldSpec(m_variableFieldId).getLength();
  if(length == PeridigmField::SCALAR)
    m_variableLength = 1;
  else if(length == PeridigmField::VECTOR)
    m_variableLength = 3;
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  Nearest_Point_Data compute class can be called only for SCALAR or VECTOR data.\n");

  if(m_variableLength == 1){
    m_outputFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, m_outputLabel);
    m_fieldIds.push_back(m_outputFieldId);
  }
  else if(m_variableLength == 3){
    m_outputXFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, m_outputLabel+"_X");
    m_outputYFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, m_outputLabel+"_Y");
    m_outputZFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, m_outputLabel+"_Z");
    m_fieldIds.push_back(m_outputXFieldId);
    m_fieldIds.push_back(m_outputYFieldId);
    m_fieldIds.push_back(m_outputZFieldId);
  }

  PeridigmField::Temporal temporal = fieldManager.getFieldSpec(m_variableFieldId).getTemporal();
  m_variableIsStated = false;
  if(temporal == PeridigmField::TWO_STEP)
    m_variableIsStated = true;
}

//! Destructor.
PeridigmNS::Compute_Nearest_Point_Data::~Compute_Nearest_Point_Data(){}

void PeridigmNS::Compute_Nearest_Point_Data::initialize( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) {

  // Find the node that is closest to the specified location
  double minDistanceSquared = DBL_MAX;
  m_elementId = INT_MAX;

  // Keep track of any ties that are found
  bool foundTies(false);

  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    double *x, *id;
    blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
    blockIt->getData(m_elementIdFieldId, PeridigmField::STEP_NONE)->ExtractView(&id);

    double distanceSquared;
    int numOwnedPoints = blockIt->getNeighborhoodData()->NumOwnedPoints();

    for(int iID=0 ; iID<numOwnedPoints ; ++iID){
      distanceSquared = (x[3*iID] - m_positionX)*(x[3*iID] - m_positionX) 
        + (x[3*iID+1] - m_positionY)*(x[3*iID+1] - m_positionY) 
        + (x[3*iID+2] - m_positionZ)*(x[3*iID+2] - m_positionZ);
      if(distanceSquared <= minDistanceSquared){
        int globalId = static_cast<int>(id[iID]);
        // If there are ties, choose the element with the lowest global id.
        if(distanceSquared == minDistanceSquared)
          foundTies = true;
        else
          foundTies = false;
        if(distanceSquared < minDistanceSquared || globalId < m_elementId){
          m_elementId = globalId;
          minDistanceSquared = distanceSquared;
        }
      }
    }
  }

  // Parallel communication to find closest point across all processors.
  vector<double> localMinDistanceSquared(1), globalMinDistanceSquared(1);
  localMinDistanceSquared[0] = minDistanceSquared;
  epetraComm()->MinAll(&localMinDistanceSquared[0], &globalMinDistanceSquared[0], 1);
  if(minDistanceSquared != globalMinDistanceSquared[0])
    m_elementId = INT_MAX;

  // If there are ties, choose the element with the lowest global id.
  vector<int> localData(1), globalData(1);
  localData[0] = m_elementId;
  epetraComm()->MinAll(&localData[0], &globalData[0], 1);
  m_elementId = globalData[0];
  if(localData[0] != INT_MAX && localData[0] != m_elementId)
    foundTies = true;

  // To make tracking more efficient, determine which block has the element.
  localData[0] = 0;
  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    if(blockIt->getOwnedScalarPointMap()->LID(m_elementId) != -1)
      localData[0] = blockIt->getID();
  }
  epetraComm()->SumAll(&localData[0], &globalData[0], 1);
  m_blockId = globalData[0];

  // If ties were found, print a warning
  localData[0] = static_cast<int>(foundTies);
  epetraComm()->SumAll(&localData[0], &globalData[0], 1);
  if(globalData[0] > 0 && epetraComm->MyPID() == 0){
    cout << "**** Warning:  The Nearest_Neighbor_Data compute class found multiple nearest neighbors." << endl;
    cout << "****           The element with the smallest global ID will be selected for tracking.\n" << endl;
  }

  // If verbose flag is set, write output to screen
  if(m_verbose){
    
    // Find the coordinates of the element that will be tracked
    vector<double> localValues(3), globalValues(3);
    for(int i=0 ; i<3 ; ++i)
      localValues[0] = 0.0;
    for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      int localId = blockIt->getOwnedScalarPointMap()->LID(m_elementId);
      if(localId != -1){
        Teuchos::RCP<Epetra_Vector> modelCoordinates = blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE);
        localValues[0] = (*modelCoordinates)[3*localId];
        localValues[1] = (*modelCoordinates)[3*localId+1];
        localValues[2] = (*modelCoordinates)[3*localId+2];
      }
    }
    epetraComm()->SumAll(&localValues[0], &globalValues[0], 3);

    // Write to the root processor
    if(epetraComm->MyPID() == 0){
      stringstream ss;
      ss << "Nearest Point Data Compute Class:" << endl;
      ss << "  Requested variable: " << m_variable << endl;
      ss << "  Requested location: " << m_positionX << ", " << m_positionY << ", " << m_positionZ << endl;
      ss << "  Closest Element Id: " << m_elementId << endl;
      ss << "  Closest Element Position: " << globalValues[0] << ", " << globalValues[1] << ", " << globalValues[2] << endl;
      cout << ss.str() << endl;
    }
  }
}

int PeridigmNS::Compute_Nearest_Point_Data::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  PeridigmField::Step step = PeridigmField::STEP_NONE;
  if(m_variableIsStated)
    step = PeridigmField::STEP_NP1;
  
  vector<double> localData(3), globalData(3);
  localData[0] = 0.0;

  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    if(blockIt->getID() == m_blockId){
      int localId = blockIt->getOwnedScalarPointMap()->LID(m_elementId);
      if(localId != -1){
        Teuchos::RCP<Epetra_Vector> data = blockIt->getData(m_variableFieldId, step);
        if(m_variableLength == 1){
          localData[0] = (*data)[localId];
        }
        else if(m_variableLength == 3){
          localData[0] = (*data)[3*localId];
          localData[1] = (*data)[3*localId+1];
          localData[2] = (*data)[3*localId+2];
        }
      }
    }
  }

  if(m_variableLength == 1){
    epetraComm()->SumAll(&localData[0], &globalData[0], 1);
    double& outputData = blocks->begin()->getGlobalData(m_outputFieldId);
    outputData = globalData[0];
  }
  else if(m_variableLength == 3){
    epetraComm()->SumAll(&localData[0], &globalData[0], 3);
    double& outputDataX = blocks->begin()->getGlobalData(m_outputXFieldId);
    outputDataX = globalData[0];
    double& outputDataY = blocks->begin()->getGlobalData(m_outputYFieldId);
    outputDataY = globalData[1];
    double& outputDataZ = blocks->begin()->getGlobalData(m_outputZFieldId);
    outputDataZ = globalData[2];
  }

  return 0;
}
