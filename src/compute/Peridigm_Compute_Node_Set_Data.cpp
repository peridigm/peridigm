/*! \file Peridigm_Compute_Node_Set_Data.cpp */

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

#include "Peridigm_Compute_Node_Set_Data.hpp"
#include "Peridigm_Discretization.hpp"
#include "Peridigm_Field.hpp"

//! Standard constructor.
PeridigmNS::Compute_Node_Set_Data::Compute_Node_Set_Data(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                         Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                         Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_), m_calculationType(UNDEFINED_CALCULATION),
    m_variableFieldId(-1), m_outputFieldId(-1)
{
  m_nodeSetName = params->get<string>("Node Set");
  m_variable = params->get<string>("Variable");
  m_outputLabel = params->get<string>("Output Label");

  // Obtain the global ids for all the locally-owned nodes in the specified node set
  Teuchos::RCP<Discretization> disc = *( computeClassGlobalData_->get<Teuchos::RCP<Discretization>*>("discretization") );
  Teuchos::RCP< std::map< std::string, std::vector<int> > > nodeSets = disc->getNodeSets();

  // Make sure the requested node set exists (it's really easy for a user to get confused between nodeset_1 and nodelist_1, etc.)
  std::string msg = "**** Error:  Invalid \"Node Set\" in Node_Set_Data compute class, specified node set not found.\n";
              msg += "             Requested node set: " + m_nodeSetName + "\n";
              msg += "             Available node sets:";
  for(std::map< std::string, std::vector<int> >::iterator it=nodeSets->begin() ; it!=nodeSets->end() ; it++)
    msg += " " + it->first;
  msg += "\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(m_nodeSetName) == nodeSets->end(), msg);

  std::vector<int> nodeSet = (*nodeSets)[m_nodeSetName];
  for(unsigned int i=0 ; i<nodeSet.size() ; ++i)
    m_nodeSet.insert(nodeSet[i]);

  string calculationType = params->get<string>("Calculation Type");
  if(calculationType == "Minimum"){
    m_calculationType = MINIMUM;
  }
  else if(calculationType == "Maximum"){
    m_calculationType = MAXIMUM;
  }
  else if(calculationType == "Sum"){
    m_calculationType = SUM;
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, 
     "**** Error:  invalid \"Calculation Type\" in Node_Set_Data compute class, must be \"Minimum\", \"Maximum\", or \"Sum\".\n");
  }

  FieldManager& fieldManager = FieldManager::self();
  m_variableFieldId = fieldManager.getFieldId(m_variable);
  m_fieldIds.push_back(m_variableFieldId);

  PeridigmField::Length length = fieldManager.getFieldSpec(m_variableFieldId).getLength();
  if(length == PeridigmField::SCALAR)
    m_variableLength = 1;
  else if(length == PeridigmField::VECTOR)
    m_variableLength = 3;
  else
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error:  Node_Set_Data compute class can be called only for SCALAR or VECTOR data.\n");

  if(m_variableLength == 1)
    m_outputFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, m_outputLabel);
  else if(m_variableLength == 3)
    m_outputFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::VECTOR, PeridigmField::CONSTANT, m_outputLabel);

  m_fieldIds.push_back(m_outputFieldId);

  PeridigmField::Temporal temporal = fieldManager.getFieldSpec(m_variableFieldId).getTemporal();
  m_variableIsStated = false;
  if(temporal == PeridigmField::TWO_STEP)
    m_variableIsStated = true;
}

void PeridigmNS::Compute_Node_Set_Data::initialize( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) {
  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    std::string blockName = blockIt->getName();
    m_blockLocalIds[blockName] = std::vector<int>();
    int numOwnedPoints = blockIt->getNeighborhoodData()->NumOwnedPoints();
    Teuchos::RCP<const Epetra_BlockMap> map = blockIt->getOwnedScalarPointMap();
    for(int i=0 ; i<numOwnedPoints ; ++i){
      int globalId = map->GID(i);
      if( m_nodeSet.find(globalId) != m_nodeSet.end() )
        m_blockLocalIds[blockName].push_back(i);
    }
  }
}

int PeridigmNS::Compute_Node_Set_Data::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  PeridigmField::Step step = PeridigmField::STEP_NONE;
  if(m_variableIsStated)
    step = PeridigmField::STEP_NP1;
  
  std::vector<double> localData(3), globalData(3);

  if(m_calculationType == MINIMUM){
    for(int i=0 ; i<3 ; ++i)
      localData[i] = DBL_MAX;
  }
  else if(m_calculationType == MAXIMUM){
    for(int i=0 ; i<3 ; ++i)
      localData[i] = DBL_MIN;
  }
  else if(m_calculationType == SUM){
    for(int i=0 ; i<3 ; ++i)
      localData[i] = 0.0;
  }

  for(std::vector<Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    double *data;
    blockIt->getData(m_variableFieldId, step)->ExtractView(&data);
    std::string blockName = blockIt->getName();
    const std::vector<int>& localNodeIds = m_blockLocalIds.at(blockName);
    int localId;

    if(m_calculationType == MINIMUM){
      for(unsigned int i=0 ; i<localNodeIds.size() ; ++i){
        localId = localNodeIds[i];
        if(m_variableLength == 1){
          if(data[localId] < localData[0])
            localData[0] = data[localId];
        }
        else if(m_variableLength == 3){
          if(data[3*localId] < localData[0])
            localData[0] = data[3*localId];
          if(data[3*localId+1] < localData[1])
            localData[1] = data[3*localId+1];
          if(data[3*localId+2] < localData[2])
            localData[2] = data[3*localId+2];
        }
      }
    }
    else if(m_calculationType == MAXIMUM){
      for(unsigned int i=0 ; i<localNodeIds.size() ; ++i){
        localId = localNodeIds[i];
        if(m_variableLength == 1){
          if(data[localId] > localData[0])
            localData[0] = data[localId];
        }
        else if(m_variableLength == 3){
          if(data[3*localId] > localData[0])
            localData[0] = data[3*localId];
          if(data[3*localId+1] > localData[1])
            localData[1] = data[3*localId+1];
          if(data[3*localId+2] > localData[2])
            localData[2] = data[3*localId+2];
        }
      }
    }
    else if(m_calculationType == SUM){
      for(unsigned int i=0 ; i<localNodeIds.size() ; ++i){
        localId = localNodeIds[i];
        if(m_variableLength == 1){
          localData[0] += data[localId];
        }
        else if(m_variableLength == 3){
          localData[0] += data[3*localId];
          localData[1] += data[3*localId+1];
          localData[2] += data[3*localId+2];
        }
      }
    }
  }

  Teuchos::RCP<Epetra_Vector> outputData = blocks->begin()->getData(m_outputFieldId, PeridigmField::STEP_NONE);
  if(m_variableLength == 1){
    if(m_calculationType == MINIMUM)
      epetraComm()->MinAll(&localData[0], &globalData[0], 1);
    else if(m_calculationType == MAXIMUM)
      epetraComm()->MaxAll(&localData[0], &globalData[0], 1);
    else if(m_calculationType == SUM)
      epetraComm()->SumAll(&localData[0], &globalData[0], 1);
    (*outputData)[0] = globalData[0];
  }
  else if(m_variableLength == 3){
    if(m_calculationType == MINIMUM)
      epetraComm()->MinAll(&localData[0], &globalData[0], 3);
    else if(m_calculationType == MAXIMUM)
      epetraComm()->MaxAll(&localData[0], &globalData[0], 3);
    else if(m_calculationType == SUM)
      epetraComm()->SumAll(&localData[0], &globalData[0], 3);
    (*outputData)[0] = globalData[0];
    (*outputData)[1] = globalData[1];
    (*outputData)[2] = globalData[2];
  }

  return 0;
}
