/*! \file Peridigm_DataLoader.cpp */

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

#include "Peridigm_DataLoader.hpp"
#include <Epetra_MpiComm.h>
#include "Teuchos_Assert.hpp"
#include <exodusII.h>
#include <sstream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <unistd.h>

PeridigmNS::DataLoader::DataLoader(const Teuchos::ParameterList& contactParams,
                                   Teuchos::RCP<const Epetra_BlockMap> epetraMap)
  : fileName_("none"), fieldName_("none"), fieldId_(-1), exodusName_("none"), exodusVariableIndex_(-1),
    numRanks_(-1), myRank_(-1), time_(-1.0), time_1_(-1.0), time_2_(-1.0)
{
  fileName_ = contactParams.get<std::string>("File Name");
  fieldName_ = contactParams.get<std::string>("Field Name");
  numRanks_ = epetraMap->Comm().NumProc();
  myRank_ = epetraMap->Comm().MyPID();
  int vecLength = epetraMap->NumMyElements();
  data_1_.resize(vecLength);
  data_2_.resize(vecLength);
  scratch_ = Teuchos::rcp(new Epetra_Vector(*epetraMap));

  // Assume field is nodal, scalar, two-step
  fieldId_ = PeridigmNS::FieldManager::self().getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, fieldName_);

  // Append processor id information to the file name, if necessary
  exodusName_ = fileName_;
  if(numRanks_ != 1){
    int width = 1;
    if(numRanks_ > 9)
      width = 2;
    if(numRanks_ > 99)
      width = 3;
    if(numRanks_ > 999)
      width = 4;
    if(numRanks_ > 9999)
      width = 5;
    if(numRanks_ > 99999)
      width = 6;
    if(numRanks_ > 999999)
      width = 7;
    if(numRanks_ > 9999999)
      width = 8;
    if(numRanks_ > 99999999)
      width = 9;
    std::stringstream ss;
    ss << "." << numRanks_ << "." << std::setfill('0') << std::setw(width) << myRank_;
    exodusName_ += ss.str();
  }

  // Open the exodus file
  int compWordSize = sizeof(double);
  int ioWordSize = 0;
  float exodusVersion;
  int exodusFileId = ex_open(exodusName_.c_str(), EX_READ, &compWordSize, &ioWordSize, &exodusVersion);
  if(exodusFileId < 0){
    std::cout << "\n****Error on processor " << myRank_ << ": unable to open file " << exodusName_.c_str() << "\n" << std::endl;
    reportExodusError(exodusFileId, "DataLoader()", "ex_open");
  }

  // Read the initialization parameters
  int numDim, numNodes, numElem, numElemBlocks, numNodeSets, numSideSets;
  char title[MAX_LINE_LENGTH];
  int retval = ex_get_init(exodusFileId, title, &numDim, &numNodes, &numElem, &numElemBlocks, &numNodeSets, &numSideSets);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_init");

  // Global node numbering
  std::vector<int> nodeIdMap(numNodes);
  retval = ex_get_id_map(exodusFileId, EX_NODE_MAP, &nodeIdMap[0]);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_id_map");
  for(int i=0 ; i<numNodes ; ++i)
    nodeIdMap[i] -= 1; // Note the switch from 1-based indexing to 0-based indexing

  // Check for auxiliary node maps
  int numNodeMaps, numElemMaps;
  retval = ex_get_map_param(exodusFileId, &numNodeMaps, &numElemMaps);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_map_param");

  // DJL
  // This block of code handles the case where an extra elem or node map
  // called "original_global_id_map" is supplied.  This can be the case
  // for parallel decompositions created with decomp or loadbal.
  // If there is an auxiliary map provided that has a different name, throw an
  // error because I don't know what to do with it.
  if(numNodeMaps > 0){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numNodeMaps > 1,
                                "**** Error in DataLoader(), genesis file contains invalid number of auxiliary node maps (>1).\n");
    char mapName[MAX_STR_LENGTH];
    retval = ex_get_name(exodusFileId, EX_NODE_MAP, 1, mapName);
    if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_name");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(std::string(mapName) != std::string("original_global_id_map"),
                                "**** Error in DataLoader::DataLoader(), unknown exodus EX_NODE_MAP: " + std::string(mapName) + ".\n");
    std::vector<int> auxMap(numNodes);
    retval = ex_get_num_map(exodusFileId, EX_NODE_MAP, 1, &auxMap[0]);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_num_map");
    for(int i=0 ; i<numNodes ; ++i)
      auxMap[i] -= 1; // Note the switch from 1-based indexing to 0-based indexing
    // Use original_global_id_map instead of the map returned by ex_get_id_map()
    nodeIdMap = auxMap;
  }

  // Ensure that the node map in the data file matches the node map in the mesh file
  for (int i=0 ; i<numNodes; ++i) {
    int epetraMapGlobalId = epetraMap->GID(i);
    int dataFileGlobalId = nodeIdMap[i];
    TEUCHOS_TEST_FOR_EXCEPT_MSG(epetraMapGlobalId != dataFileGlobalId, "**** Error in DataLoad(), node map in mesh file and node map in data file do not match! (Accidental misuse of SEACAS tools decomp/epu/algebra?)\n");
  }

  int numNodalVariables;
  retval = ex_get_variable_param(exodusFileId, EX_NODAL, &numNodalVariables);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_variable_param");

  // Determine the index for the requested variable
  exodusVariableIndex_ = -1;
  for (int varIndex=1 ; varIndex<numNodalVariables+1 ; ++varIndex){
    char variableName[MAX_STR_LENGTH+1];
    retval = ex_get_variable_name(exodusFileId, EX_NODAL, varIndex, variableName);
    if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_variable_name");
    if(std::string(variableName) == fieldName_){
      exodusVariableIndex_ = varIndex;
    }
  }

  if (exodusVariableIndex_ == -1){
    std::string msg = "**** Error in DataLoader::DataLoader(), requested variable name " + fieldName_ + " not found.\n";
    msg += "**** Nodal variables in " + fileName_ + " are:\n";
    for (int varIndex=1 ; varIndex<numNodalVariables+1 ; ++varIndex){
      char variableName[MAX_STR_LENGTH+1];
      retval = ex_get_variable_name(exodusFileId, EX_NODAL, varIndex, variableName);
      if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_variable_name");
      msg += "****   " + std::string(variableName);
    }
    msg += "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(exodusVariableIndex_ == -1, msg);
  }

  // Close the genesis file
  retval = ex_close(exodusFileId);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_close");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(numNodes != vecLength,
                              "**** Error in DataLoader::DataLoader(), unexpected array length.\n");
}

std::vector<int> PeridigmNS::DataLoader::getFieldIds() const
{
  std::vector<int> fieldIds;
  fieldIds.push_back(fieldId_);
  return fieldIds;
}

void PeridigmNS::DataLoader::loadData(double time,
                                      Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks)
{
  if (time != time_) {

    time_ = time;

    if (time_ < time_1_ || time_ > time_2_) {
      // data must be read from disk
      //waitForData();
      loadDataFromFile();
    }

    // interpolate data and copy to scratch
    double c = (time_ - time_1_) / (time_2_ - time_1_);
    for (unsigned int i=0 ; i<data_1_.size() ; ++i) {
      (*scratch_)[i] = (1.0 - c)*data_1_[i] + c*data_2_[i];
    }
  }

  // copy data into the block DataManagers
  for(std::vector<PeridigmNS::Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++) {
    blockIt->importData(scratch_, fieldId_, PeridigmField::STEP_NP1, Insert);
  }
}

void PeridigmNS::DataLoader::waitForData()
{
  double pause_time = 5.0; // seconds
  double total_pause_time = 0.0;
  double max_pause_time = 12.0;
  bool success = false;
  bool time_limit = false;

  std::string ready_file_name("PFLOTRAN_DONE");
  std::ifstream ready_file(ready_file_name);
  if (ready_file) {
    std::cout << "FILE EXISTS" << std::endl;
    ready_file.close();
    remove(ready_file_name.c_str());
  }
  else {
    std::cout << "FILE DOES NOT EXIST" << std::endl;
  }


  std::string line;
  while(true) {
    std::cout << "about to read" << std::endl;
    std::ifstream steps_file("Notched_Plate_Data_Loader_Steps.txt");
    while(std::getline(steps_file, line)) {
      std::stringstream ss(line);
      double time;
      ss >> time;
      std::cout << "Looking for time_ " << time_ << " currently have " << time << std::endl;
      if (time > time_) {
        success = true;
        break;
      }
    }
    steps_file.close();
    std::cout << "end of read" << std::endl;
    if (success) {
      break;
    }
    else {
      if (total_pause_time > max_pause_time) {
        time_limit = true;
        break;
      }
      total_pause_time += pause_time;
      // microseconds
      double pause_time_microseconds = pause_time * 1.0e6;
      std::cout << "pausing" << std::endl;
      usleep(pause_time_microseconds);
    }
  }

  if(time_limit) {
    std::stringstream ss;
    ss << "**** Error in DataLoader::waitForData(), exceeded maximum wait time of " << max_pause_time << " seconds.\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(time_limit, ss.str());
  }
}

void PeridigmNS::DataLoader::loadDataFromFile()
{
  int compWordSize = sizeof(double);
  int ioWordSize = 0;
  float exodusVersion;
  int exodusFileId = ex_open(exodusName_.c_str(), EX_READ, &compWordSize, &ioWordSize, &exodusVersion);
  if(exodusFileId < 0){
    std::cout << "\n****Error on processor " << myRank_ << ": unable to open file " << exodusName_.c_str() << "\n" << std::endl;
    reportExodusError(exodusFileId, "loadDataFromFile()", "ex_open");
  }

  int intReturnValue;
  float floatReturnValue;
  char charReturnValue;
  int retval = ex_inquire(exodusFileId, EX_INQ_TIME, &intReturnValue, &floatReturnValue, &charReturnValue);
  if (retval != 0) reportExodusError(retval, "loadDataFromFile()", "ex_inquiry");
  int numTimeSteps = intReturnValue;
  std::vector<double> times(numTimeSteps);
  retval = ex_get_all_times(exodusFileId, times.data());
  if (retval != 0) reportExodusError(retval, "loadDataFromFile()", "ex_get_all_times");

  int step_1(-1), step_2(-1);
  for (int i=0 ; i<numTimeSteps-1 ; i++) {
    if (time_ >= times[i] && time_ <= times[i+1]) {
      step_1 = i+1;
      step_2 = i+2;
      time_1_ = times[i];
      time_2_ = times[i+1];
      break;
    }
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG((step_1 == -1 || step_2 == -1), "Error in DataLoader::loadDataFromFile(), requested time is outside the range of times in the data file.\n");

  int objId = 0;
  retval = ex_get_var(exodusFileId, step_1, EX_NODAL, exodusVariableIndex_, objId, data_1_.size(), data_1_.data());
  if (retval != 0) reportExodusError(retval, "loadDataFromFile()", "ex_get_var");

  retval = ex_get_var(exodusFileId, step_2, EX_NODAL, exodusVariableIndex_, objId, data_2_.size(), data_2_.data());
  if (retval != 0) reportExodusError(retval, "loadDataFromFile()", "ex_get_var");

  // Close the genesis file
  retval = ex_close(exodusFileId);
  if (retval != 0) reportExodusError(retval, "loadDataFromFile()", "ex_close");
}

void PeridigmNS::DataLoader::reportExodusError(int errorCode, const char *methodName, const char*exodusMethodName)
{
  std::stringstream ss;
  if (errorCode < 0) { // error
    if (numRanks_ > 1) ss << "Error on PID #" << myRank_ << ": ";
    ss << "PeridigmNS::DataLoader" << methodName << "() -- Error code: " << errorCode << " (" << exodusMethodName << ")";
    TEUCHOS_TEST_FOR_EXCEPTION(1, std::invalid_argument, ss.str());
  }
  else {
    if (numRanks_ > 1) ss << "Warning on PID #" << myRank_ << ": ";
    ss << "PeridigmNS::DataLoader::" << methodName << "() -- Warning code: " << errorCode << " (" << exodusMethodName << ")";
    std::cout << ss.str() << std::endl;
  }
}
