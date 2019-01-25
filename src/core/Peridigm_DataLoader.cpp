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

PeridigmNS::DataLoader::DataLoader(const Teuchos::ParameterList& contactParams,
                                   Teuchos::RCP<const Epetra_BlockMap> epetraMap)
  : fileName("none"), fieldName("none"), exodusName("none"), exodusVariableIndex(-1), numRanks(-1), myRank(-1)
{
  fileName = contactParams.get<std::string>("File Name");
  fieldName = contactParams.get<std::string>("Field Name");
  numRanks = epetraMap->Comm().NumProc();
  myRank = epetraMap->Comm().MyPID();
  int vecLength = epetraMap->NumMyElements();
  scratchArray.resize(vecLength);
  scratchVector = Teuchos::rcp(new Epetra_Vector(*epetraMap));

  // Assume field is nodal, scalar, two-step
  fieldId = PeridigmNS::FieldManager::self().getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, fieldName);

  // Append processor id information to the file name, if necessary
  exodusName = fileName;
  if(numRanks != 1){
    int width = 1;
    if(numRanks > 9)
      width = 2;
    if(numRanks > 99)
      width = 3;
    if(numRanks > 999)
      width = 4;
    if(numRanks > 9999)
      width = 5;
    if(numRanks > 99999)
      width = 6;
    if(numRanks > 999999)
      width = 7;
    if(numRanks > 9999999)
      width = 8;
    if(numRanks > 99999999)
      width = 9;
    std::stringstream ss;
    ss << "." << numRanks << "." << std::setfill('0') << std::setw(width) << myRank;
    exodusName += ss.str();
  }

  // Open the exodus file
  int compWordSize = sizeof(double);
  int ioWordSize = 0;
  float exodusVersion;
  int exodusFileId = ex_open(exodusName.c_str(), EX_READ, &compWordSize, &ioWordSize, &exodusVersion);
  if(exodusFileId < 0){
    std::cout << "\n****Error on processor " << myRank << ": unable to open file " << exodusName.c_str() << "\n" << std::endl;
    reportExodusError(exodusFileId, "DataLoader()", "ex_open");
  }

  // Read the initialization parameters
  int numDim, numNodes, numElem, numElemBlocks, numNodeSets, numSideSets;
  char title[MAX_LINE_LENGTH];
  int retval = ex_get_init(exodusFileId, title, &numDim, &numNodes, &numElem, &numElemBlocks, &numNodeSets, &numSideSets);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_init");

  int numNodalVariables;
  retval = ex_get_variable_param(exodusFileId, EX_NODAL, &numNodalVariables);
  if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_variable_param");

  exodusVariableIndex = -1;
  for (int varIndex=1 ; varIndex<numNodalVariables+1 ; ++varIndex){
    char variableName[MAX_STR_LENGTH+1];
    retval = ex_get_variable_name(exodusFileId, EX_NODAL, varIndex, variableName);
    if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_variable_name");
    if(std::string(variableName) == fieldName){
      exodusVariableIndex = varIndex;
    }
  }

  if (exodusVariableIndex == -1){
    std::string msg = "**** Error in DataLoader::DataLoader(), requested variable name " + fieldName + " not found.\n";
    msg += "**** Nodal variables in " + fileName + " are:\n";
    for (int varIndex=1 ; varIndex<numNodalVariables+1 ; ++varIndex){
      char variableName[MAX_STR_LENGTH+1];
      retval = ex_get_variable_name(exodusFileId, EX_NODAL, varIndex, variableName);
      if (retval != 0) reportExodusError(retval, "DataLoader()", "ex_get_variable_name");
      msg += "****   " + std::string(variableName);
    }
    msg += "\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(exodusVariableIndex == -1, msg);
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
  fieldIds.push_back(fieldId);
  return fieldIds;
}

void PeridigmNS::DataLoader::loadDataFromFile(int step)
{
  int compWordSize = sizeof(double);
  int ioWordSize = 0;
  float exodusVersion;
  int exodusFileId = ex_open(exodusName.c_str(), EX_READ, &compWordSize, &ioWordSize, &exodusVersion);
  if(exodusFileId < 0){
    std::cout << "\n****Error on processor " << myRank << ": unable to open file " << exodusName.c_str() << "\n" << std::endl;
    reportExodusError(exodusFileId, "loadData()", "ex_open");
  }

  int objId = 0;
  int retval = ex_get_var(exodusFileId, step, EX_NODAL, exodusVariableIndex, objId, scratchArray.size(), scratchArray.data());
  if (retval != 0) reportExodusError(retval, "loadData()", "ex_get_var");

  for (unsigned int i=0 ; i<scratchArray.size() ; i++) {
    (*scratchVector)[i] = scratchArray[i];
  }

  // Close the genesis file
  retval = ex_close(exodusFileId);
  if (retval != 0) reportExodusError(retval, "loadData()", "ex_close");
}

void PeridigmNS::DataLoader::copyDataToDataManagers(Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks)
{
  for(std::vector<PeridigmNS::Block>::iterator blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*scratchVector, fieldId, PeridigmField::STEP_NP1, Insert);
  }
}

void PeridigmNS::DataLoader::reportExodusError(int errorCode, const char *methodName, const char*exodusMethodName)
{
  std::stringstream ss;
  if (errorCode < 0) { // error
    if (numRanks > 1) ss << "Error on PID #" << myRank << ": ";
    ss << "PeridigmNS::DataLoader" << methodName << "() -- Error code: " << errorCode << " (" << exodusMethodName << ")";
    TEUCHOS_TEST_FOR_EXCEPTION(1, std::invalid_argument, ss.str());
  }
  else {
    if (numRanks > 1) ss << "Warning on PID #" << myRank << ": ";
    ss << "PeridigmNS::DataLoader::" << methodName << "() -- Warning code: " << errorCode << " (" << exodusMethodName << ")";
    std::cout << ss.str() << std::endl;
  }
}
