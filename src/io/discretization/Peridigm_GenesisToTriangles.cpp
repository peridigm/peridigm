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

#include "Peridigm_GenesisToTriangles.hpp"
#include "Peridigm_Enums.hpp"
#include <exodusII.h>
#include <sstream>
#include <iostream>

void PeridigmNS::GenesisToTriangles(std::string genesis_file_name,
                                    std::vector< std::vector< std::vector<double> > > & triangles) {

  // Open the genesis file
  int compWordSize = sizeof(double);
  int ioWordSize = 0;
  float exodusVersion;
  int exodusFileId = ex_open(genesis_file_name.c_str(), EX_READ, &compWordSize, &ioWordSize, &exodusVersion);
  if(exodusFileId < 0){
    std::cout << "\n****Error: unable to open file " << genesis_file_name << "\n" << std::endl;
    report_exodus_error(exodusFileId, "ex_open");
  }

  // Read the initialization parameters
  int numDim, numNodes, numElem, numElemBlocks, numNodeSets, numSideSets;
  char title[MAX_LINE_LENGTH];
  int retval = ex_get_init(exodusFileId, title, &numDim, &numNodes, &numElem, &numElemBlocks, &numNodeSets, &numSideSets);
  if (retval != 0) report_exodus_error(retval, "ex_get_init");

  // Node coordinates
  std::vector<double> exodusNodeCoordX(numNodes), exodusNodeCoordY(numNodes), exodusNodeCoordZ(numNodes);
  retval = ex_get_coord(exodusFileId, &exodusNodeCoordX[0], &exodusNodeCoordY[0], &exodusNodeCoordZ[0]);
  if (retval != 0) report_exodus_error(retval, "ex_get_coord");

  // Process the element blocks
  std::vector<int> elemBlockIds(numElemBlocks);
  retval = ex_get_elem_blk_ids(exodusFileId, &elemBlockIds[0]);
  if (retval != 0) report_exodus_error(retval, "ex_get_elem_blk_ids");

  for(int iElemBlock=0 ; iElemBlock<numElemBlocks ; iElemBlock++){

    int elemBlockId = elemBlockIds[iElemBlock];

    // Get the block parameters and the element connectivity
    char elemType[MAX_STR_LENGTH];
    int numElemThisBlock, numNodesPerElem, numAttributes;
    retval = ex_get_elem_block(exodusFileId, elemBlockId, elemType, &numElemThisBlock, &numNodesPerElem, &numAttributes);
    if (retval != 0) report_exodus_error(retval, "ex_get_elem_block");
    std::vector<int> conn;
    if(numElemThisBlock > 0){
      std::string elemTypeString(elemType);
      to_upper(elemTypeString);
      std::stringstream ss;
      ss << "**** Error in GenesisToTriangles(), invalid element type: " << elemTypeString << ".\n";
      TEUCHOS_TEST_FOR_EXCEPTION(elemTypeString != "TRI" && elemTypeString != "TRI3", std::invalid_argument, ss.str());
      conn.resize(numElemThisBlock*numNodesPerElem);
      retval = ex_get_elem_conn(exodusFileId, elemBlockId, &conn[0]);
      if (retval != 0) report_exodus_error(retval, "ex_get_elem_conn");
    }

    // Loop over the elements in the block

    for(int iElem=0 ; iElem<numElemThisBlock ; iElem++){

      // Node coordinates for all the nodes in this element
      std::vector< std::vector<double> > triangle;
      for(int i=0 ; i<numNodesPerElem ; ++i){
        int nodeId = conn[iElem*numNodesPerElem + i] - 1;
        std::vector<double> point(3);
        point[0] = exodusNodeCoordX[nodeId];
        point[1] = exodusNodeCoordY[nodeId];
        point[2] = exodusNodeCoordZ[nodeId];
        triangle.push_back(point);
      }
      triangles.push_back(triangle);
    }
  }

  // Close the genesis file
  retval = ex_close(exodusFileId);
  if (retval != 0) report_exodus_error(retval, "ex_close");
}

void PeridigmNS::report_exodus_error(int errorCode, const char*exodusMethodName)
{
  std::stringstream ss;
  if (errorCode < 0) { // error
    ss << "GenesisToTriangles() -- Error code: " << errorCode << " (" << exodusMethodName << ")";
    TEUCHOS_TEST_FOR_EXCEPTION(1, std::invalid_argument, ss.str());
  }
  else {
    ss << "GenesisToTriangles() -- Warning code: " << errorCode << " (" << exodusMethodName << ")";
    std::cout << ss.str() << std::endl;
  }
}
