/*! \file Peridigm_InterfaceData.hpp */

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

#ifndef PERIDIGM_INTERFACEDATA_HPP
#define PERIDIGM_INTERFACEDATA_HPP

#include <vector>
#include <Teuchos_RCP.hpp>
#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_BlockMap.h>

namespace PeridigmNS {

class InterfaceData {

public:

  InterfaceData(): numOwnedPoints(0), ownedIDs(0), elementLeft(0), elementRight(0), numNodes(0), exoid(0), numQuads(0), numTris(0){}

  ~InterfaceData(){
    if(ownedIDs != 0)
      delete[] ownedIDs;
    if(elementLeft != 0)
      delete[] elementLeft;
    if(elementRight != 0)
      delete[] elementRight;
    if(numNodes != 0)
      delete[] numNodes;
  }

  void Initialize(std::vector<int> leftElements, std::vector<int> rightElements, std::vector<int> numNodesPerElem,
    std::vector<std::vector<int> > interfaceNodesVec, const Teuchos::RCP<const Epetra_Comm> & Comm);

  void InitializeExodusOutput(Teuchos::RCP<Epetra_Vector> exodusMeshElementConnectivity, Teuchos::RCP<Epetra_Vector> exodusMeshNodePositions);

  void WriteExodusOutput(int timeStep, const float & timeValue, Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector> y);

  int NumOwnedPoints() const{
  return numOwnedPoints;
  }

  int* OwnedIDs() const{
  return ownedIDs;
  }

  int* ElementLeft() const{
  return elementLeft;
  }

  int* ElementRight() const{
  return elementRight;
  }

  Teuchos::RCP<Epetra_BlockMap> Map() const{
    return interfaceMap;
  }

protected:
  int numOwnedPoints;
  int* ownedIDs;
  int* elementLeft;
  int* elementRight;
  int* numNodes;
  int exoid;
  std::stringstream filename;
  int numQuads;
  int numTris;
  Teuchos::RCP<const Epetra_Comm> comm;
  Teuchos::RCP<Epetra_BlockMap> interfaceMap;
  Teuchos::RCP<Epetra_BlockMap> interfaceNodesMap;
  Teuchos::RCP<Epetra_BlockMap> elemOverlapMap;
  Teuchos::RCP<Epetra_Vector> interfaceAperture;
  Teuchos::RCP<Epetra_Vector> interfaceNodes;

};

}

#endif // PERIDIGM_INTERFACEDATA_HPP
