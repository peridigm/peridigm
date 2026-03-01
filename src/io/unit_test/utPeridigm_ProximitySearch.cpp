/*! \file utPeridigm_ProximitySearch.cpp */

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

#include <set>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>
#include "Peridigm_ProximitySearch.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <vector>

using namespace Teuchos;
using namespace PeridigmNS;
using std::vector;
using std::map;
using std::pair;
using std::cout;
using std::set;


TEUCHOS_UNIT_TEST(ProximitySearch, TwoPointProblem) {

  Teuchos::RCP<Epetra_Comm> comm;
#ifdef HAVE_MPI
  comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
  int numProc = comm->NumProc();

  // This test cannot be run on more than 2 processors
  if(numProc > 2){
    std::cerr << "Unit test runtime ERROR: utPeridigm_Compute_Force only makes sense on 1 to 2 processors." << std::endl;
    return;
  }


  Epetra_BlockMap map(2, 3, 0, *comm);
  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(map));

  int numMyElements = map.NumMyElements();

  if(numProc == 1)
    TEST_ASSERT(numMyElements == 2);
  if(numProc == 2)
    TEST_ASSERT(numMyElements == 1);

  vector<double> node(3);
  std::map<int, vector<double> > nodes;
  node[0] =  1.2 ; node[1] =  1.3 ; node[2] =  6.0 ; nodes[0] = node;
  node[0] =  4.2 ; node[1] = -2.1 ; node[2] = -3.8 ; nodes[1] = node;

  double distance = sqrt( (nodes[0][0] - nodes[1][0])*(nodes[0][0] - nodes[1][0]) + 
                          (nodes[0][1] - nodes[1][1])*(nodes[0][1] - nodes[1][1]) + 
                          (nodes[0][2] - nodes[1][2])*(nodes[0][2] - nodes[1][2]) );

  for(int i=0 ; i<numMyElements ; ++i){
    int globalId = map.GID(i);
    (*x)[3*i]   = nodes[globalId][0];
    (*x)[3*i+1] = nodes[globalId][1];
    (*x)[3*i+2] = nodes[globalId][2];
  }

  // These are filled by the proximity search
  Teuchos::RCP<Epetra_BlockMap> overlapMap;
  int neighborListSize(0);
  int* neighborList(0);

  // ---- Call the proximity search with a radius just below the distance ----

  double searchRadius = distance - 1.0e-10;
  Epetra_BlockMap oneDimensionalMap(map.NumGlobalElements(), map.NumMyElements(), map.MyGlobalElements(), 1, 0, *comm);
  Teuchos::RCP<Epetra_Vector> searchRadii = Teuchos::rcp(new Epetra_Vector(oneDimensionalMap));
  searchRadii->PutScalar(searchRadius);

  ProximitySearch::GlobalProximitySearch(x, searchRadii, overlapMap, neighborListSize, neighborList);

  int neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    TEST_ASSERT(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      TEST_ASSERT(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 1){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
  }

  // ---- Call the proximity search with a radius just above the distance ----

  searchRadius = distance + 1.0e-10;
  searchRadii->PutScalar(searchRadius);

  ProximitySearch::GlobalProximitySearch(x, searchRadii, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    TEST_ASSERT(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      TEST_ASSERT(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      TEST_ASSERT(numNeighbors == 1);
      TEST_ASSERT(neighborGlobalIds.size() == 1);
      TEST_ASSERT(neighborGlobalIds[0] == 1);
    }
    else if(nodeGlobalId == 1){
      TEST_ASSERT(numNeighbors == 1);
      TEST_ASSERT(neighborGlobalIds.size() == 1);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
    }
  }

}

TEUCHOS_UNIT_TEST(ProximitySearch, FivePointProblem) {

  Teuchos::RCP<Epetra_Comm> comm;
#ifdef HAVE_MPI
  comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
  int numProc = comm->NumProc();

  // This test cannot be run on more than 5 processors
  if(numProc > 5){
    std::cerr << "Unit test runtime ERROR: utPeridigm_Compute_Force only makes sense on 1 to 5 processors." << std::endl;
    return;
  }


  Epetra_BlockMap map(5, 3, 0, *comm);
  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(map));

  int numMyElements = map.NumMyElements();

  if(numProc == 1)
    TEST_ASSERT(numMyElements == 5);
  if(numProc == 5)
    TEST_ASSERT(numMyElements == 1);

  vector<double> node(3);
  std::map<int, vector<double> > nodes;
  node[0] = 0.0 ; node[1] = 0.0 ; node[2] = 0.0 ; nodes[0] = node;
  node[0] = 1.0 ; node[1] = 0.0 ; node[2] = 0.0 ; nodes[1] = node;
  node[0] = 0.0 ; node[1] = 1.0 ; node[2] = 0.0 ; nodes[2] = node;
  node[0] = 0.0 ; node[1] = 0.0 ; node[2] = 1.0 ; nodes[3] = node;
  node[0] = 1.0 ; node[1] = 1.0 ; node[2] = 1.0 ; nodes[4] = node;

  for(int i=0 ; i<numMyElements ; ++i){
    int globalId = map.GID(i);
    (*x)[3*i]   = nodes[globalId][0];
    (*x)[3*i+1] = nodes[globalId][1];
    (*x)[3*i+2] = nodes[globalId][2];
  }

  // These are filled by the proximity search
  Teuchos::RCP<Epetra_BlockMap> overlapMap;
  int neighborListSize(0);
  int* neighborList(0);

  // ---- Call the proximity search with a radius of 1.1 ----

  double searchRadius = 1.1;
  Epetra_BlockMap oneDimensionalMap(map.NumGlobalElements(), map.NumMyElements(), map.MyGlobalElements(), 1, 0, *comm);
  Teuchos::RCP<Epetra_Vector> searchRadii = Teuchos::rcp(new Epetra_Vector(oneDimensionalMap));
  searchRadii->PutScalar(searchRadius);

  ProximitySearch::GlobalProximitySearch(x, searchRadii, overlapMap, neighborListSize, neighborList);

  int neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    TEST_ASSERT(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      TEST_ASSERT(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      TEST_ASSERT(numNeighbors == 3);
      TEST_ASSERT(neighborGlobalIds.size() == 3);
      TEST_ASSERT(neighborGlobalIds[0] == 1);
      TEST_ASSERT(neighborGlobalIds[1] == 2);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
    }
    else if(nodeGlobalId == 1){
      TEST_ASSERT(numNeighbors == 1);
      TEST_ASSERT(neighborGlobalIds.size() == 1);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
    }
    else if(nodeGlobalId == 2){
      TEST_ASSERT(numNeighbors == 1);
      TEST_ASSERT(neighborGlobalIds.size() == 1);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
    }
    else if(nodeGlobalId == 3){
      TEST_ASSERT(numNeighbors == 1);
      TEST_ASSERT(neighborGlobalIds.size() == 1);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
    }
    else if(nodeGlobalId == 4){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
  }

  // ---- Call the proximity search with a radius of 0.0001 ----

  searchRadius = 0.0001;
  searchRadii->PutScalar(searchRadius);

  ProximitySearch::GlobalProximitySearch(x, searchRadii, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    TEST_ASSERT(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      TEST_ASSERT(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 1){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 2){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 3){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 4){
      TEST_ASSERT(numNeighbors == 0);
      TEST_ASSERT(neighborGlobalIds.size() == 0);
    }
  }

  // ---- Call the proximity search with a radius of 1.415 ----

  searchRadius = 1.415;
  searchRadii->PutScalar(searchRadius);
  ProximitySearch::GlobalProximitySearch(x, searchRadii, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    TEST_ASSERT(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      TEST_ASSERT(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      TEST_ASSERT(numNeighbors == 3);
      TEST_ASSERT(neighborGlobalIds.size() == 3);
      TEST_ASSERT(neighborGlobalIds[0] == 1);
      TEST_ASSERT(neighborGlobalIds[1] == 2);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
    }
    else if(nodeGlobalId == 1){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 2);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 2){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 1);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 3){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 1);
      TEST_ASSERT(neighborGlobalIds[2] == 2);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 4){
      TEST_ASSERT(numNeighbors == 3);
      TEST_ASSERT(neighborGlobalIds.size() == 3);
      TEST_ASSERT(neighborGlobalIds[0] == 1);
      TEST_ASSERT(neighborGlobalIds[1] == 2);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
    }
  }


  // ---- Call the proximity search with a radius of 1.733 ----

  searchRadius = 1.733;
  searchRadii->PutScalar(searchRadius);

  ProximitySearch::GlobalProximitySearch(x, searchRadii, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    TEST_ASSERT(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      TEST_ASSERT(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 1);
      TEST_ASSERT(neighborGlobalIds[1] == 2);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 1){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 2);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 2){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 1);
      TEST_ASSERT(neighborGlobalIds[2] == 3);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 3){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 1);
      TEST_ASSERT(neighborGlobalIds[2] == 2);
      TEST_ASSERT(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 4){
      TEST_ASSERT(numNeighbors == 4);
      TEST_ASSERT(neighborGlobalIds.size() == 4);
      TEST_ASSERT(neighborGlobalIds[0] == 0);
      TEST_ASSERT(neighborGlobalIds[1] == 1);
      TEST_ASSERT(neighborGlobalIds[2] == 2);
      TEST_ASSERT(neighborGlobalIds[3] == 3);
    }
  }

}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
