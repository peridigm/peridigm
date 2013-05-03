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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>
#include <Epetra_ConfigDefs.h> // used to define HAVE_MPI
#ifdef HAVE_MPI
  #include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>
#include "Peridigm_ProximitySearch.hpp"
#include <vector>

using namespace boost::unit_test;
using namespace Teuchos;
using namespace PeridigmNS;
using namespace std;

void fivePointProblem()
{
  Teuchos::RCP<Epetra_Comm> comm;
#ifdef HAVE_MPI
  comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  comm = Teuchos::rcp(new Epetra_SerialComm);
#endif
  int numProc = comm->NumProc();
  // int myPID = comm->MyPID();

  Epetra_BlockMap map(5, 3, 0, *comm);
  Epetra_Vector x(map);

  int numMyElements = map.NumMyElements();

  if(numProc == 1)
    BOOST_CHECK(numMyElements == 5);
  if(numProc == 5)
    BOOST_CHECK(numMyElements == 1);

  vector<int> node(3);
  std::map<int, vector<int> > nodes;
  node[0] = 0.0 ; node[1] = 0.0 ; node[2] = 0.0 ; nodes[0] = node;
  node[0] = 1.0 ; node[1] = 0.0 ; node[2] = 0.0 ; nodes[1] = node;
  node[0] = 0.0 ; node[1] = 1.0 ; node[2] = 0.0 ; nodes[2] = node;
  node[0] = 0.0 ; node[1] = 0.0 ; node[2] = 1.0 ; nodes[3] = node;
  node[0] = 1.0 ; node[1] = 1.0 ; node[2] = 1.0 ; nodes[4] = node;

  for(int i=0 ; i<numMyElements ; ++i){
    int globalId = map.GID(i);
    x[3*i]   = nodes[globalId][0];
    x[3*i+1] = nodes[globalId][1];
    x[3*i+2] = nodes[globalId][2];
  }

  // These are filled by the proximity search
  Teuchos::RCP<Epetra_BlockMap> overlapMap;
  int neighborListSize(0);
  int* neighborList(0);

  // ---- Call the proximity search with a radius of 1.1 ----

  double searchRadius = 1.1;
  ProximitySearch::GlobalProximitySearch(x, searchRadius, overlapMap, neighborListSize, neighborList);

  int neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    BOOST_CHECK(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      BOOST_CHECK(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      BOOST_CHECK(numNeighbors == 3);
      BOOST_CHECK(neighborGlobalIds.size() == 3);
      BOOST_CHECK(neighborGlobalIds[0] == 1);
      BOOST_CHECK(neighborGlobalIds[1] == 2);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
    }
    else if(nodeGlobalId == 1){
      BOOST_CHECK(numNeighbors == 1);
      BOOST_CHECK(neighborGlobalIds.size() == 1);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
    }
    else if(nodeGlobalId == 2){
      BOOST_CHECK(numNeighbors == 1);
      BOOST_CHECK(neighborGlobalIds.size() == 1);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
    }
    else if(nodeGlobalId == 3){
      BOOST_CHECK(numNeighbors == 1);
      BOOST_CHECK(neighborGlobalIds.size() == 1);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
    }
    else if(nodeGlobalId == 4){
      BOOST_CHECK(numNeighbors == 0);
      BOOST_CHECK(neighborGlobalIds.size() == 0);
    }
  }

  // ---- Call the proximity search with a radius of 0.0001 ----

  searchRadius = 0.0001;
  ProximitySearch::GlobalProximitySearch(x, searchRadius, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    BOOST_CHECK(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      BOOST_CHECK(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      BOOST_CHECK(numNeighbors == 0);
      BOOST_CHECK(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 1){
      BOOST_CHECK(numNeighbors == 0);
      BOOST_CHECK(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 2){
      BOOST_CHECK(numNeighbors == 0);
      BOOST_CHECK(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 3){
      BOOST_CHECK(numNeighbors == 0);
      BOOST_CHECK(neighborGlobalIds.size() == 0);
    }
    else if(nodeGlobalId == 4){
      BOOST_CHECK(numNeighbors == 0);
      BOOST_CHECK(neighborGlobalIds.size() == 0);
    }
  }

  // ---- Call the proximity search with a radius of 1.415 ----

  searchRadius = 1.415;
  ProximitySearch::GlobalProximitySearch(x, searchRadius, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    BOOST_CHECK(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      BOOST_CHECK(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      BOOST_CHECK(numNeighbors == 3);
      BOOST_CHECK(neighborGlobalIds.size() == 3);
      BOOST_CHECK(neighborGlobalIds[0] == 1);
      BOOST_CHECK(neighborGlobalIds[1] == 2);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
    }
    else if(nodeGlobalId == 1){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 2);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 2){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 1);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 3){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 1);
      BOOST_CHECK(neighborGlobalIds[2] == 2);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 4){
      BOOST_CHECK(numNeighbors == 3);
      BOOST_CHECK(neighborGlobalIds.size() == 3);
      BOOST_CHECK(neighborGlobalIds[0] == 1);
      BOOST_CHECK(neighborGlobalIds[1] == 2);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
    }
  }


  // ---- Call the proximity search with a radius of 1.733 ----

  searchRadius = 1.733;
  ProximitySearch::GlobalProximitySearch(x, searchRadius, overlapMap, neighborListSize, neighborList);

  neighborListIndex = 0;
  for(int i=0 ; i<numMyElements ; ++i){
    int nodeGlobalId = map.GID(i);
    BOOST_CHECK(neighborListIndex < neighborListSize);
    int numNeighbors = neighborList[neighborListIndex++];
    vector<int> neighborGlobalIds;
    for(int j=0 ; j<numNeighbors ; ++j){
      BOOST_CHECK(neighborListIndex < neighborListSize);
      int neighborLocalId = neighborList[neighborListIndex++];
      int neighborGlobalId = overlapMap->GID(neighborLocalId);
      neighborGlobalIds.push_back(neighborGlobalId);
    }

    sort(neighborGlobalIds.begin(), neighborGlobalIds.end());

    if(nodeGlobalId == 0){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 1);
      BOOST_CHECK(neighborGlobalIds[1] == 2);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 1){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 2);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 2){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 1);
      BOOST_CHECK(neighborGlobalIds[2] == 3);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 3){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 1);
      BOOST_CHECK(neighborGlobalIds[2] == 2);
      BOOST_CHECK(neighborGlobalIds[3] == 4);
    }
    else if(nodeGlobalId == 4){
      BOOST_CHECK(numNeighbors == 4);
      BOOST_CHECK(neighborGlobalIds.size() == 4);
      BOOST_CHECK(neighborGlobalIds[0] == 0);
      BOOST_CHECK(neighborGlobalIds[1] == 1);
      BOOST_CHECK(neighborGlobalIds[2] == 2);
      BOOST_CHECK(neighborGlobalIds[3] == 3);
    }
  }

}

bool init_unit_test_suite()
{
	// Add a suite for each processor in the test
	bool success = true;

	test_suite* proc = BOOST_TEST_SUITE("utPeridigm_ProximitySearch");
	proc->add(BOOST_TEST_CASE(&fivePointProblem));
	framework::master_test_suite().add(proc);

	return success;
}

bool init_unit_test()
{
	init_unit_test_suite();
	return true;
}

int main
(int argc, char* argv[])
{
  int numProcs = 1;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

  int returnCode = -1;
  returnCode = unit_test_main(init_unit_test, argc, argv);
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return returnCode;
}
