/*! \file utPeridigm_SearchTree.cpp */

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

#include "Peridigm_JAMSearchTree.hpp"
#include "Peridigm_ZoltanSearchTree.hpp"
#include <Epetra_SerialComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <vector>
#include <sstream>
#include <fstream>

using std::vector;
using std::string;
using namespace PeridigmNS;

//! Generate simple eight-point cube mesh
void eightPointMesh(vector<double>& mesh)
{
  mesh = vector<double>(8*3);
  mesh[0]  = 0.0 ; mesh[1]  = 0.0 ; mesh[2]  = 0.0;
  mesh[3]  = 1.0 ; mesh[4]  = 0.0 ; mesh[5]  = 0.0;
  mesh[6]  = 0.0 ; mesh[7]  = 1.0 ; mesh[8]  = 0.0;
  mesh[9]  = 1.0 ; mesh[10] = 1.0 ; mesh[11] = 0.0;
  mesh[12] = 0.0 ; mesh[13] = 0.0 ; mesh[14] = 1.0;
  mesh[15] = 1.0 ; mesh[16] = 0.0 ; mesh[17] = 1.0;
  mesh[18] = 0.0 ; mesh[19] = 1.0 ; mesh[20] = 1.0;
  mesh[21] = 1.0 ; mesh[22] = 1.0 ; mesh[23] = 1.0;
}


//! Test a search tree created for the eight-point cube mesh

void testEightPointMesh(vector<double>& mesh, PeridigmNS::SearchTree* searchTree, vector<int>& neighborList, int searchPointIndex, int degreesOfFreedom, double searchRadius){

 double* meshPtr = &mesh[0];
 neighborList.clear();
 searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
 if (neighborList.size() > 1) 
    sort(neighborList.begin(), neighborList.end());

}


//! Zoltan eight-point test

TEUCHOS_UNIT_TEST(SearchTree, ZoltanEightPointMesh) {

  vector<double> mesh;
  eightPointMesh(mesh);
  vector<int> neighborList;
  int searchPointIndex, degreesOfFreedom(3);
  double searchRadius;
  PeridigmNS::SearchTree* searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);

  // This search should find all the other points
  
  searchPointIndex = 2;
  searchRadius = 3.015;
  testEightPointMesh(mesh,searchTree, neighborList, searchPointIndex, degreesOfFreedom, searchRadius);
  TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 8);

  for(int i=0 ; i<4 ; ++i)
    TEST_EQUALITY(neighborList[i], i);
  
  

 // This search should find three neighbors
  
  searchPointIndex = 0;
  searchRadius = 1.015;
  testEightPointMesh(mesh,searchTree, neighborList, searchPointIndex, degreesOfFreedom, searchRadius);
  TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 4);

  TEST_EQUALITY_CONST(neighborList[0], 0);
  TEST_EQUALITY_CONST(neighborList[1], 1);
  TEST_EQUALITY_CONST(neighborList[2], 2);
  TEST_EQUALITY_CONST(neighborList[3], 4);
   
  
  
 // This search should find no neighbors
 
  searchPointIndex = 0;
  searchRadius = 0.015;
  testEightPointMesh(mesh,searchTree, neighborList, searchPointIndex, degreesOfFreedom, searchRadius);
  TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 1);

  delete searchTree;
}

//! JAM eight-point test

TEUCHOS_UNIT_TEST(SearchTree, JAMEightPointMesh) {

  vector<double> mesh;
  eightPointMesh(mesh);

  vector<int> neighborList;
  int searchPointIndex, degreesOfFreedom(3);
  double searchRadius;
  PeridigmNS::SearchTree* searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);

  // This search should find all the other points
  
  searchPointIndex = 2;
  searchRadius = 3.015;
  testEightPointMesh(mesh,searchTree, neighborList, searchPointIndex, degreesOfFreedom, searchRadius);
  TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 8);
  
  for(int i=0 ; i<8 ; ++i)
    TEST_EQUALITY(neighborList[i], i);

 // This search should find three neighbors
  
  searchPointIndex = 0;
  searchRadius = 1.015;
  testEightPointMesh(mesh,searchTree, neighborList, searchPointIndex, degreesOfFreedom, searchRadius);
  TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 4);
   
  TEST_EQUALITY_CONST(neighborList[0], 0);
  TEST_EQUALITY_CONST(neighborList[1], 1);
  TEST_EQUALITY_CONST(neighborList[2], 2);
  TEST_EQUALITY_CONST(neighborList[3], 4);

  
 // This search should find no neighbors
 
  searchPointIndex = 0;
  searchRadius = 0.015;
  testEightPointMesh(mesh,searchTree, neighborList, searchPointIndex, degreesOfFreedom, searchRadius);
  TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 1);
 
  delete searchTree;
}



// //! Tests the search tree associated with the equally-spaced 1000-point cube mesh
// void testEquallySpacedCubeMesh1000(vector<double>& mesh, PeridigmNS::SearchTree* searchTree)
// {
//   double* meshPtr = &mesh[0];
//   vector<int> neighborList;
//   int searchPointIndex, degreesOfFreedom(3);
//   double searchRadius;

//   // This search should find three other points
//   neighborList.clear();
//   searchPointIndex = 0;
//   searchRadius = 1.015;
//   searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
//   TEST_EQUALITY_CONST(static_cast<int>(neighborList.size()), 4);
// }

// //! Zoltan 1000-point test
// void testZoltanEquallySpacedCubeMesh1000()
// {
//   vector<double> mesh;
//   string fileName("./input_files/cube_1000.txt");
//   readMeshFromTextFile(mesh, fileName);
//   PeridigmNS::SearchTree* searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);
//   testEquallySpacedCubeMesh1000(mesh, searchTree);
//   delete searchTree;
// }

// //! JAM 1000-point test
// void testJAMEquallySpacedCubeMesh1000()
// {
//   vector<double> mesh;
//   string fileName("./input_files/cube_1000.txt");
//   readMeshFromTextFile(mesh, fileName);
//   PeridigmNS::SearchTree* searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);
//   testEquallySpacedCubeMesh1000(mesh, searchTree);
//   delete searchTree;
// }



int main
(int argc, char* argv[])
{
  // This is a serial test, but it requires MPI to be running (Zoltan dependency?)

  int numProcs = 1;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int returnCode = -1;
  if(numProcs == 1){
    // Run the tests
    
    returnCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  }
  else{
    std::cerr << "Unit test runtime ERROR: utPeridigm_State only makes sense on 1 processor." << std::endl;
  }
  
  return returnCode;
}
