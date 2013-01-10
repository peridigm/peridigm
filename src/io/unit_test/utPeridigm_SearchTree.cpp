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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/unit_test.hpp>
#include "Peridigm_JAMSearchTree.hpp"
#include "Peridigm_ZoltanSearchTree.hpp"
#include <Epetra_SerialComm.h>
#include <vector>
#include <sstream>
#include <fstream>

using namespace boost::unit_test;
using namespace std;
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
void testEightPointMesh(vector<double>& mesh, PeridigmNS::SearchTree* searchTree)
{
  double* meshPtr = &mesh[0];
  vector<int> neighborList;
  int searchPointIndex, degreesOfFreedom(3);
  double searchRadius;

  // This search should find all the other points
  neighborList.clear();
  searchPointIndex = 2;
  searchRadius = 3.015;
  searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
  BOOST_CHECK_EQUAL(static_cast<int>(neighborList.size()), 8);
  sort(neighborList.begin(), neighborList.end());
  for(int i=0 ; i<8 ; ++i)
    BOOST_CHECK_EQUAL(neighborList[i], i);

  // This search should find thee neighbors
  neighborList.clear();
  searchPointIndex = 0;
  searchRadius = 1.015;
  searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
  BOOST_CHECK_EQUAL(static_cast<int>(neighborList.size()), 4);
  sort(neighborList.begin(), neighborList.end());
  BOOST_CHECK_EQUAL(neighborList[0], 0);
  BOOST_CHECK_EQUAL(neighborList[1], 1);
  BOOST_CHECK_EQUAL(neighborList[2], 2);
  BOOST_CHECK_EQUAL(neighborList[3], 4);

  // This search should find no neighbors
  neighborList.clear();
  searchPointIndex = 0;
  searchRadius = 0.015;
  searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
  BOOST_CHECK_EQUAL(static_cast<int>(neighborList.size()), 1);
}

//! Zoltan eight-point test
void testZoltanEightPointMesh()
{
  vector<double> mesh;
  eightPointMesh(mesh);
  PeridigmNS::SearchTree* searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);
  testEightPointMesh(mesh, searchTree);
  delete searchTree;
}

//! JAM eight-point test
void testJAMEightPointMesh()
{
  vector<double> mesh;
  eightPointMesh(mesh);
  PeridigmNS::SearchTree* searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);
  testEightPointMesh(mesh, searchTree);
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
//   BOOST_CHECK_EQUAL(static_cast<int>(neighborList.size()), 4);
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

bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPeridigm_SearchTree");
  proc->add(BOOST_TEST_CASE(&testZoltanEightPointMesh));
  proc->add(BOOST_TEST_CASE(&testJAMEightPointMesh));
  // proc->add(BOOST_TEST_CASE(&testZoltanEquallySpacedCubeMesh1000));
  // proc->add(BOOST_TEST_CASE(&testJAMEquallySpacedCubeMesh1000));
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
  // This is a serial test, but it requires MPI to be running (Zoltan dependency?)

  int numProcs = 1;
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

  int returnCode = -1;
  if(numProcs == 1){
    // Run the tests
    returnCode = unit_test_main(init_unit_test, argc, argv);
  }
  else{
    std::cerr << "Unit test runtime ERROR: utPeridigm_State only makes sense on 1 processor." << std::endl;
  }
  
#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  
  return returnCode;
}
