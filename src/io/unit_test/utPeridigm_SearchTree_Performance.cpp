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

#define PERFORMANCE_TESTS 0

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

//! Read a mesh from a text file
void readMeshFromTextFile(vector<double>& mesh, string fileName)
{
  ifstream inFile(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  BOOST_CHECK_EQUAL(inFile.is_open(), true);
  while(inFile.good()){
    string str;
    getline(inFile, str);
    // Ignore comment lines, otherwise parse
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
      istringstream iss(str);
      vector<double> data;
      copy(istream_iterator<double>(iss),
           istream_iterator<double>(),
           back_inserter<vector<double> >(data));
      // Check for obvious problems with the data
      BOOST_CHECK_EQUAL(static_cast<int>(data.size()), 5);
      // Store the coordinates
      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);
    }
  }
  inFile.close();
}

//! Tests the search tree associated with the equally-spaced 1000-point cube mesh
void testEquallySpacedCubeMesh1000(vector<double>& mesh, PeridigmNS::SearchTree* searchTree)
{
  double* meshPtr = &mesh[0];
  vector<int> neighborList;
  int searchPointIndex, degreesOfFreedom(3);
  double searchRadius;

  // This search should find three other points
  neighborList.clear();
  searchPointIndex = 0;
  searchRadius = 1.015;
  searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
  BOOST_CHECK_EQUAL(static_cast<int>(neighborList.size()), 4);
}

//! Zoltan 1000-point test
void testZoltanEquallySpacedCubeMesh1000()
{
  vector<double> mesh;
  string fileName("./input_files/cube_1000.txt");
  readMeshFromTextFile(mesh, fileName);
  PeridigmNS::SearchTree* searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);
  testEquallySpacedCubeMesh1000(mesh, searchTree);
  delete searchTree;
}

//! JAM 1000-point test
void testJAMEquallySpacedCubeMesh1000()
{
  vector<double> mesh;
  string fileName("./input_files/cube_1000.txt");
  readMeshFromTextFile(mesh, fileName);
  PeridigmNS::SearchTree* searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), &mesh[0]);
  testEquallySpacedCubeMesh1000(mesh, searchTree);
  delete searchTree;
}

//! Zoltan performance test
void testZoltanPerformance()
{
  vector<int> neighborList;
  int searchPointIndex, degreesOfFreedom(3);
  double searchRadius;
  vector<double> mesh;
  string fileName;
  double* meshPtr;
  PeridigmNS::SearchTree* searchTree;
  unsigned int numBonds, totalBonds, maxBonds;

  // Create a 1000-point discretization and find the neighbors of all the points
  mesh.clear();
  fileName = "./input_files/cube_1000.txt";
  readMeshFromTextFile(mesh, fileName);
  meshPtr = &mesh[0];
  searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), meshPtr);
  totalBonds = 0;
  maxBonds = 0;
  neighborList.resize(130);
  for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    searchRadius = 1.0*3.015;
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
  }
  BOOST_CHECK_EQUAL(totalBonds, static_cast<unsigned int>(84288));
  BOOST_CHECK_EQUAL(maxBonds, static_cast<unsigned int>(122));
  delete searchTree;

  // Create a 8000-point discretization and find the neighbors of all the points
  mesh.clear();
  fileName = "./input_files/cube_8000.txt";
  readMeshFromTextFile(mesh, fileName);
  meshPtr = &mesh[0];
  searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), meshPtr);
  totalBonds = 0;
  maxBonds = 0;
  neighborList.resize(130);
  for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    searchRadius = 0.5*3.015;
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
  }
  BOOST_CHECK_EQUAL(totalBonds, static_cast<unsigned int>(816728));
  BOOST_CHECK_EQUAL(maxBonds, static_cast<unsigned int>(122));
  delete searchTree;

  // Create a 27000-point discretization and find the neighbors of all the points
  mesh.clear();
  fileName = "./input_files/cube_27000.txt";
  readMeshFromTextFile(mesh, fileName);
  meshPtr = &mesh[0];
  searchTree = new PeridigmNS::ZoltanSearchTree(static_cast<int>(mesh.size()/3), meshPtr);
  totalBonds = 0;
  maxBonds = 0;
  neighborList.resize(130);
  for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    searchRadius = (1.0/3.0)*3.015;
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
  }
  BOOST_CHECK_EQUAL(totalBonds, static_cast<unsigned int>(2929168));
  BOOST_CHECK_EQUAL(maxBonds, static_cast<unsigned int>(122));
  delete searchTree;
}

//! JAM performance test
void testJAMPerformance()
{
  vector<int> neighborList;
  int searchPointIndex, degreesOfFreedom(3);
  double searchRadius;
  vector<double> mesh;
  string fileName;
  double* meshPtr;
  PeridigmNS::SearchTree* searchTree;
  unsigned int numBonds, totalBonds, maxBonds;

  // Create a 1000-point discretization and find the neighbors of all the points
  mesh.clear();
  fileName = "./input_files/cube_1000.txt";
  readMeshFromTextFile(mesh, fileName);
  meshPtr = &mesh[0];
  searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), meshPtr);
  totalBonds = 0;
  maxBonds = 0;
  neighborList.resize(130);
  for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    searchRadius = 1.0*3.015;
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
  }
  BOOST_CHECK_EQUAL(totalBonds, static_cast<unsigned int>(84288));
  BOOST_CHECK_EQUAL(maxBonds, static_cast<unsigned int>(122));
  delete searchTree;

  // Create a 8000-point discretization and find the neighbors of all the points
  mesh.clear();
  fileName = "./input_files/cube_8000.txt";
  readMeshFromTextFile(mesh, fileName);
  meshPtr = &mesh[0];
  searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), meshPtr);
  totalBonds = 0;
  maxBonds = 0;
  neighborList.resize(130);
  for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    searchRadius = 0.5*3.015;
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
  }
  BOOST_CHECK_EQUAL(totalBonds, static_cast<unsigned int>(816728));
  BOOST_CHECK_EQUAL(maxBonds, static_cast<unsigned int>(122));
  delete searchTree;

  // Create a 27000-point discretization and find the neighbors of all the points
  mesh.clear();
  fileName = "./input_files/cube_27000.txt";
  readMeshFromTextFile(mesh, fileName);
  meshPtr = &mesh[0];
  searchTree = new PeridigmNS::JAMSearchTree(static_cast<int>(mesh.size()/3), meshPtr);
  totalBonds = 0;
  maxBonds = 0;
  neighborList.resize(130);
  for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    searchRadius = (1.0/3.0)*3.015;
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
  }
  BOOST_CHECK_EQUAL(totalBonds, static_cast<unsigned int>(2929168));
  BOOST_CHECK_EQUAL(maxBonds, static_cast<unsigned int>(122));
  delete searchTree;
}

bool init_unit_test_suite()
{
  // Add a suite for each processor in the test
  bool success = true;

  test_suite* proc = BOOST_TEST_SUITE("utPeridigm_SearchTree");
  proc->add(BOOST_TEST_CASE(&testJAMPerformance));
  proc->add(BOOST_TEST_CASE(&testZoltanPerformance));
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

    // Generate the test meshes
    system("cd input_files ; python CreateTestMeshes.py");

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
