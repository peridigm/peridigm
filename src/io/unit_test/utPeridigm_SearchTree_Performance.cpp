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


#include "Peridigm_Timer.hpp"
#include "Peridigm_JAMSearchTree.hpp"
#include "Peridigm_ZoltanSearchTree.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include <Epetra_SerialComm.h>
#include <vector>
#include <sstream>
#include <fstream>


using namespace std;
using namespace PeridigmNS;




PeridigmNS::SearchTree* createTree(string treeType, int numPoints, double* coordinates)
{
  PeridigmNS::SearchTree* tree(NULL);
  if(treeType == "Zoltan")
    tree = new PeridigmNS::ZoltanSearchTree(numPoints, coordinates);
  else if(treeType == "JAM")
    tree = new PeridigmNS::JAMSearchTree(numPoints, coordinates);
  return tree;
}


//! Performance tests

void testPerformance( vector<int> neighborList, double searchRadius, vector<double> mesh, string testName, string treeType, PeridigmNS::SearchTree* searchTree, unsigned int &totalBonds, unsigned int &maxBonds, unsigned int &minBonds){

   double* meshPtr = &mesh[0];
   totalBonds = 0;
   maxBonds = 0;
   minBonds = 1000000;

   int searchPointIndex;
   int degreesOfFreedom(3);
   //neighborList.clear();
   
   searchTree = createTree(treeType, static_cast<int>(mesh.size()/3), meshPtr);
   neighborList.resize(130);

   for(unsigned int i=0 ; i<mesh.size()/3 ; i++){
    neighborList.clear();
    searchPointIndex = i;
    
    searchTree->FindPointsWithinRadius(&meshPtr[searchPointIndex*degreesOfFreedom], searchRadius, neighborList);
    unsigned int numBonds = neighborList.size() - 1;
    totalBonds += numBonds;
    if(numBonds > maxBonds)
      maxBonds = numBonds;
    if(numBonds < minBonds)
      minBonds = numBonds;
  }
}



//! Performance tests


TEUCHOS_UNIT_TEST(SearchTree_Performance, ZoltanTest) {

  vector<int> neighborList;
  double searchRadius;
  vector<double> mesh;
  string fileName, testName, treeType;
  PeridigmNS::SearchTree* searchTree;
  unsigned int totalBonds, maxBonds, minBonds;
  string str;
  vector<double> data;
  double num;
  ifstream inFile;
  treeType = "Zoltan";

  
// Create a 8022-point discretization shaped like a dumbbell and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/dumbbell.txt";
  searchRadius = (1.0/3.0)*3.015;

  //! Read a mesh from a text file

  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);

  while(inFile.good()){
   
       getline(inFile, str);
       // Ignore comment lines, otherwise parse

       if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
           istringstream iss(str);
 
           while ( iss >> num) data.push_back(num);

          // Check for obvious problems with the data

          TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

         // Store the coordinates
          mesh.push_back(data[0]);
          mesh.push_back(data[1]);
          mesh.push_back(data[2]);

        data.clear();
    }
  }
  inFile.close();

  testName = treeType + " test 1)  Dumbbell mesh with 8022 points";
  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList,  searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(8630086));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(1934));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(52));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);


 
  // Create a random, 8000-point discretization and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/random.txt";
  //! Read a mesh from a text file
 
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
    TEST_EQUALITY(inFile.is_open(), true);

    while(inFile.good()){
    
          getline(inFile, str);
    
          if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
              istringstream iss(str);
      
              while ( iss >> num) data.push_back(num);

              TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

              mesh.push_back(data[0]);
              mesh.push_back(data[1]);
              mesh.push_back(data[2]);

              data.clear();
    }
  }
  inFile.close();

  testName = treeType + " test 2)  Random mesh with 8000 points";
  searchRadius = 3.0;
  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList,  searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(5005818));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(963));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(127));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);



  // Create a 27000-point discretization and find the neighbors of all the points

 
  mesh.clear();
  fileName = "./input_files/cube_27000.txt";

  //! Read a mesh from a text file
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
    getline(inFile, str);
    
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
        istringstream iss(str);
      
        while ( iss >> num) data.push_back(num);

        TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

        mesh.push_back(data[0]);
        mesh.push_back(data[1]);
        mesh.push_back(data[2]);

        data.clear();
    }
  }
  inFile.close();

  
  testName = treeType + " test 3)  Equally-Spaced Cube with 27000 points";
  searchRadius = (1.0/3.0)*3.015;
  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList, searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);
  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(2929168));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(122));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(28));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);


  // Create a 8000-point discretization and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/cube_8000.txt";

  //! Read a mesh from a text file
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
    getline(inFile, str);
    
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
       istringstream iss(str);
      
      while ( iss >> num) data.push_back(num);

      TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);

      data.clear();
    }
  }
  inFile.close();

  testName = treeType + " test 4)  Equally-Spaced Cube with 8000 points";
  searchRadius = 0.5*3.015;
  PeridigmNS::Timer::self().startTimer(testName);

  testPerformance( neighborList, searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(816728));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(122));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(28));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);


  // Create a 1000-point discretization and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/cube_1000.txt";

  //! Read a mesh from a text file
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
        getline(inFile, str);
    
        if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
            istringstream iss(str);
      
        while ( iss >> num) data.push_back(num);

        TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

        mesh.push_back(data[0]);
        mesh.push_back(data[1]);
        mesh.push_back(data[2]);

        data.clear();
    }
  }
  inFile.close();

  testName = treeType + " test 5)  Equally-Spaced Cube with 1000 points";
  searchRadius = 1.0*3.015;

  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList, searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(84288));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(122));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(28));
  delete searchTree;
  
  PeridigmNS::Timer::self().stopTimer(testName);

  
}

TEUCHOS_UNIT_TEST(SearchTree_Performance, JAMTest) {

  vector<int> neighborList;
  double searchRadius;
  vector<double> mesh;
  string fileName, testName, treeType;
  PeridigmNS::SearchTree* searchTree;
  unsigned int totalBonds, maxBonds, minBonds;
  string str;
  vector<double> data;
  double num;
  ifstream inFile;


  treeType = "JAM";

  // Create a 8022-point discretization shaped like a dumbbell and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/dumbbell.txt";
  
  searchRadius = (1.0/3.0)*3.015;
  //! Read a mesh from a text file

  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
   
    getline(inFile, str);
    // Ignore comment lines, otherwise parse
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
      
      istringstream iss(str);
      

      while ( iss >> num) data.push_back(num);

      // Check for obvious problems with the data

      TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

      // Store the coordinates
      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);

      data.clear();

      
    }
  }
  inFile.close();

  testName = treeType + " test 1)  Dumbbell mesh with 8022 points";

  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList,  searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(8630086));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(1934));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(52));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);

  // Create a random, 8000-point discretization and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/random.txt";
  //! Read a mesh from a text file
 
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
    getline(inFile, str);
    
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
       istringstream iss(str);
      
      while ( iss >> num) data.push_back(num);

      TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);

      data.clear();

      
    }
  }
  inFile.close();

  testName = treeType + " test 2)  Random mesh with 8000 points";
  searchRadius = 3.0;

  PeridigmNS::Timer::self().startTimer(testName);
 
  testPerformance( neighborList,  searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(5005818));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(963));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(127));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);

  // Create a 27000-point discretization and find the neighbors of all the points

 
  mesh.clear();
  fileName = "./input_files/cube_27000.txt";

  //! Read a mesh from a text file
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
    getline(inFile, str);
    
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
       istringstream iss(str);
      
      while ( iss >> num) data.push_back(num);

      TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);

      data.clear();

      
    }
  }
  inFile.close();

  
  testName = treeType + " test 3)  Equally-Spaced Cube with 27000 points";
  searchRadius = (1.0/3.0)*3.015;

  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList, searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);
  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(2929168));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(122));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(28));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);

  // Create a 8000-point discretization and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/cube_8000.txt";

  //! Read a mesh from a text file
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
    getline(inFile, str);
    
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
       istringstream iss(str);
      
      while ( iss >> num) data.push_back(num);

      TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);

      data.clear();
    }
  }
  inFile.close();

  testName = treeType + " test 4)  Equally-Spaced Cube with 8000 points";
  searchRadius = 0.5*3.015;

  PeridigmNS::Timer::self().startTimer(testName);

  testPerformance( neighborList, searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(816728));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(122));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(28));
  delete searchTree;
  PeridigmNS::Timer::self().stopTimer(testName);

  // Create a 1000-point discretization and find the neighbors of all the points

  mesh.clear();
  fileName = "./input_files/cube_1000.txt";

  //! Read a mesh from a text file
  inFile.open(fileName.c_str());
  if(!inFile.is_open())
    cout << "\n**** Warning:  This test can only be run from the directory where it resides (otherwise it won't find the input files) ****\n" << endl;
  TEST_EQUALITY(inFile.is_open(), true);
  while(inFile.good()){
    
    getline(inFile, str);
    
    if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
       istringstream iss(str);
      
      while ( iss >> num) data.push_back(num);

      TEST_EQUALITY_CONST(static_cast<int>(data.size()), 5);

      mesh.push_back(data[0]);
      mesh.push_back(data[1]);
      mesh.push_back(data[2]);

      data.clear();
    }
  }
  inFile.close();

  testName = treeType + " test 5)  Equally-Spaced Cube with 1000 points";
  searchRadius = 1.0*3.015;

  PeridigmNS::Timer::self().startTimer(testName);
  testPerformance( neighborList, searchRadius, mesh, testName, treeType, searchTree, totalBonds, maxBonds, minBonds);

  TEST_EQUALITY_CONST(totalBonds, static_cast<unsigned int>(84288));
  TEST_EQUALITY_CONST(maxBonds, static_cast<unsigned int>(122));
  TEST_EQUALITY_CONST(minBonds, static_cast<unsigned int>(28));
  delete searchTree;
  
  PeridigmNS::Timer::self().stopTimer(testName);

}



int main
(int argc, char* argv[])
{
  // This is a serial test, but it requires MPI to be running (Zoltan dependency?)

  int numProcs = 1;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);


  int returnCode = -1;
  if(numProcs == 1){

    // Generate the test meshes
    system("cd input_files ; python CreateTestMeshes.py");

    // Run the tests
    returnCode = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  }
  else{
    std::cerr << "Unit test runtime ERROR: utPeridigm_State only makes sense on 1 processor." << std::endl;
  }

  std::cout << endl;
  PeridigmNS::Timer::self().printTimingData(cout);
  

  return returnCode;
}
