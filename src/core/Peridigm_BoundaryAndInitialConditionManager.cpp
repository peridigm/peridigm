/*! \file Peridigm_BoundaryAndInitialConditionManager.cpp */

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

#include "Peridigm_BoundaryAndInitialConditionManager.hpp"

#include "muParser/muParserPeridigmFunctions.h"

using namespace std;

PeridigmNS::BoundaryAndInitialConditionManager::BoundaryAndInitialConditionManager(const Teuchos::ParameterList& boundaryAndInitialConditionParams)
  : params(boundaryAndInitialConditionParams),
    muParserX(0.0), muParserY(0.0), muParserZ(0.0), muParserT(0.0)
{
  // Set up muParser
  try {
    muParser.DefineVar("x", &muParserX);
    muParser.DefineVar("y", &muParserY);
    muParser.DefineVar("z", &muParserZ);
    muParser.DefineVar("t", &muParserT);
    muParser.DefineFun(_T("rnd"), mu::Rnd, false);
  } 
  catch (mu::Parser::exception_type &e)
    TEST_FOR_EXCEPT_MSG(1, e.GetMsg());
}

void PeridigmNS::BoundaryAndInitialConditionManager::initialize(Teuchos::RCP<AbstractDiscretization> discretization)
{
  nodeSets = Teuchos::rcp(new map< string, vector<int> >());

  // Load node sets defined in the input deck into the nodeSets container
  for(Teuchos::ParameterList::ConstIterator it = params.begin() ; it != params.end() ; it++){
	const string& name = it->first;
	size_t position = name.find("Node Set");
	if(position != string::npos){
	  stringstream ss(Teuchos::getValue<string>(it->second));
      TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
	  vector<int>& nodeList = (*nodeSets)[name];
	  int nodeID;
	  while(ss.good()){
		ss >> nodeID;
		nodeList.push_back(nodeID);
	  }
	}
  }

  // Load node sets defined in the mesh file into the nodeSets container
  Teuchos::RCP< map< string, vector<int> > > discretizationNodeSets = discretization->getNodeSets();
  for(map< string, vector<int> >::iterator it=discretizationNodeSets->begin() ; it!=discretizationNodeSets->end() ; it++){
    string name = it->first;
    TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
    vector<int>& nodeList = it->second;
    (*nodeSets)[name] = nodeList;
  }
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyInitialDisplacements(Teuchos::RCP<Epetra_Vector> x,
                                                                               Teuchos::RCP<Epetra_Vector> u,
                                                                               Teuchos::RCP<Epetra_Vector> y)
{
  const Epetra_BlockMap& threeDimensionalMap = x->Map();

  // apply the initial conditions
  Teuchos::ParameterList::ConstIterator it;
  for(it = params.begin() ; it != params.end() ; it++){
    const string & name = it->first;
    size_t position = name.find("Initial Displacement");
    if(position != string::npos){
      Teuchos::ParameterList & boundaryConditionParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
      string nodeSet = boundaryConditionParams.get<string>("Node Set");
      string type = boundaryConditionParams.get<string>("Type");
      string coordinate = boundaryConditionParams.get<string>("Coordinate");
      string function = boundaryConditionParams.get<string>("Value");

      int coord = 0;
      if(coordinate == "y" || coordinate == "Y")
        coord = 1;
      if(coordinate == "z" || coordinate == "Z")
        coord = 2;

      try{
        muParser.SetExpr(function);
      }
      catch (mu::Parser::exception_type &e)
        TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

      if (nodeSet == "All") { // Apply initial position to all locally-owned nodes
        for(int localNodeID = 0; localNodeID < x->MyLength(); localNodeID++) {
          muParserX = (*x)[localNodeID*3];
          muParserY = (*x)[localNodeID*3 + 1];
          muParserZ = (*x)[localNodeID*3 + 2];
          try {
            (*u)[localNodeID*3 + coord] = muParser.Eval();
          }
          catch (mu::Parser::exception_type &e)
          TEST_FOR_EXCEPT_MSG(1, e.GetMsg());
        }
      }
      else { // Apply initial position to specific node set
        // apply initial displacement boundary conditions to locally-owned nodes
        TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSet) == nodeSets->end(), "**** Node set not found: " + name + "\n");
        vector<int> & nodeList = (*nodeSets)[nodeSet];
        for(unsigned int i=0 ; i<nodeList.size() ; i++){
          int localNodeID = threeDimensionalMap.LID(nodeList[i]);
          if(localNodeID != -1) {
            muParserX = (*x)[localNodeID*3];
            muParserY = (*x)[localNodeID*3 + 1];
            muParserZ = (*x)[localNodeID*3 + 2];
            try {
              (*u)[localNodeID*3 + coord] = muParser.Eval();
            }
            catch (mu::Parser::exception_type &e)
              TEST_FOR_EXCEPT_MSG(1, e.GetMsg());
          }
        }
      }
    }
  }

  // Update curcoord field to be consistent with initial displacement
  y->Update(1.0, *x, 1.0, *u, 0.0);
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyInitialVelocities(Teuchos::RCP<const Epetra_Vector> x,
                                                                            Teuchos::RCP<Epetra_Vector> v)
{
  const Epetra_BlockMap& threeDimensionalMap = v->Map();

  // apply the initial conditions
  Teuchos::ParameterList::ConstIterator it;
  for(it = params.begin() ; it != params.end() ; it++){
    const string & name = it->first;
    
    size_t position = name.find("Initial Velocity");
    if(position != string::npos){ // user wants to assign velocity using function
      
      Teuchos::ParameterList & boundaryConditionParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
      string nodeSet = boundaryConditionParams.get<string>("Node Set");
      string type = boundaryConditionParams.get<string>("Type");
      string coordinate = boundaryConditionParams.get<string>("Coordinate");
      string function = boundaryConditionParams.get<string>("Value");

      try{
        muParser.SetExpr(function);
      }
      catch (mu::Parser::exception_type &e)
        TEST_FOR_EXCEPT_MSG(1, e.GetMsg());

      int coord = 0;
      if(coordinate == "y" || coordinate == "Y")
        coord = 1;
      if(coordinate == "z" || coordinate == "Z")
        coord = 2;
      
      if (nodeSet == "All") { // Apply initial velocity to all locally-owned nodes
        for(int localNodeID = 0; localNodeID < x->MyLength(); localNodeID++) {
          muParserX = (*x)[localNodeID*3];
          muParserY = (*x)[localNodeID*3 + 1];
          muParserZ = (*x)[localNodeID*3 + 2];
          try{
            (*v)[localNodeID*3 + coord] = muParser.Eval();
          }
          catch (mu::Parser::exception_type &e)
            TEST_FOR_EXCEPT_MSG(1, e.GetMsg());
        }
      }
      else { // Apply initial velocity to specific node set
        // apply initial velocity boundary conditions to locally-owned nodes
        TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSet) == nodeSets->end(), "**** Node set not found: " + name + "\n");
        vector<int> & nodeList = (*nodeSets)[nodeSet];
        for(unsigned int i=0 ; i<nodeList.size() ; i++){
          int localNodeID = threeDimensionalMap.LID(nodeList[i]);
          if(localNodeID != -1) {
            muParserX = (*x)[localNodeID*3];
            muParserY = (*x)[localNodeID*3 + 1];
            muParserZ = (*x)[localNodeID*3 + 2];
            try {
              (*v)[localNodeID*3 + coord] = muParser.Eval();
            }
            catch (mu::Parser::exception_type &e)
              TEST_FOR_EXCEPT_MSG(1, e.GetMsg());
          }
        }
      }
    }
  }
}
