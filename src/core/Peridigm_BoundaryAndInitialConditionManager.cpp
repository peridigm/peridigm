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
#include <sstream>
#include <fstream>
#include <set>
#include <boost/algorithm/string/trim.hpp>
#include "Peridigm_Timer.hpp"
#include "Peridigm_Enums.hpp"
#include "Peridigm.hpp"

using namespace std;

PeridigmNS::BoundaryAndInitialConditionManager::BoundaryAndInitialConditionManager(const Teuchos::ParameterList& boundaryAndInitialConditionParams, Peridigm * peridigm_)
  : params(boundaryAndInitialConditionParams),
    peridigm(peridigm_){
}

void PeridigmNS::BoundaryAndInitialConditionManager::initialize(Teuchos::RCP<Discretization> discretization)
{
  bool hasPrescDisp = false;
  bool hasPrescVel = false;

  // Load node sets defined in the input deck into the nodeSets container
  for(Teuchos::ParameterList::ConstIterator it = params.begin() ; it != params.end() ; it++){
    string name = it->first;
    tidy_string(name);
    size_t position = name.find("NODE_SET");
    if(position!=string::npos) continue; // skip node set definitions, etc.
    if(!it->second.isList()){ // ensure that no other parameters are accidentally specified
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: Unknown parameter in boundary conditions specification " + name);
    }
    Teuchos::RCP<BoundaryCondition> bcPtr;
    Teuchos::ParameterList & bcParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
    const Boundary_Condition_Type bcType = to_boundary_condition_type(bcParams);

    // TODO: move this to a factory class
    switch(bcType)
    {
      case INITIAL_DISPLACEMENT:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getU();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm,false));
        initialConditions.push_back(bcPtr);
      }
      break;
      case INITIAL_VELOCITY :
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getV();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm,false));
        initialConditions.push_back(bcPtr);
      }
      break;
      case PRESCRIBED_DISPLACEMENT:
      {
        hasPrescDisp = true;
        // a prescribed displacement boundary condition will automatically update deltaU and the velocity vector
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getDeltaU();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVector,peridigm,false,1.0,0.0));
        boundaryConditions.push_back(bcPtr);
        Teuchos::RCP<Epetra_Vector> toVectorV = peridigm->getV();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVectorV,peridigm,false,0.0,1.0));
        boundaryConditions.push_back(bcPtr);
      }
      break;
      case INITIAL_TEMPERATURE:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getDeltaTemperature();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm,false));
        initialConditions.push_back(bcPtr);
      }
      break;
      case PRESCRIBED_TEMPERATURE:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getDeltaTemperature();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVector,peridigm,false,1.0,0.0));
        boundaryConditions.push_back(bcPtr);
      }
      break;
      case BODY_FORCE:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getExternalForce();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm,true));
        forceContributions.push_back(bcPtr);
      }
      break;
      case NO_SUCH_BOUNDARY_CONDITION_TYPE:
        break;
      default :
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: unknown boundary condition type " + to_string(bcType));
        break;
    }
  }
  TEUCHOS_TEST_FOR_EXCEPTION(hasPrescDisp&&hasPrescVel,std::logic_error,"Error: cannot specify prescribed displacement and velocity boundary conditions");

//      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"ERROR: the following parameter was not used in the initialization of boundary or initial conditions: " + it->first);

  initializeNodeSets(discretization);
}

void PeridigmNS::BoundaryAndInitialConditionManager::initializeNodeSets(Teuchos::RCP<Discretization> discretization)
{
  nodeSets = Teuchos::rcp(new map< string, vector<int> >());

  // Load node sets defined in the input deck into the nodeSets container
  for(Teuchos::ParameterList::ConstIterator it = params.begin() ; it != params.end() ; it++){
	string name = it->first;
	tidy_string(name);
	size_t position = name.find("NODE_SET");
	if(position != string::npos){

      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
      vector<int>& nodeList = (*nodeSets)[name];

      // Determine if the string is the name of a file or a list of node numbers
      string fileName = nodeSetStringToFileName(Teuchos::getValue<string>(it->second));

      if(fileName.size() == 0){
        // The string is a list of node numbers
        stringstream ss(Teuchos::getValue<string>(it->second));
        int nodeID;
        while(ss.good()){
          ss >> nodeID;
          // Convert from 1-based node numbering (Exodus II) to 0-based node numbering (Epetra and all the rest of Peridigm)
          TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeID < 1, "**** Error:  Node number 0 detected in nodeset definition; node numbering must begin with 1.\n");
          nodeList.push_back(nodeID - 1);
        }
      }
      else{
        // The string is the name of a file containing the node numbers
        ifstream inFile(fileName.c_str());
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!inFile.is_open(), "**** Error opening node set text file: " + fileName + "\n");
        while(inFile.good()){
          string str;
          getline(inFile, str);
          boost::trim(str);
          // Ignore comment lines, otherwise parse
          if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
            istringstream iss(str);
            vector<int> nodeNumbers;
            copy(istream_iterator<int>(iss),
                 istream_iterator<int>(),
                 back_inserter<vector<int> >(nodeNumbers));
            for(unsigned int i=0 ; i<nodeNumbers.size() ; ++i){
              // Convert from 1-based node numbering (Exodus II) to 0-based node numbering (Epetra and all the rest of Peridigm)
              TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeNumbers[i] < 1, "**** Error:  Node number 0 detected in nodeset file; node numbering must begin with 1.\n");
              nodeList.push_back(nodeNumbers[i] - 1);
            }
          }
        }
        inFile.close();
      }
    }
  }

  // Load node sets defined in the mesh file into the nodeSets container
  Teuchos::RCP< map< string, vector<int> > > discretizationNodeSets = discretization->getNodeSets();
  for(map< string, vector<int> >::iterator it=discretizationNodeSets->begin() ; it!=discretizationNodeSets->end() ; it++){
    string name = it->first;
    tidy_string(name);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
    vector<int>& nodeList = it->second;
    (*nodeSets)[name] = nodeList;
  }

  // Cull any off-processor nodes from the node lists
  Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap = discretization->getGlobalOwnedMap(1);
  for(map< string, vector<int> >::iterator it = nodeSets->begin() ; it != nodeSets->end() ; it++){
    vector<int>& nodeSet = it->second;
    vector<int>::iterator nIt = nodeSet.begin();
    while(nIt != nodeSet.end()){
      if(oneDimensionalMap->LID(*nIt) == -1)
        nIt = nodeSet.erase(nIt);
      else
        ++nIt;
    }
  }
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyInitialConditions(){
  for(unsigned i=0;i<initialConditions.size();++i){
    initialConditions[i]->apply(nodeSets);
    if(initialConditions[i]->getType() == INITIAL_DISPLACEMENT)
      updateCurrentCoordinates();
  }
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyBoundaryConditions(const double & timeCurrent, const double & timePrevious){
  for(unsigned i=0;i<boundaryConditions.size();++i){
    boundaryConditions[i]->apply(nodeSets,timeCurrent,timePrevious);
  }
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyForceContributions(const double & timeCurrent, const double & timePrevious){
  clearForceContributions();
  for(unsigned i=0;i<forceContributions.size();++i){
    forceContributions[i]->apply(nodeSets,timeCurrent,timePrevious);
  }
}

void PeridigmNS::BoundaryAndInitialConditionManager::clearForceContributions(){
  Teuchos::RCP<Epetra_Vector> externalForce = peridigm->getExternalForce();
  externalForce->PutScalar(0.0);
}

void PeridigmNS::BoundaryAndInitialConditionManager::updateCurrentCoordinates(){
  peridigm->getY()->Update(1.0, *(peridigm->getX()), 1.0, *(peridigm->getU()), 0.0);
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyKinematicBC_ComputeReactions(Teuchos::RCP<const Epetra_Vector> force,
                                                                                       Teuchos::RCP<Epetra_Vector> reaction)
{
  PeridigmNS::Timer::self().startTimer("Apply Boundary Conditions");
  reaction->PutScalar(0.0);
  const Epetra_BlockMap& threeDimensionalMap = force->Map();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(threeDimensionalMap.ElementSize() != 3, "**** applyKinematicBC_ComputeReactions() must be called with map having element size = 3.\n");


  // apply the boundary conditions
  for(unsigned i=0;i<boundaryConditions.size();++i)
  {
    Teuchos::RCP<BoundaryCondition> boundaryCondition = boundaryConditions[i];
    if(boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT)
    {
      const int coord = boundaryCondition->getCoord();
      const Set_Definition setDef = to_set_definition(boundaryCondition->getNodeSetName());
      // apply the bc to every element in the entire domain
      if(setDef==FULL_DOMAIN)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: Dirichlet conditions on the displacement cannot be prescribed over the full domain.");
      }
      // apply the bc only to specific node sets
      else{
        std::map< std::string, std::vector<int> >::iterator itBegin;
        std::map< std::string, std::vector<int> >::iterator itEnd;
        if (setDef == ALL_SETS){
          itBegin = nodeSets->begin();
          itEnd = nodeSets->end();
        }
        else{
          TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(boundaryCondition->getNodeSetName()) == nodeSets->end(), "**** Node set not found: " + boundaryCondition->getNodeSetName() + "\n");
          itBegin = nodeSets->find(boundaryCondition->getNodeSetName());
          itEnd = itBegin; itEnd++;
        }
        for(std::map<std::string,std::vector<int> > ::iterator setIt=itBegin;setIt!=itEnd;++setIt){
          vector<int> & nodeList = setIt->second;
          for(unsigned int i=0 ; i<nodeList.size() ; i++){
            int localNodeID = force->Map().LID(nodeList[i]);
            if(!force.is_null() && localNodeID != -1)
              (*reaction)[3*localNodeID + coord] = (*force)[3*localNodeID + coord];
          }
        }
      }
    }
  }
  PeridigmNS::Timer::self().stopTimer("Apply Boundary Conditions");
}

////! Add the external force to the residual
//void PeridigmNS::BoundaryAndInitialConditionManager::addExternalForceToResidual(Teuchos::RCP<Epetra_Vector> residual){
//  Teuchos::RCP<Epetra_Vector> forceVec = peridigm->getExternalForce();
//  TEUCHOS_TEST_FOR_EXCEPTION(residual->MyLength()!=forceVec->MyLength(),std::logic_error,"ERROR: the residual and external force vectors should be the same length.");
//  for(int i=0 ; i<residual->MyLength(); ++i)
//     (*residual)[i] += (*forceVec)[i];
//}

void PeridigmNS::BoundaryAndInitialConditionManager::applyKinematicBC_InsertZeros(Teuchos::RCP<Epetra_Vector> vec)
{
  PeridigmNS::Timer::self().startTimer("Apply Boundary Conditions");

  const Epetra_BlockMap& oneDimensionalMap = vec->Map();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(oneDimensionalMap.ElementSize() != 1, "**** applyKinematicBC_InsertZeros() must be called with map having element size = 1.\n");

  for(unsigned i=0;i<boundaryConditions.size();++i)
  {
    Teuchos::RCP<BoundaryCondition> boundaryCondition = boundaryConditions[i];
    if(boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT)
    {
      const int coord = boundaryCondition->getCoord();
      const Set_Definition setDef = to_set_definition(boundaryCondition->getNodeSetName());
      // apply the bc to every element in the entire domain
      if(setDef==FULL_DOMAIN)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: Dirichlet conditions on the displacement cannot be prescribed over the entire domain.");
      }
      // apply the bc only to specific node sets
      else{
        std::map< std::string, std::vector<int> >::iterator itBegin;
        std::map< std::string, std::vector<int> >::iterator itEnd;
        if (setDef == ALL_SETS){
          itBegin = nodeSets->begin();
          itEnd = nodeSets->end();
        }
        else{
          TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(boundaryCondition->getNodeSetName()) == nodeSets->end(), "**** Node set not found: " + boundaryCondition->getNodeSetName() + "\n");
          itBegin = nodeSets->find(boundaryCondition->getNodeSetName());
          itEnd = itBegin; itEnd++;
        }
        for(std::map<std::string,std::vector<int> > ::iterator setIt=itBegin;setIt!=itEnd;++setIt){
          vector<int> & nodeList = setIt->second;
          for(unsigned int i=0 ; i<nodeList.size() ; i++){
            int localNodeID = oneDimensionalMap.LID(3*nodeList[i]);
            if(!vec.is_null() && localNodeID != -1)
              (*vec)[localNodeID + coord] = 0.0;
          }
        }
      }
    }
  }
  PeridigmNS::Timer::self().stopTimer("Apply Boundary Conditions");
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyKinematicBC_InsertZerosAndSetDiagonal(Teuchos::RCP<Epetra_FECrsMatrix> mat)
{
  PeridigmNS::Timer::self().startTimer("Apply Boundary Conditions");

  // determine the L2 norm of the diagonal
  // this will be used to scale the diagonal entry for kinematic B.C.s
  Epetra_Vector diagonal(mat->Map());
  mat->ExtractDiagonalCopy(diagonal);
  double diagonalNorm1;
  diagonal.Norm1(&diagonalNorm1);
  double diagonalEntry = -1.0*diagonalNorm1/diagonal.GlobalLength();

  for(unsigned i=0;i<boundaryConditions.size();++i)
  {
    Teuchos::RCP<BoundaryCondition> boundaryCondition = boundaryConditions[i];
    if(boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT)
    {
      const int coord = boundaryCondition->getCoord();
      const Set_Definition setDef = to_set_definition(boundaryCondition->getNodeSetName());
      // apply the bc to every element in the entire domain
      if(setDef==FULL_DOMAIN)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,"ERROR: Dirichlet conditions on the displacement cannot be prescribed over the full domain.");
      }
      // apply the bc only to specific node sets
      else{
        std::map< std::string, std::vector<int> >::iterator itBegin;
        std::map< std::string, std::vector<int> >::iterator itEnd;
        if (setDef == ALL_SETS){
          itBegin = nodeSets->begin();
          itEnd = nodeSets->end();
        }
        else{
          TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(boundaryCondition->getNodeSetName()) == nodeSets->end(), "**** Node set not found: " + boundaryCondition->getNodeSetName() + "\n");
          itBegin = nodeSets->find(boundaryCondition->getNodeSetName());
          itEnd = itBegin; itEnd++;
        }
        for(std::map<std::string,std::vector<int> > ::iterator setIt=itBegin;setIt!=itEnd;++setIt){
          vector<int> & nodeList = setIt->second;

          // create data structures for inserting values into jacobian
          // an upper bound on the number of entries to set to zero is given by mat->NumMyCols()
          vector<double> jacobianValues(mat->NumMyCols(), 0.0);
          vector<int> jacobianIndices(mat->NumMyCols());
          vector<int> jacobianColIndices(mat->NumMyCols());
          for(unsigned int i=0 ; i<jacobianIndices.size() ; ++i)
            jacobianIndices[i] = i;

          // zero out the columns associated with kinematic boundary conditions:
          // create the list of columns only once:
          int columnIndex(0);
          for(unsigned int i=0 ; i<nodeList.size() ; i++){
            const int globalID = 3*nodeList[i] + coord;
            const int localColID = mat->LCID(globalID);
            if(localColID != -1)
              jacobianColIndices[columnIndex++] = localColID;
          }
          // iterate the local rows and set the approprate column values to 0
          int numEntriesToSetToZero = columnIndex;
          for(int iRow=0 ; iRow<mat->NumMyRows() ; ++iRow)
            mat->ReplaceMyValues(iRow, numEntriesToSetToZero, &jacobianValues[0], &jacobianColIndices[0]);

          for(unsigned int i=0 ; i<nodeList.size() ; i++){

            // zero out the row and put diagonalEntry on the diagonal
            int globalID = 3*nodeList[i] + coord;
            int localRowID = mat->LRID(globalID);
            int localColID = mat->LCID(globalID);

            // zero out the row and put diagonalEntry on the diagonal
            if(localRowID != -1){
              if(localColID != -1)
                jacobianValues[localColID] = diagonalEntry;
              // From Epetra_CrsMatrix documentation:
              // If a value is not already present for the specified location in the matrix, the
              // input value will be ignored and a positive warning code will be returned.
              // \todo Do the bookkeeping to send in data only for locations that actually exist in the matrix structure.
              mat->ReplaceMyValues(localRowID, mat->NumMyCols(), &jacobianValues[0], &jacobianIndices[0]);
              if(localColID != -1)
                jacobianValues[localColID] = 0.0;
            }
          }
        }
      }
    }
  }
  PeridigmNS::Timer::self().stopTimer("Apply Boundary Conditions");
}

string PeridigmNS::BoundaryAndInitialConditionManager::nodeSetStringToFileName(string str)
{
  // This function differentiates between two possible types of node set strings:
  // 1)  A list of nonzero integers (i.e., a list of node numbers)
  // 2)  The name of a file

  // If the string is determined to be the name of a file, the file name is returned
  // If the string is determined to be a list of node numbers, an empty string is returned.

  string emptyString = "";
  string whitespace = " \t";
  
  boost::trim(str);
  
  // If there is any whitespace then str is a list, not a file name
  if(str.find_first_of(whitespace) != std::string::npos)
    return emptyString;

  // If the string is not an integer, assume it is a file name
  char* pEnd;
  long int integerValue = strtol(str.c_str(), &pEnd, 0);
  if(integerValue != 0 || str == "0" || str == "+0" || str == "-0")
    return emptyString;
  return str;
}
