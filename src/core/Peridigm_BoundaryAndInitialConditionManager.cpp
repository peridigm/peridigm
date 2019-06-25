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
#include "Peridigm_DegreesOfFreedomManager.hpp"
#include <sstream>
#include <fstream>
#include <set>
#include <iterator>
#include "Peridigm_Timer.hpp"
#include "Peridigm_Enums.hpp"
#include "Peridigm.hpp"

using namespace std;

PeridigmNS::BoundaryAndInitialConditionManager::BoundaryAndInitialConditionManager(const Teuchos::ParameterList& boundaryAndInitialConditionParams, Peridigm * peridigm_)
  : params(boundaryAndInitialConditionParams), peridigm(peridigm_), createRankDeficientNodesNodeSet(false)
{
  if(params.isParameter("Create Node Set For Rank Deficient Nodes")){
    createRankDeficientNodesNodeSet = params.get<bool>("Create Node Set For Rank Deficient Nodes");
    // Remove the parameter to avoid glitches downstream in the processing of node sets
    params.remove("Create Node Set For Rank Deficient Nodes");
  }
}

void PeridigmNS::BoundaryAndInitialConditionManager::initialize(Teuchos::RCP<Discretization> discretization)
{
  bool hasPrescDisp = false;
  bool hasPrescVel = false;

  if(createRankDeficientNodesNodeSet)
    createRankDeficientBC();

  initializeNodeSets(discretization);

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

    switch(bcType)
    {
      case INITIAL_DISPLACEMENT:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getU();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm));
        initialConditions.push_back(bcPtr);
      }
      break;
      case INITIAL_VELOCITY :
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getV();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm));
        initialConditions.push_back(bcPtr);
      }
      break;
      case PRESCRIBED_DISPLACEMENT:
      {
        hasPrescDisp = true;
        // a prescribed displacement boundary condition will automatically update deltaU and the velocity vector
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getDeltaU();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVector,peridigm,1.0,0.0));
        boundaryConditions.push_back(bcPtr);
        Teuchos::RCP<Epetra_Vector> toVectorV = peridigm->getV();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVectorV,peridigm,0.0,1.0));
        boundaryConditions.push_back(bcPtr);

      }
      break;
      case PRESCRIBED_FLUID_PRESSURE_U:
      {
				// a prescribed fluid pressure boundary condition will automatically update deltaFluidPressureU
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getFluidPressureDeltaU();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVector,peridigm,1.0,0.0));
        boundaryConditions.push_back(bcPtr);

				Teuchos::RCP<Epetra_Vector> toVectorV = peridigm->getFluidPressureV();
				bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVectorV, peridigm,0.0,1.0));
				boundaryConditions.push_back(bcPtr);
      }
      break;
      case INITIAL_FLUID_PRESSURE_U:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getFluidPressureU();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm));
        initialConditions.push_back(bcPtr);
      }
      break;
      case INITIAL_TEMPERATURE:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getTemperature();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm));
        initialConditions.push_back(bcPtr);
      }
      break;
      case PRESCRIBED_TEMPERATURE:
      {
        Teuchos::RCP<Epetra_Vector> toVectorTemperature = peridigm->getTemperature();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVectorTemperature,peridigm));
        boundaryConditions.push_back(bcPtr);

        // Create a BC to evaluate the change in temperature (e.g., for thermal strains)
        bool computeChangeRelativeToInitialValue = true;
        Teuchos::RCP<Epetra_Vector> toVectorDeltaTemperature = peridigm->getDeltaTemperature();
        bcPtr = Teuchos::rcp(new DirichletIncrementBC(name,bcParams,toVectorDeltaTemperature,peridigm,1.0,0.0,computeChangeRelativeToInitialValue));
        boundaryConditions.push_back(bcPtr);
      }
      break;
      case THERMAL_FLUX:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getTemperature();
        bcPtr = Teuchos::rcp(new NeumannBC(name,bcParams,toVector,peridigm,nodeSets));
        boundaryConditions.push_back(bcPtr);
      }
      break;
      case BODY_FORCE:
      {
        Teuchos::RCP<Epetra_Vector> toVector = peridigm->getExternalForce();
        bcPtr = Teuchos::rcp(new DirichletBC(name,bcParams,toVector,peridigm));
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
}

void PeridigmNS::BoundaryAndInitialConditionManager::createRankDeficientBC()
{
  // create a placeholder boundary condition to dump nodes into when they have become rank deficient (this helps with convergence for quasistatic bond breaking)
  Teuchos::ParameterList containerParams;
  containerParams.set("Type","Prescribed_Displacement");
  containerParams.set("Name","Rank Deficient Nodes");
  containerParams.set("Node Set","RANK_DEFICIENT_NODES");
  containerParams.set("Value","0.0");
  containerParams.set("Coordinate","x");
  Teuchos::RCP<Epetra_Vector> toVector = peridigm->getDeltaU();
  Teuchos::RCP<BoundaryCondition> bcPtr = Teuchos::rcp(new DirichletIncrementBC(containerParams.get<string>("Name"),containerParams,toVector,peridigm,1.0,0.0));
  boundaryConditions.push_back(bcPtr);
  containerParams.set("Coordinate","y");
  bcPtr = Teuchos::rcp(new DirichletIncrementBC(containerParams.get<string>("Name"),containerParams,toVector,peridigm,1.0,0.0));
  containerParams.set("Coordinate","z");
  boundaryConditions.push_back(bcPtr);
  bcPtr = Teuchos::rcp(new DirichletIncrementBC(containerParams.get<string>("Name"),containerParams,toVector,peridigm,1.0,0.0));
  boundaryConditions.push_back(bcPtr);
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
          str = trim(str);
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

  // Create the bogus node set that will hold any rank deficient nodes that show up
  if(createRankDeficientNodesNodeSet)
    (*nodeSets)["RANK_DEFICIENT_NODES"] = vector<int>();

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
    if(initialConditions[i]->getType() == INITIAL_FLUID_PRESSURE_U)
      updateFluidPressureY();
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

void PeridigmNS::BoundaryAndInitialConditionManager::updateFluidPressureY(){
  peridigm->getFluidPressureY()->Update(1.0, *(peridigm->getFluidPressureU()), 0.0);
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyKinematicBC_ComputeReactions(Teuchos::RCP<const Epetra_Vector> force,
                                                                                       Teuchos::RCP<Epetra_Vector> reaction)
{
  PeridigmNS::Timer::self().startTimer("Apply Boundary Conditions");

  TEUCHOS_TEST_FOR_EXCEPT_MSG(force->Map().ElementSize() != reaction->Map().ElementSize(), "**** In applyKinematicBC_ComputeReactions(), incompatible vector maps.");

  reaction->PutScalar(0.0);

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  int numDof = dofManager.totalNumberOfDegreesOfFreedom();

  for(unsigned i=0;i<boundaryConditions.size();++i) {

    Teuchos::RCP<BoundaryCondition> boundaryCondition = boundaryConditions[i];

    if ((boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT && dofManager.displacementTreatedAsUnknown()) ||
       (boundaryCondition->getType() == PRESCRIBED_FLUID_PRESSURE_U && dofManager.pressureTreatedAsUnknown()) ||
       (boundaryCondition->getType() == PRESCRIBED_TEMPERATURE && dofManager.temperatureTreatedAsUnknown()) ||
       (boundaryCondition->getType() == THERMAL_FLUX && dofManager.temperatureTreatedAsUnknown())) {

      int offset = 0;
      if (boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT) {
        offset = dofManager.displacementDofOffset();
      }
      else if (boundaryCondition->getType() == PRESCRIBED_FLUID_PRESSURE_U) {
        offset = dofManager.pressureDofOffset();
      }
      else if (boundaryCondition->getType() == PRESCRIBED_TEMPERATURE || boundaryCondition->getType() == THERMAL_FLUX) {
        offset = dofManager.temperatureDofOffset();
      }

      int coord = boundaryCondition->getCoord();

      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(boundaryCondition->getNodeSetName()) == nodeSets->end(),
                                  "**** Error in applyKinematicBC_ComputeReactions(), node set not found: " + boundaryCondition->getNodeSetName() + "\n");
      std::map< std::string, std::vector<int> >::iterator setIt = nodeSets->find(boundaryCondition->getNodeSetName());
      vector<int> & nodeList = setIt->second;
      for(unsigned int i=0 ; i<nodeList.size() ; i++){
        int localNodeID = force->Map().LID(nodeList[i]);
        if(!force.is_null() && localNodeID != -1)
          (*reaction)[numDof * localNodeID + offset + coord] = (*force)[numDof * localNodeID + offset + coord];
      }
    }
  }
  PeridigmNS::Timer::self().stopTimer("Apply Boundary Conditions");
}

void PeridigmNS::BoundaryAndInitialConditionManager::applyKinematicBC_InsertZeros(Teuchos::RCP<Epetra_Vector> vec)
{
  PeridigmNS::Timer::self().startTimer("Apply Boundary Conditions");

  const Epetra_BlockMap& oneDimensionalMap = vec->Map();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(oneDimensionalMap.ElementSize() != 1, "**** applyKinematicBC_InsertZeros() must be called with map having element size = 1.\n");

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  int numDof = dofManager.totalNumberOfDegreesOfFreedom();

  for (unsigned i=0;i<boundaryConditions.size();++i) {

    Teuchos::RCP<BoundaryCondition> boundaryCondition = boundaryConditions[i];

    if ((boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT && dofManager.displacementTreatedAsUnknown()) ||
       (boundaryCondition->getType() == PRESCRIBED_FLUID_PRESSURE_U && dofManager.pressureTreatedAsUnknown()) ||
       (boundaryCondition->getType() == PRESCRIBED_TEMPERATURE && dofManager.temperatureTreatedAsUnknown()) ||
       (boundaryCondition->getType() == THERMAL_FLUX && dofManager.temperatureTreatedAsUnknown())) {

      int offset = 0;
      if (boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT) {
        offset = dofManager.displacementDofOffset();
      }
      else if (boundaryCondition->getType() == PRESCRIBED_FLUID_PRESSURE_U) {
        offset = dofManager.pressureDofOffset();
      }
      else if (boundaryCondition->getType() == PRESCRIBED_TEMPERATURE || boundaryCondition->getType() == THERMAL_FLUX) {
        offset = dofManager.temperatureDofOffset();
      }

      int coord = boundaryCondition->getCoord();

      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(boundaryCondition->getNodeSetName()) == nodeSets->end(),
                                  "**** Error in applyKinematicBC_ComputeReactions(), node set not found: " + boundaryCondition->getNodeSetName() + "\n");
      std::map< std::string, std::vector<int> >::iterator setIt = nodeSets->find(boundaryCondition->getNodeSetName());
      vector<int> & nodeList = setIt->second;
      for(unsigned int i=0 ; i<nodeList.size() ; i++){
        int localNodeID = oneDimensionalMap.LID(numDof * nodeList[i]);
        if(!vec.is_null() && localNodeID != -1)
          (*vec)[localNodeID + offset + coord] = 0.0;
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

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  int numDof = dofManager.totalNumberOfDegreesOfFreedom();

  for (unsigned i=0;i<boundaryConditions.size();++i) {

    Teuchos::RCP<BoundaryCondition> boundaryCondition = boundaryConditions[i];

    if ((boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT && dofManager.displacementTreatedAsUnknown()) ||
       (boundaryCondition->getType() == PRESCRIBED_FLUID_PRESSURE_U && dofManager.pressureTreatedAsUnknown()) ||
       (boundaryCondition->getType() == PRESCRIBED_TEMPERATURE && dofManager.temperatureTreatedAsUnknown()) ||
       (boundaryCondition->getType() == THERMAL_FLUX && dofManager.temperatureTreatedAsUnknown())) {

      int offset = 0;
      if (boundaryCondition->getType() == PRESCRIBED_DISPLACEMENT) {
        offset = dofManager.displacementDofOffset();
      }
      else if (boundaryCondition->getType() == PRESCRIBED_FLUID_PRESSURE_U) {
        offset = dofManager.pressureDofOffset();
      }
      else if (boundaryCondition->getType() == PRESCRIBED_TEMPERATURE || boundaryCondition->getType() == THERMAL_FLUX) {
        offset = dofManager.temperatureDofOffset();
      }

      int coord = boundaryCondition->getCoord();

      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(boundaryCondition->getNodeSetName()) == nodeSets->end(),
                                  "**** Error in applyKinematicBC_InsertZerosAndSetDiagonal(), node set not found: " + boundaryCondition->getNodeSetName() + "\n");
      std::map< std::string, std::vector<int> >::iterator setIt = nodeSets->find(boundaryCondition->getNodeSetName());
      vector<int> & nodeList = setIt->second;

      // create data structures for inserting values into jacobian
      // an upper bound on the number of entries to set to zero is given by mat->NumMyCols()
      vector<double> jacobianValues(mat->NumMyCols(), 0.0);
      vector<int> jacobianIndices(mat->NumMyCols());
      vector<int> jacobianColIndices(mat->NumMyCols());
      for(unsigned int i=0 ; i<jacobianIndices.size() ; ++i){
        jacobianIndices[i] = i;
      }

      // zero out the columns associated with kinematic boundary conditions:
      // create the list of columns only once:
      int columnIndex(0);
      for(unsigned int i=0 ; i<nodeList.size() ; i++){
        const int globalID = numDof * nodeList[i] + offset + coord;
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
        int globalID = numDof * nodeList[i] + offset + coord;
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

  str = trim(str);

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
