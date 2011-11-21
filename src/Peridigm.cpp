/*! \file Peridigm.cpp
 *
 * File containing main class for Peridigm: A parallel, multi-physics,
 * peridynamics simulation code.
 */

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

#include <iostream>
#include <vector>
#include <map>

#include <boost/math/special_functions/fpclassify.hpp>

#include <AztecOO.h>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>
#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <EpetraExt_BlockMapOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Teuchos_VerboseObject.hpp>

#include "Peridigm.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "Peridigm_PdQuickGridDiscretization.hpp"
#include "Peridigm_OutputManager_VTK_XML.hpp"
#include "Peridigm_ComputeManager.hpp"
#include "materials/Peridigm_MaterialFactory.hpp"
#include "contact/Peridigm_ContactModelFactory.hpp"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "mesh_input/quick_grid/QuickGridData.h"
#include "pdneigh/PdZoltan.h"
#include "pdneigh/NeighborhoodList.h"
#include "InitialCondition.hpp"
#include "Peridigm_Timer.hpp"

using namespace std;

PeridigmNS::Peridigm::Peridigm(const Teuchos::RCP<const Epetra_Comm>& comm,
                   const Teuchos::RCP<Teuchos::ParameterList>& params)
  : numBlocks(1),
    analysisHasRebalance(false),
    rebalanceFrequency(1),
    analysisHasContact(false),
    contactRebalanceFrequency(0),
    contactSearchRadius(0.0)
{
  peridigmComm = comm;
  peridigmParams = params;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // Instantiate materials and associate them with the blocks
  instantiateMaterials();

  // Read mesh from disk or generate using geometric primatives.
  Teuchos::RCP<Teuchos::ParameterList> discParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Problem", true).sublist("Discretization", true) );
  DiscretizationFactory discFactory(discParams);
  Teuchos::RCP<AbstractDiscretization> peridigmDisc = discFactory.create(peridigmComm);
  initializeDiscretization(peridigmDisc);

  // Load node sets from input deck and/or input mesh file into nodeSets container
  Teuchos::RCP<Teuchos::ParameterList> bcParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Problem", true).sublist("Boundary Conditions") );
  initializeNodeSets(bcParams, peridigmDisc);

  // The PeridigmNS::DataManager contained in each PeridigmNS::Block allocates space for field data.
  // Keep track of the data storage required by Peridigm and the ComputeManager and pass the associated
  // field specs when calling Block::initializeDataManager().
  Teuchos::RCP< std::vector<Field_NS::FieldSpec> > fieldSpecs = Teuchos::rcp(new std::vector<Field_NS::FieldSpec>);

  // Load the Peridigm specs into fieldSpecs
  std::vector<Field_NS::FieldSpec> peridigmSpecs = this->getFieldSpecs();
  fieldSpecs->insert(fieldSpecs->end(), peridigmSpecs.begin(), peridigmSpecs.end());

  // Initialize compute manager
  initializeComputeManager();

  // Load the compute manager's field specs into fieldSpecs
  std::vector<Field_NS::FieldSpec> computeSpecs = computeManager->getFieldSpecs();
  fieldSpecs->insert(fieldSpecs->end(), computeSpecs.begin(), computeSpecs.end());

  // Instantiate the blocks
  blocks = Teuchos::rcp(new std::vector<PeridigmNS::Block>());
  std::vector<std::string> blockNames = peridigmDisc->getBlockNames();
  for(unsigned int iBlock=0 ; iBlock<blockNames.size() ; ++iBlock){
    std::string blockName = blockNames[iBlock];
    PeridigmNS::Block block(blockName, iBlock+1);
    blocks->push_back(block);
  }

  // Load the material models into the blocks
  // \todo Hook up plumbing for multiple materials in the input deck, for now all blocks get the same material.
  Teuchos::RCP<const Teuchos::ParameterList> materialParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Problem", true).sublist("Material", true) );
  MaterialFactory materialFactory;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->setMaterialModel( materialFactory.create(materialParams) );

  // Load the auxiliary field specs into the blocks (they will be
  // combined with material model and contact model specs when allocating
  // space in the Block's DataManager)
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->setAuxiliaryFieldSpecs(fieldSpecs);    

  // Initialize the blocks (creates maps, neighborhoods, DataManager)
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->initialize(peridigmDisc->getGlobalOwnedMap(1),
                        peridigmDisc->getGlobalOverlapMap(1),
                        peridigmDisc->getGlobalOwnedMap(3),
                        peridigmDisc->getGlobalOverlapMap(3),
                        peridigmDisc->getGlobalBondMap(),
                        blockIDs,
                        globalNeighborhoodData);

  // Load initial data into the blocks
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*(peridigmDisc->getCellVolume()), Field_NS::VOLUME,     Field_ENUM::STEP_NONE, Insert);
    blockIt->importData(*(peridigmDisc->getInitialX()),   Field_NS::COORD3D,    Field_ENUM::STEP_NONE, Insert);
    blockIt->importData(*(peridigmDisc->getInitialX()),   Field_NS::CURCOORD3D, Field_ENUM::STEP_N,    Insert);
    blockIt->importData(*(peridigmDisc->getInitialX()),   Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1,  Insert);
  }

  // apply initial velocities
  applyInitialVelocities();

  // apply initial displacements
  applyInitialDisplacements();

  // Initialize material models
  // Initialization functions require valid initial values, e.g. velocities and displacements.
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->initializeMaterialModel();

  // Setup contact
  initializeContact();

  // Initialize the workset
  initializeWorkset();

  // Create the model evaluator
  modelEvaluator = Teuchos::rcp(new PeridigmNS::ModelEvaluator(analysisHasContact));

  // Initialize output manager
  initializeOutputManager();

  // Call rebalance function if analysis has contact
  // this is required to set up proper contact neighbor list
  if(analysisHasContact)
    rebalance();
}

void PeridigmNS::Peridigm::instantiateMaterials() {

  // Extract problem parameters sublist
  Teuchos::ParameterList& problemParams = peridigmParams->sublist("Problem", true);
  Teuchos::ParameterList& materialParams = problemParams.sublist("Material", true);

  // Instantiate material objects
  // \todo This will break for multiple material blocks.
  MaterialFactory materialFactory;
  materialModels = Teuchos::rcp(new std::vector< Teuchos::RCP<const PeridigmNS::Material> >()); 
  for(Teuchos::ParameterList::ConstIterator it = materialParams.begin() ; it != materialParams.end() ; it++){
    Teuchos::RCP<const Teuchos::ParameterList> matParams = Teuchos::rcpFromRef(materialParams);
    materialModels->push_back( materialFactory.create(matParams) );
  }
  TEST_FOR_EXCEPT_MSG(materialModels->size() == 0, "No material models created!");
}

void PeridigmNS::Peridigm::initializeDiscretization(Teuchos::RCP<AbstractDiscretization> peridigmDisc) {

  // oneDimensionalMap
  // used for cell volumes and scalar constitutive data
  oneDimensionalMap = peridigmDisc->getGlobalOwnedMap(1); 

  // oneDimensionalOverlapMap
  // used for cell volumes and scalar constitutive data
  // includes ghosts
  oneDimensionalOverlapMap = peridigmDisc->getGlobalOverlapMap(1);

  // threeDimensionalMap
  // used for positions, displacements, velocities and vector constitutive data
  threeDimensionalMap = peridigmDisc->getGlobalOwnedMap(3);

  // bondConstitutiveDataMap
  // a non-overlapping map used for storing constitutive data on bonds
  bondMap = peridigmDisc->getGlobalBondMap();

  // Create mothership vectors

  oneDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*oneDimensionalMap, 2));
  blockIDs = Teuchos::rcp((*oneDimensionalMothership)(0), false);        // block ID
  volume = Teuchos::rcp((*oneDimensionalMothership)(1), false);          // cell volume

  // \todo Do not allocate space for the contact force, residual, and deltaU if not needed.
  threeDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*threeDimensionalMap, 10));
  // Set ref-count pointers for each of the global vectors
  x = Teuchos::rcp((*threeDimensionalMothership)(0), false);             // initial positions
  u = Teuchos::rcp((*threeDimensionalMothership)(1), false);             // displacement
  y = Teuchos::rcp((*threeDimensionalMothership)(2), false);             // current positions
  v = Teuchos::rcp((*threeDimensionalMothership)(3), false);             // velocities
  a = Teuchos::rcp((*threeDimensionalMothership)(4), false);             // accelerations
  force = Teuchos::rcp((*threeDimensionalMothership)(5), false);         // force
  contactForce = Teuchos::rcp((*threeDimensionalMothership)(6), false);  // contact force (used only for contact simulations)
  deltaU = Teuchos::rcp((*threeDimensionalMothership)(7), false);        // increment in displacement (used only for implicit time integration)
  residual = Teuchos::rcp((*threeDimensionalMothership)(8), false);      // residual (used only for implicit time integration)
  scratch = Teuchos::rcp((*threeDimensionalMothership)(9), false);       // scratch space

  // Set the block IDs
  double* bID;
  peridigmDisc->getBlockID()->ExtractView(&bID);
  double* blockIDsPtr;
  blockIDs->ExtractView(&blockIDsPtr);
  blas.COPY(blockIDs->MyLength(), bID, blockIDsPtr);

  // Set the volumes
  double* vol;
  peridigmDisc->getCellVolume()->ExtractView(&vol);
  double* volumePtr;
  volume->ExtractView(&volumePtr);
  blas.COPY(volume->MyLength(), vol, volumePtr);

  // Set the initial positions
  double* initialX;
  peridigmDisc->getInitialX()->ExtractView(&initialX);
  double* xPtr;
  x->ExtractView(&xPtr);
  blas.COPY(x->MyLength(), initialX, xPtr);
  double* yPtr;
  y->ExtractView(&yPtr);
  blas.COPY(y->MyLength(), initialX, yPtr);

  // get the neighborlist from the discretization
  globalNeighborhoodData = peridigmDisc->getNeighborhoodData();
}

void PeridigmNS::Peridigm::initializeNodeSets(Teuchos::RCP<Teuchos::ParameterList>& bcParams,
                                              Teuchos::RCP<AbstractDiscretization> peridigmDisc) {

  nodeSets = Teuchos::rcp(new map< string, vector<int> >());

  // Load node sets defined in the input deck into the nodeSets container
  for(Teuchos::ParameterList::ConstIterator it = bcParams->begin() ; it != bcParams->end() ; it++){
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
  Teuchos::RCP< map< string, vector<int> > > discretizationNodeSets = peridigmDisc->getNodeSets();
  for(map< string, vector<int> >::iterator it=discretizationNodeSets->begin() ; it!=discretizationNodeSets->end() ; it++){
    string name = it->first;
    TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
    vector<int>& nodeList = it->second;
    (*nodeSets)[name] = nodeList;
  }
}

void PeridigmNS::Peridigm::applyInitialVelocities() {

  TEST_FOR_EXCEPT_MSG(!threeDimensionalMap->SameAs(v->Map()), 
                      "Peridigm::applyInitialVelocities():  Inconsistent velocity vector map.\n");

  /*
   * UNCOMMENT the following two lines
   */

//  RCP<InitialConditionsNS::InitialCondition> icOperator = InitialConditionsNS::getInstance(*peridigmParams);
//  icOperator->apply(*x,*u,*v);

  /*
   * COMMENT OUT ALL OF BELOW TO RUN new Initial Condition Capability
   */
  // apply the initial conditions
  Teuchos::ParameterList& bcParams = peridigmParams->sublist("Problem").sublist("Boundary Conditions");
  Teuchos::ParameterList::ConstIterator it;
  for(it = bcParams.begin() ; it != bcParams.end() ; it++){
	const string & name = it->first;
	size_t position = name.find("Initial Velocity");
	if(position != string::npos){
	  Teuchos::ParameterList & boundaryConditionParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
	  string nodeSet = boundaryConditionParams.get<string>("Node Set");
	  string type = boundaryConditionParams.get<string>("Type");
	  string coordinate = boundaryConditionParams.get<string>("Coordinate");
	  double value = boundaryConditionParams.get<double>("Value");

	  int coord = 0;
	  if(coordinate == "y" || coordinate == "Y")
		coord = 1;
	  if(coordinate == "z" || coordinate == "Z")
		coord = 2;

	  // apply initial velocity boundary conditions
	  // to locally-owned nodes
      TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSet) == nodeSets->end(), "**** Node set not found: " + name + "\n");
	  vector<int> & nodeList = (*nodeSets)[nodeSet];
	  for(unsigned int i=0 ; i<nodeList.size() ; i++){
		int localNodeID = threeDimensionalMap->LID(nodeList[i]);
		if(localNodeID != -1)
		  (*v)[localNodeID*3 + coord] = value;
	  }
	}
  }

  // Fill the dataManager with initial velocity data
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->importData(*v, Field_NS::VELOC3D, Field_ENUM::STEP_N, Insert);
}

void PeridigmNS::Peridigm::applyInitialDisplacements() {

  TEST_FOR_EXCEPT_MSG(!threeDimensionalMap->SameAs(u->Map()),
                      "Peridigm::applyInitialDisplacements():  Inconsistent displacement vector map.\n");

  Teuchos::ParameterList& problemParams = peridigmParams->sublist("Problem");
  Teuchos::ParameterList& bcParams = problemParams.sublist("Boundary Conditions");
  Teuchos::ParameterList::ConstIterator it;

  // apply the initial conditions
  for(it = bcParams.begin() ; it != bcParams.end() ; it++){
        const string & name = it->first;
        size_t position = name.find("Initial Displacement");
        if(position != string::npos){
          Teuchos::ParameterList & boundaryConditionParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
          string nodeSet = boundaryConditionParams.get<string>("Node Set");
          string type = boundaryConditionParams.get<string>("Type");
          string coordinate = boundaryConditionParams.get<string>("Coordinate");
          double value = boundaryConditionParams.get<double>("Value");

          int coord = 0;
          if(coordinate == "y" || coordinate == "Y")
                coord = 1;
          if(coordinate == "z" || coordinate == "Z")
                coord = 2;

          // apply initial displacement boundary conditions
          // to locally-owned nodes
          TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSet) == nodeSets->end(), "**** Node set not found: " + name + "\n");
          vector<int> & nodeList = (*nodeSets)[nodeSet];
          for(unsigned int i=0 ; i<nodeList.size() ; i++){
                int localNodeID = threeDimensionalMap->LID(nodeList[i]);
                if(localNodeID != -1)
                  (*u)[localNodeID*3 + coord] = value;
          }
        }
  }

  // Update curcoord field to be consistent with initial displacement
  y->Update(1.0, *x, 1.0, *u, 0.0);

  // Fill the dataManagers with initial displacement, position data
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, Field_NS::DISPL3D, Field_ENUM::STEP_N, Insert);
    blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_N, Insert);
  }
}

void PeridigmNS::Peridigm::initializeContact() {

  // Extract problem parameters sublist
  Teuchos::RCP<Teuchos::ParameterList> problemParams = Teuchos::rcp(&(peridigmParams->sublist("Problem")),false);

  // Extract discretization parameters sublist
  Teuchos::RCP<Teuchos::ParameterList> discParams = Teuchos::rcp(&(problemParams->sublist("Discretization")), false);

  // Assume no contact
  analysisHasContact = false;
  contactSearchRadius = 0.0;
  contactRebalanceFrequency = 0;

  // Set up global contact parameters for rebalance and proximity search
  if(problemParams->isSublist("Contact")){
    Teuchos::ParameterList & contactParams = problemParams->sublist("Contact");
    analysisHasContact = true;
    if(!contactParams.isParameter("Search Radius"))
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Radius\" not specified.");
    contactSearchRadius = contactParams.get<double>("Search Radius");
    if(!contactParams.isParameter("Search Frequency"))
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Frequency\" not specified.");
    contactRebalanceFrequency = contactParams.get<int>("Search Frequency");
  }

  // Instantiate contact models
  ContactModelFactory contactModelFactory;
  contactModels = Teuchos::rcp(new std::vector<Teuchos::RCP<const PeridigmNS::ContactModel> >);
  if(analysisHasContact){
    Teuchos::ParameterList& contactParams = problemParams->sublist("Contact");
    Teuchos::ParameterList::ConstIterator it;
    for(it = contactParams.begin() ; it != contactParams.end() ; it++){
      const string & name = it->first;
      // assume that the sublist is a set of parameters for a contact model
      if(contactParams.isSublist(name)){
        Teuchos::ParameterList contactModelParams;
        contactModelParams.set(name, contactParams.sublist(name));
        // Add the horizon to the contact model parameters, if needed
        if(!contactModelParams.sublist(name).isParameter("Horizon"))
          contactModelParams.sublist(name).set("Horizon", discParams->get<double>("Horizon"));

        // \todo This will break for multiple material blocks, need scheme to associate blocks with contact models
        for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
          blockIt->setContactModel( contactModelFactory.create(contactModelParams) );
      }
    }
  }
}

void PeridigmNS::Peridigm::initializeWorkset() {
  workset = Teuchos::rcp(new PHAL::Workset);
  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  *timeStep = 0.0;
  workset->timeStep = timeStep;
  workset->jacobian = overlapJacobian;
  workset->blocks = blocks;
  workset->myPID = -1;
}

void PeridigmNS::Peridigm::initializeComputeManager() {

  Teuchos::RCP<Teuchos::ParameterList> outputParams;

  if (peridigmParams->isSublist("Output")) {
    outputParams  = Teuchos::rcp(&(peridigmParams->sublist("Output")),false);
  }

  computeManager = Teuchos::rcp( new PeridigmNS::ComputeManager( outputParams, this  ) );

}


void PeridigmNS::Peridigm::initializeOutputManager() {

  bool active = false;
  Teuchos::RCP<Teuchos::ParameterList> outputParams;

  if (peridigmParams->isSublist("Output")) {
    active = true;
    // Get reference to existing sublit
    Teuchos::ParameterList& masterList = peridigmParams->sublist("Output");
    // Make copy of master list
    outputParams  = Teuchos::rcp( new Teuchos::ParameterList(masterList) );
    // Add proc id data to copied list
    outputParams->set("NumProc", (int)(peridigmComm->NumProc()));
    outputParams->set("MyPID", (int)(peridigmComm->MyPID()));
  }

  if (active) {
    // Make the default format "VTK_XML"
    string outputFormat = outputParams->get("Output File Type", "VTK_XML");
    TEST_FOR_EXCEPTION( outputFormat != "VTK_XML",
                        std::invalid_argument,
                        "PeridigmNS::Peridigm: \"Output File Type\" must be \"VTK_XML\".");
    if (outputFormat == "VTK_XML")
       outputManager = Teuchos::rcp(new PeridigmNS::OutputManager_VTK_XML( outputParams, this, blocks ));
    else
      TEST_FOR_EXCEPTION( true, std::invalid_argument,"PeridigmNS::Peridigm::initializeOutputManager: \"Output File Type\" must be \"VTK_XML\".");

    // Initialize current time in this parameterlist
    Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::rcp(&(peridigmParams->sublist("Solver")),false);
    //double timeInitial = solverParams->get("Initial Time", 0.0);
    // Initial conditions to disk written by time integrators before taking first step
    //this->synchDataManagers();
    //outputManager->write(dataManager,neighborhoodData,timeInitial);
  }
  else { // no output requested
    outputManager = Teuchos::rcp(new PeridigmNS::OutputManager_VTK_XML( outputParams, this, blocks ));
  }

  //  verbose = problemParams->get("Verbose", false);

}

void PeridigmNS::Peridigm::execute() {

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);

  // allowable explicit time integration schemes:  Verlet
  if(solverParams->isSublist("Verlet"))
    executeExplicit();

  // allowable implicit time integration schemes:  Implicit, QuasiStatic
  else if(solverParams->isSublist("QuasiStatic"))    
    executeQuasiStatic();
  else if(solverParams->isSublist("Implicit"))    
    executeImplicit();
}

void PeridigmNS::Peridigm::executeExplicit() {

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, Field_NS::DISPL3D,    Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*v, Field_NS::VELOC3D,    Field_ENUM::STEP_NP1, Insert);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  Teuchos::RCP<Teuchos::ParameterList> verletParams = sublist(solverParams, "Verlet", true);
  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal   = solverParams->get("Final Time", 1.0);
  double timeCurrent = timeInitial;
  double dt = verletParams->get("Fixed dt", 1.0);
  *timeStep = dt;
  double dt2 = dt/2.0;
  int nsteps = (int)floor((timeFinal-timeInitial)/dt);
  if(solverParams->isSublist("Rebalance")){
    Teuchos::RCP<Teuchos::ParameterList> rebalanceParams = sublist(solverParams, "Rebalance", true);
    analysisHasRebalance = true;
    rebalanceFrequency = rebalanceParams->get("Rebalance Frequency", 1);
  }
  // Pointer index into sub-vectors for use with BLAS
  double *xptr, *uptr, *yptr, *vptr, *aptr;
  x->ExtractView( &xptr );
  u->ExtractView( &uptr );
  y->ExtractView( &yptr );
  v->ExtractView( &vptr );
  a->ExtractView( &aptr );
  int length = a->MyLength();

  // Evaluate force in initial configuration for use in first timestep

  // \todo The velocity copied into the DataManager is actually the midstep velocity, not the NP1 velocity; this can be fixed by creating a midstep velocity field in the DataManager and setting the NP1 value as invalid.

  // Update forces based on initial displacements
  PeridigmNS::Timer::self().startTimer("Model Evaluator");
  modelEvaluator->evalModel(workset);
  PeridigmNS::Timer::self().stopTimer("Model Evaluator");

  // Copy force from the data manager to the mothership vector
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  force->PutScalar(0.0);
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    scratch->PutScalar(0.0);
    blockIt->exportData(*scratch, Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
    force->Update(1.0, *scratch, 1.0);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  if(analysisHasContact){
    // Copy contact force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    contactForce->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(*scratch, Field_NS::CONTACT_FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
      contactForce->Update(1.0, *scratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
    // Add contact forces to forces
    force->Update(1.0, *contactForce, 1.0);
  }

  // fill the acceleration vector
  (*a) = (*force);
  // \todo Possibly move this functionality into ModelEvaluator.
  // \todo This will break for multiple materials
  double density = (*materialModels)[0]->Density();
  a->Scale(1.0/density);

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  this->synchDataManagers();

  outputManager->write(blocks,timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  int displayTrigger = nsteps/100;
  if(displayTrigger == 0)
    displayTrigger = 1;

  for(int step=1; step<=nsteps; step++){

    if((step-1)%displayTrigger==0)
      displayProgress("Explicit time integration", (step-1)*100.0/nsteps);

    // rebalance, if requested
    if( (analysisHasRebalance && step%rebalanceFrequency == 0) || (analysisHasContact && step%contactRebalanceFrequency == 0) ){
      PeridigmNS::Timer::self().startTimer("Rebalance");
      rebalance();
      PeridigmNS::Timer::self().stopTimer("Rebalance");
      x->ExtractView( &xptr );
      u->ExtractView( &uptr );
      y->ExtractView( &yptr );
      v->ExtractView( &vptr );
      a->ExtractView( &aptr );
      length = a->MyLength();
    }

    // Do one step of velocity-Verlet

    // V^{n+1/2} = V^{n} + (dt/2)*A^{n}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt2, aptr, vptr, 1, 1);

    // Y^{n+1} = X_{o} + U^{n} + (dt)*V^{n+1/2}
    // \todo Replace with blas call
    for(int i=0 ; i<y->MyLength() ; ++i)
      yptr[i] = xptr[i] + uptr[i] + dt*vptr[i];

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt, vptr, uptr, 1, 1);

    // \todo The velocity copied into the DataManager is actually the midstep velocity, not the NP1 velocity; this can be fixed by creating a midstep velocity field in the DataManager and setting the NP1 value as invalid.

    // Copy data from mothership vectors to overlap vectors in data manager
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(*u, Field_NS::DISPL3D,    Field_ENUM::STEP_NP1, Insert);
      blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1, Insert);
      blockIt->importData(*v, Field_NS::VELOC3D,    Field_ENUM::STEP_NP1, Insert);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Model Evaluator");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Model Evaluator");

    // Copy force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    force->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(*scratch, Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
      force->Update(1.0, *scratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");    

    // Check for NaNs in force evaluation
    // We'd like to know now because a NaN will likely cause a difficult-to-unravel crash downstream.
    for(int i=0 ; i<force->MyLength() ; ++i)
      TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*force)[i]), "**** NaN returned by force evaluation.\n");

    if(analysisHasContact){
      // Copy contact force from the data manager to the mothership vector
      PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      contactForce->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scratch->PutScalar(0.0);
        blockIt->exportData(*scratch, Field_NS::CONTACT_FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
        contactForce->Update(1.0, *scratch, 1.0);
      }
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Check for NaNs in contact force evaluation
      for(int i=0 ; i<contactForce->MyLength() ; ++i)
        TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*contactForce)[i]), "**** NaN returned by contact force evaluation.\n");

      // Add contact forces to forces
      force->Update(1.0, *contactForce, 1.0);
    }

    // fill the acceleration vector
    (*a) = (*force);
    // \todo Possibly move this functionality into ModelEvaluator.
    // \todo This will break for multiple materials.
    double density = (*materialModels)[0]->Density();
    a->Scale(1.0/density);

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt2, aptr, vptr, 1, 1);

    timeCurrent = timeInitial + (step*dt);

    PeridigmNS::Timer::self().startTimer("Output");
    this->synchDataManagers();

    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->updateState();
  }
  displayProgress("Explicit time integration", 100.0);
  *out << "\n\n";
}

void PeridigmNS::Peridigm::executeQuasiStatic() {

  // Allocate memory for non-zeros in global Jacobain and lock in the structure
  allocateJacobian();

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  Teuchos::RCP<Teuchos::ParameterList> implicitParams = sublist(solverParams, "QuasiStatic", true);
  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal = solverParams->get("Final Time", 1.0);
  double timeCurrent = timeInitial;
  int numLoadSteps = implicitParams->get("Number of Load Steps", 10);
  double absoluteTolerance = implicitParams->get("Absolute Tolerance", 1.0e-6);
  double maximumSolverIterations = implicitParams->get("Maximum Solver Iterations", 10);

  // Pointer index into sub-vectors for use with BLAS
  double *xptr, *uptr, *yptr, *vptr, *aptr;
  x->ExtractView( &xptr );
  u->ExtractView( &uptr );
  y->ExtractView( &yptr );
  v->ExtractView( &vptr );
  a->ExtractView( &aptr );

  // \todo Put in mothership.
  Epetra_Vector lhs(*residual);

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  this->synchDataManagers();

  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  for(int step=0; step<numLoadSteps ; step++){

    double loadIncrement = 1.0/double(numLoadSteps);
    double dt = (timeFinal - timeInitial)*loadIncrement;
    timeCurrent = timeInitial + step*dt;
    *timeStep = dt;
    
    if(peridigmComm->MyPID() == 0)
      cout << "Load step " << step+1 << ", load increment = " << loadIncrement << ", time step = " << dt << ", current time = " << timeCurrent << endl;

    // Update nodal positions for nodes with kinematic B.C.
    deltaU->PutScalar(0.0);
    applyKinematicBC(loadIncrement, deltaU, Teuchos::RCP<Epetra_FECrsMatrix>());

    // Set the current position
    // \todo We probably want to rework this so that the material models get valid x, u, and y values
    // Currently the u values are from the previous load step (and if we update u here we'll be unable to properly undo a time step, which we'll need to adaptive time stepping).
    for(int i=0 ; i<y->MyLength() ; ++i)
      (*y)[i] = (*x)[i] + (*u)[i] + (*deltaU)[i];

    // compute the residual
    double residualNorm = computeQuasiStaticResidual();

    int solverIteration = 1;
    while(residualNorm > absoluteTolerance && solverIteration <= maximumSolverIterations){

      if(peridigmComm->MyPID() == 0)
        cout << "  residual = " << residualNorm << endl;

      // Compute the tangent
      tangent->PutScalar(0.0);
      PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
      modelEvaluator->evalJacobian(workset);
      tangent->GlobalAssemble();
      PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");
      applyKinematicBC(0.0, residual, tangent);
      residual->Scale(-1.0);

      // Solve linear system
      PeridigmNS::Timer::self().startTimer("Solve Linear System");
      Epetra_LinearProblem linearProblem;
      AztecOO solver(linearProblem);
      solver.SetAztecOption(AZ_precond, AZ_Jacobi);
      stringstream ss;
      solver.SetOutputStream(ss);
      int maxAztecIterations = 500;
      double aztecTolerance = 1.0e-6;
      lhs.PutScalar(0.0);
      solver.Iterate(tangent.get(), &lhs, residual.get(), maxAztecIterations, aztecTolerance);
      PeridigmNS::Timer::self().stopTimer("Solve Linear System");

      // Apply increment to nodal positions
      for(int i=0 ; i<y->MyLength() ; ++i)
        (*deltaU)[i] += lhs[i];
      for(int i=0 ; i<y->MyLength() ; ++i)
        (*y)[i] = (*x)[i] + (*u)[i] + (*deltaU)[i];
      
      // Compute residual
      residualNorm = computeQuasiStaticResidual();

      solverIteration++;
    }

    if(peridigmComm->MyPID() == 0)
      cout << "  residual = " << residualNorm << endl;

    // Add the converged displacement increment to the displacement
    for(int i=0 ; i<u->MyLength() ; ++i)
      (*u)[i] += (*deltaU)[i];

    // Write output for completed load step
    PeridigmNS::Timer::self().startTimer("Output");
    this->synchDataManagers();

    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->updateState();

    cout << endl;
  }
}

void PeridigmNS::Peridigm::executeImplicit() {

  // Allocate memory for non-zeros in global Jacobain and lock in the structure
  allocateJacobian();

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  Teuchos::RCP<Teuchos::ParameterList> implicitParams = sublist(solverParams, "Implicit", true);
  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal = solverParams->get("Final Time", 1.0);
  double timeCurrent = timeInitial;
  double absoluteTolerance       = implicitParams->get("Absolute Tolerance", 1.0e-6);
  double maximumSolverIterations = implicitParams->get("Maximum Solver Iterations", 10);
  double dt                      = implicitParams->get("Fixed dt", 1.0);
  double beta                    = implicitParams->get("Beta", 0.25);
  double gamma                   = implicitParams->get("Gamma", 0.50);
  *timeStep = dt;
  double dt2 = dt*dt;
  int nsteps = (int)floor((timeFinal-timeInitial)/dt);

  // Pointer index into sub-vectors for use with BLAS
  double *xptr, *uptr, *yptr, *vptr, *aptr;
  x->ExtractView( &xptr );
  u->ExtractView( &uptr );
  y->ExtractView( &yptr );
  v->ExtractView( &vptr );
  a->ExtractView( &aptr );

  // Data for linear solver object
  Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> linearProblem;
  Teuchos::ParameterList belosList;
  belosList.set( "Block Size", 1 );                                // Use single-vector iteration
  belosList.set( "Maximum Iterations", tangent->NumGlobalRows() ); // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", 1.e-10 );                // Relative convergence tolerance requested
  belosList.set( "Output Frequency", -1 );
  //int verbosity = Belos::Errors + Belos::Warnings + Belos::StatusTestDetails;
  int verbosity = Belos::Errors + Belos::Warnings;
  belosList.set( "Verbosity", verbosity );
  belosList.set( "Output Style", Belos::Brief );
  Teuchos::RCP< Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> > belosSolver
    = Teuchos::rcp( new Belos::BlockCGSolMgr<double,Epetra_MultiVector,Epetra_Operator>(Teuchos::rcp(&linearProblem,false), Teuchos::rcp(&belosList,false)) );

  // Create owned (e.g., "mothership") vectors for data at timestep n
  // MLP: Move this to "Peridigm" level, along with residual vector, deltaU, etc.
  Teuchos::RCP<Epetra_Vector> un = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> vn = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> an = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  // Create temporary owned (e.g., "mothership") vectors for data at timestep n
  // to be used in Newmark integration
  Teuchos::RCP<Epetra_Vector> u2 = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> v2 = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  // \todo Possibly move this functionality into ModelEvaluator.
  // \todo Generalize this for multiple materials
  double density = (*materialModels)[0]->Density();

  // \todo Put in mothership.
  Teuchos::RCP<Epetra_Vector> deltaU = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  this->synchDataManagers();

  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  for(int step=0; step<nsteps ; step++){

    *un = *u;
    *vn = *v;
    *an = *a;

    //u2 = un + dt*vn + 0.5*dt*dt*(1-2*beta)*an
    u2->Update(1.0,*un,0.0);
    u2->Update(dt, *vn, 0.5*dt2*(1.0-2*beta), *an, 1.0);
    //v2 = vn + dt*(1-gamma)*an
    v2->Update(1.0, *vn, dt*(1.0-gamma), *an, 0.0);

    // Fill the owned vectors with probe data
    // Assign predictor (use u2)
    u->Update(1.0,*u2,0.0);
    // a = (1.0/(beta*dt*dt))*(u_np1 - u2);
    //  a will be zero unless a different predictor is used, so do the following computation anyway
    a->Update(1.0, *u, -1.0, *u2, 0.0);
    a->Scale(1.0/(beta*dt2));
    // v = v2 + dt*gamma*an
    v->Update(1.0, *v2, dt*gamma, *a, 0.0);
    // Update y to be consistent with u
    y->Update(1.0, *x, 1.0, *u, 0.0);

    // Copy data from mothership vectors to overlap vectors in data manager
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(*u, Field_NS::DISPL3D,    Field_ENUM::STEP_NP1, Insert);
      blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1, Insert);
      blockIt->importData(*v, Field_NS::VELOC3D,    Field_ENUM::STEP_NP1, Insert);
      blockIt->importData(*v, Field_NS::ACCEL3D,    Field_ENUM::STEP_NP1, Insert);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Model Evaluator");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Model Evaluator");

    // Copy force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->exportData(*force, Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Compute the residual
    // residual = beta*dt*dt*(M*a - force)
    residual->Update(density,*a,0.0); // This computes M*a, since mass matrix is a multiple of the identity
    residual->Update(-1.0,*force,1.0);
    residual->Scale(beta*dt2);

    // Modify residual for kinematic BC
    applyKinematicBC(0.0, residual, Teuchos::RCP<Epetra_FECrsMatrix>());

    double residualNorm;
    residual->Norm2(&residualNorm);

    int NLSolverIteration = 0;
    while(residualNorm > absoluteTolerance && NLSolverIteration <= maximumSolverIterations){

      if(peridigmComm->MyPID() == 0)
        cout << "Time step " << step << ", Newton iteration = " << NLSolverIteration << ", norm(residual) = " << residualNorm << endl;

      // Fill the Jacobian
      computeImplicitJacobian(beta);

      // Modify Jacobian for kinematic BC
      applyKinematicBC(0.0, Teuchos::RCP<Epetra_Vector>(), tangent);

      // Want to solve J*deltaU = -residual
      residual->Scale(-1.0);

      // Solve linear system
      deltaU->PutScalar(0.0);
      linearProblem.setOperator(tangent);
      bool isSet = linearProblem.setProblem(deltaU, residual);
      if (isSet == false) {
        if(peridigmComm->MyPID() == 0)
          std::cout << std::endl << "ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
      }
      PeridigmNS::Timer::self().startTimer("Solve Linear System");
      /* Belos::ReturnType ret = */ belosSolver->solve();
      PeridigmNS::Timer::self().stopTimer("Solve Linear System");

      // Apply increment to nodal positions
      u->Update(1.0,*deltaU,1.0);

      // Update y to be consistent with u
      y->Update(1.0, *x, 1.0, *u, 0.0);

      // a = (1.0/(beta*dt*dt))*(u_np1 - u2);
      a->Update(1.0, *u, -1.0, *u2, 0.0);
      a->Scale(1.0/(beta*dt2));
      // v = v2 + dt*gamma*an
      v->Update(1.0, *v2, dt*gamma, *a, 0.0);

      // Copy data from mothership vectors to overlap vectors in data manager
      PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(*u, Field_NS::DISPL3D,    Field_ENUM::STEP_NP1, Insert);
        blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1, Insert);
        blockIt->importData(*v, Field_NS::VELOC3D,    Field_ENUM::STEP_NP1, Insert);
        blockIt->importData(*v, Field_NS::ACCEL3D,    Field_ENUM::STEP_NP1, Insert);
      }
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Update forces based on new positions
      PeridigmNS::Timer::self().startTimer("Model Evaluator");
      modelEvaluator->evalModel(workset);
      PeridigmNS::Timer::self().stopTimer("Model Evaluator");

      // Copy force from the data manager to the mothership vector
      PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->exportData(*force, Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Compute residual vector and its norm
      // residual = beta*dt*dt*(M*a - force)
      residual->Update(density,*a,0.0); // This computes M*a, since mass matrix is a multiple of the identity
      residual->Update(-1.0,*force,1.0);
      residual->Scale(beta*dt2);

      // Modify residual for kinematic BC
      applyKinematicBC(0.0, residual, Teuchos::RCP<Epetra_FECrsMatrix>());

      residual->Norm2(&residualNorm);

      NLSolverIteration++;
    }

    if(peridigmComm->MyPID() == 0)
      cout << "Time step " << step << ", Newton iteration = " << NLSolverIteration << ", norm(residual) = " << residualNorm << endl;

    timeCurrent = timeInitial + step*dt;

    // Write output for completed time step
    PeridigmNS::Timer::self().startTimer("Output");
    this->synchDataManagers();

    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->updateState();

    cout << endl;
  }
}

void PeridigmNS::Peridigm::allocateJacobian() {

  // Construct map for global tangent matrix
  // Note that this must be an Epetra_Map, not an Epetra_BlockMap, so we can't use threeDimensionalMap directly
  int numGlobalElements = 3*oneDimensionalMap->NumGlobalElements();
  int numMyElements = 3*oneDimensionalMap->NumMyElements();
  int* myGlobalElements = new int[numMyElements];
  int* oneDimensionalMapGlobalElements = oneDimensionalMap->MyGlobalElements();
  for(int iElem=0 ; iElem<oneDimensionalMap->NumMyElements() ; ++iElem){
    myGlobalElements[3*iElem]     = 3*oneDimensionalMapGlobalElements[iElem];
    myGlobalElements[3*iElem + 1] = 3*oneDimensionalMapGlobalElements[iElem] + 1;
    myGlobalElements[3*iElem + 2] = 3*oneDimensionalMapGlobalElements[iElem] + 2;
  }
  int indexBase = 0;
  tangentMap = Teuchos::rcp(new Epetra_Map(numGlobalElements, numMyElements, myGlobalElements, indexBase, *peridigmComm));
  delete[] myGlobalElements;

  // Create the global tangent matrix
  Epetra_DataAccess CV = Copy;
  int numEntriesPerRow = 0;  // If this is zero, allocation will take place during the insertion phase \todo Compute non-zeros instead of allocation during insertion.
  bool staticProfile = false;  // \todo Can staticProfile be set to true?  Bond breaking would alter the non-zeros, but we could just leave them there to avoid reallocation.
  tangent = Teuchos::rcp(new Epetra_FECrsMatrix(CV, *tangentMap, numEntriesPerRow, staticProfile));

  // Loop over the neighborhood for each locally-owned point and create non-zero entries in the matrix
  vector<int> globalIndicies;
  vector<double> zeros;
  int* neighborhoodList = globalNeighborhoodData->NeighborhoodList();
  int neighborhoodListIndex = 0;
  for(int LID=0 ; LID<globalNeighborhoodData->NumOwnedPoints() ; ++LID){
    int GID =  oneDimensionalOverlapMap->GID(LID);
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    unsigned int numEntries = 3*(numNeighbors+1);
    globalIndicies.resize(numEntries);
    globalIndicies[0] = 3*GID;
    globalIndicies[1] = 3*GID + 1;
    globalIndicies[2] = 3*GID + 2;
    for(int j=0 ; j<numNeighbors ; ++j){
      int neighborLocalID = neighborhoodList[neighborhoodListIndex++];
      int neighborGlobalID = oneDimensionalOverlapMap->GID(neighborLocalID);
      globalIndicies[3*j+3] = 3*neighborGlobalID;
      globalIndicies[3*j+4] = 3*neighborGlobalID + 1;
      globalIndicies[3*j+5] = 3*neighborGlobalID + 2;
    }
    if(numEntries > zeros.size())
      zeros.resize(numEntries, 0.0);
    tangent->InsertGlobalValues(3*GID,   numEntries, &zeros[0], &globalIndicies[0]); 
    tangent->InsertGlobalValues(3*GID+1, numEntries, &zeros[0], &globalIndicies[0]); 
    tangent->InsertGlobalValues(3*GID+2, numEntries, &zeros[0], &globalIndicies[0]); 
  }
  tangent->GlobalAssemble();

  // create the serial Jacobian
  overlapJacobian = Teuchos::rcp(new PeridigmNS::SerialMatrix(tangent, oneDimensionalOverlapMap));
  workset->jacobian = overlapJacobian;
}

void PeridigmNS::Peridigm::applyKinematicBC(double loadIncrement,
                                            Teuchos::RCP<Epetra_Vector> vec,
                                            Teuchos::RCP<Epetra_FECrsMatrix> mat) {
  PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");

  Teuchos::ParameterList& problemParams = peridigmParams->sublist("Problem");
  Teuchos::ParameterList& bcParams = problemParams.sublist("Boundary Conditions");
  Teuchos::ParameterList::ConstIterator it;

  // create data structures for inserting ones and zeros into jacobian
  vector<double> jacobianRow;
  vector<int> jacobianIndicies;
  if(!mat.is_null()){
    jacobianRow.resize(mat->NumMyCols(), 0.0);
    jacobianIndicies.resize(mat->NumMyCols());
    for(unsigned int i=0 ; i<jacobianIndicies.size() ; ++i)
      jacobianIndicies[i] = i;
  }

  // apply the kinematic boundary conditions
  for(it = bcParams.begin() ; it != bcParams.end() ; it++){
    const string & name = it->first;
    size_t position = name.find("Prescribed Displacement");
    if(position != string::npos){
      Teuchos::ParameterList & boundaryConditionParams = Teuchos::getValue<Teuchos::ParameterList>(it->second);
      string nodeSet = boundaryConditionParams.get<string>("Node Set");
      string type = boundaryConditionParams.get<string>("Type");
      string coordinate = boundaryConditionParams.get<string>("Coordinate");
      double value = boundaryConditionParams.get<double>("Value");

      int coord = 0;
      if(coordinate == "y" || coordinate == "Y")
        coord = 1;
      if(coordinate == "z" || coordinate == "Z")
        coord = 2;

      // apply kinematic boundary conditions to locally-owned nodes
      TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSet) == nodeSets->end(), "**** Node set not found: " + name + "\n");
      vector<int> & nodeList = (*nodeSets)[nodeSet];
      for(unsigned int i=0 ; i<nodeList.size() ; i++){
        // zero out the row and column and put a 1.0 on the diagonal
        if(!mat.is_null()){
          int globalID = 3*nodeList[i] + coord;
          int localRowID = mat->LRID(globalID);
          int localColID = mat->LCID(globalID);

          // zero out all locally-owned entries in the column associated with this dof
          // \todo Call ReplaceMyValues only for entries that actually exist in the matrix structure.
          double zero = 0.0;
          for(int iRow=0 ; iRow<mat->NumMyRows() ; ++iRow)
            mat->ReplaceMyValues(iRow, 1, &zero, &localColID);

          // zero out the row and put a 1.0 on the diagonal
          if(localRowID != -1){
            jacobianRow[localColID] = 1.0;
            // From Epetra_CrsMatrix documentation:
            // If a value is not already present for the specified location in the matrix, the
            // input value will be ignored and a positive warning code will be returned.
            // \todo Do the bookkeeping to send in data only for locations that actually exist in the matrix structure.
            mat->ReplaceMyValues(localRowID, mat->NumMyCols(), &jacobianRow[0], &jacobianIndicies[0]);
            jacobianRow[localColID] = 0.0;
          }
        }

        // set entry in residual vector equal to the displacement increment for the kinematic bc
        // this will cause the solution procedure to solve for the correct U at the bc
        int localNodeID = threeDimensionalMap->LID(nodeList[i]);
	if(!vec.is_null() && localNodeID != -1){
 	  (*vec)[localNodeID*3 + coord] = value*loadIncrement;
        }

      }

    }
  }
  PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");
}

double PeridigmNS::Peridigm::computeQuasiStaticResidual() {

  PeridigmNS::Timer::self().startTimer("Compute Residual");

  // The residual is computed as the L2 norm of the internal force vector with the
  // entries corresponding to kinematic BC zeroed out.

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, Field_NS::DISPL3D,    Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*v, Field_NS::VELOC3D,    Field_ENUM::STEP_NP1, Insert);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // Update forces based on new positions
  PeridigmNS::Timer::self().startTimer("Model Evaluator");
  modelEvaluator->evalModel(workset);
  PeridigmNS::Timer::self().stopTimer("Model Evaluator");

  // Copy force from the data manager to the mothership vector
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->exportData(*force, Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Add);
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // copy the internal force to the residual vector
  (*residual) = (*force);
    
  // zero out the rows corresponding to kinematic boundary conditions and compute the residual
  applyKinematicBC(0.0, residual, Teuchos::RCP<Epetra_FECrsMatrix>());
  double residualNorm2;
  residual->Norm2(&residualNorm2);

  PeridigmNS::Timer::self().stopTimer("Compute Residual");

  return residualNorm2;
}

void PeridigmNS::Peridigm::computeImplicitJacobian(double beta) {

  // MLP: Pass this in from integrator routine
  double dt = *(workset->timeStep);
  double dt2 = dt*dt;

  // Compute the tangent
  tangent->PutScalar(0.0);
  PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
  modelEvaluator->evalJacobian(workset);
  tangent->GlobalAssemble();
  PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");
//      applyKinematicBC(0.0, residual, tangent);

  // Code to symmeterize Jacobian

  // First, construct transpose
  bool makeDataContiguous = true;
  EpetraExt::RowMatrix_Transpose transposer( makeDataContiguous );
  Epetra_CrsMatrix & transTangent = dynamic_cast<Epetra_CrsMatrix&>(transposer(*tangent));

  // Now loop over all owned rows and average entries of both
  int numRows      = tangent->NumMyRows();
  int numRowsTrans = transTangent.NumMyRows();
  if (numRows != numRowsTrans) cout << "Number of rows mismatch!" << std::endl;
  int numEntries, numEntriesTrans;
  double *values, *valuesTrans;
  int *indices, *indicesTrans;
  // Assert here that numRows == numRowsTrans
  for(int i=0;i<numRows;i++) {
    tangent->ExtractMyRowView(i, numEntries, values, indices);
    transTangent.ExtractMyRowView(i, numEntriesTrans, valuesTrans, indicesTrans);
    if (numEntries != numEntriesTrans) cout << "Number of entries mismatch!" << std::endl;
    for (int j=0;j<numEntries;j++) {
      values[j] = 0.5*(values[j]+valuesTrans[j]);
      if (indices[j] != indicesTrans[j]) cout << "index mismatch!" << std::endl;
    }
  }

  // Now add in mass matrix contribution
  // \todo Generalize this for multiple materials
  double density = (*materialModels)[0]->Density();

  // tangent = M - beta*dt*dt*K
  tangent->Scale(-beta*dt2);

  Epetra_Vector diagonal1(tangent->RowMap());
  Epetra_Vector diagonal2(tangent->RowMap());
  tangent->ExtractDiagonalCopy(diagonal1);
  diagonal2.PutScalar(density);
  diagonal1.Update(1.0, diagonal2, 1.0);
  tangent->ReplaceDiagonalValues(diagonal1);

}

void PeridigmNS::Peridigm::synchDataManagers() {
  // Need to ensure these primal fields are synchronized: VOLUME, COORD3D, DISPL3D, CURCOORD3D, VELOC3D, FORCE_DENSITY3D, CONTACT_FORCE_DENSITY_3D

  // Copy data from mothership vectors to overlap vectors in blocks
  // VOLUME is synched during creation and rebalance, and otherwise never changes
  // COORD3D is synched during creation and rebalance, and otherwise never changes
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, Field_NS::DISPL3D, Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*y, Field_NS::CURCOORD3D, Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*v, Field_NS::VELOC3D, Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*force, Field_NS::FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Insert);
    blockIt->importData(*contactForce, Field_NS::CONTACT_FORCE_DENSITY3D, Field_ENUM::STEP_NP1, Insert);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
}

void FalseDeleter(Epetra_Comm* c) {}

void PeridigmNS::Peridigm::rebalance() {

  // \todo Handle serial case.  We don't need to rebalance, but we still want to update the contact search.
  QUICKGRID::Data rebalancedDecomp = currentConfigurationDecomp();

  Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGridDiscretization::getOwnedMap(*peridigmComm, rebalancedDecomp, 1)));
  Teuchos::RCP<const Epetra_Import> oneDimensionalMapImporter = Teuchos::rcp(new Epetra_Import(*rebalancedOneDimensionalMap, *oneDimensionalMap));

  Teuchos::RCP<Epetra_BlockMap> rebalancedThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGridDiscretization::getOwnedMap(*peridigmComm, rebalancedDecomp, 3)));
  Teuchos::RCP<const Epetra_Import> threeDimensionalMapImporter = Teuchos::rcp(new Epetra_Import(*rebalancedThreeDimensionalMap, *threeDimensionalMap));

  Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap = createRebalancedBondMap(rebalancedOneDimensionalMap, oneDimensionalMapImporter);
  Teuchos::RCP<const Epetra_Import> bondMapImporter = Teuchos::rcp(new Epetra_Import(*rebalancedBondMap, *bondMap));

  // create a list of neighbors in the rebalanced configuration
  // this list has the global ID for each neighbor of each on-processor point (that is, on processor in the rebalanced configuration)
  Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs = createRebalancedNeighborGlobalIDList(rebalancedBondMap, bondMapImporter);

  // create a list of all the off-processor IDs that will need to be ghosted
  // \todo Use set::reserve() for better memory allocation here.
  set<int> offProcessorIDs;
  for(int i=0 ; i<rebalancedNeighborGlobalIDs->MyLength() ; ++i){
    int globalID = (int)( (*rebalancedNeighborGlobalIDs)[i] );
    if(!rebalancedOneDimensionalMap->MyGID(globalID))
      offProcessorIDs.insert(globalID);
  }

  // this function does three things:
  // 1) fills the neighborhood information in rebalancedDecomp based on the contact search
  // 2) creates a list of global IDs for each locally-owned point that will need to be searched for contact (contactNeighborGlobalIDs)
  // 3) keeps track of the additional off-processor IDs that need to be ghosted as a result of the contact search (offProcessorContactIDs)
  Teuchos::RCP< map<int, vector<int> > > contactNeighborGlobalIDs = Teuchos::rcp(new map<int, vector<int> >());
  Teuchos::RCP< set<int> > offProcessorContactIDs = Teuchos::rcp(new set<int>());
  if(analysisHasContact)
    contactSearch(rebalancedOneDimensionalMap, rebalancedBondMap, rebalancedNeighborGlobalIDs, rebalancedDecomp, contactNeighborGlobalIDs, offProcessorContactIDs);

  // add the off-processor IDs required for contact to the list of points that will be ghosted
  for(set<int>::const_iterator it=offProcessorContactIDs->begin() ; it!=offProcessorContactIDs->end() ; it++){
    offProcessorIDs.insert(*it);
  }

  // construct the rebalanced overlap maps
  int numGlobalElements = -1;
  int numMyElements = rebalancedOneDimensionalMap->NumMyElements() + offProcessorIDs.size();
  int* myGlobalElements = new int[numMyElements];
  rebalancedOneDimensionalMap->MyGlobalElements(myGlobalElements);
  int offset = rebalancedOneDimensionalMap->NumMyElements();
  int index = 0;
  for(set<int>::const_iterator it=offProcessorIDs.begin() ; it!=offProcessorIDs.end() ; ++it, ++index){
    myGlobalElements[offset+index] = *it;
  }
  int indexBase = 0;
  Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalOverlapMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, *peridigmComm));
  Teuchos::RCP<Epetra_BlockMap> rebalancedThreeDimensionalOverlapMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 3, indexBase, *peridigmComm));
  delete[] myGlobalElements;

  // create a new NeighborhoodData object
  Teuchos::RCP<PeridigmNS::NeighborhoodData> rebalancedNeighborhoodData = createRebalancedNeighborhoodData(rebalancedOneDimensionalMap,
                                                                                                           rebalancedOneDimensionalOverlapMap,
                                                                                                           rebalancedBondMap,
                                                                                                           rebalancedNeighborGlobalIDs);

  // create a new NeighborhoodData object for contact
  Teuchos::RCP<PeridigmNS::NeighborhoodData> rebalancedContactNeighborhoodData;
  if(analysisHasContact)
    rebalancedContactNeighborhoodData = createRebalancedContactNeighborhoodData(contactNeighborGlobalIDs,
                                                                                rebalancedOneDimensionalMap,
                                                                                rebalancedOneDimensionalOverlapMap);

  // rebalance the global vectors (stored in the mothership multivectors)

  Teuchos::RCP<Epetra_MultiVector> rebalancedOneDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*rebalancedOneDimensionalMap, oneDimensionalMothership->NumVectors()));
  rebalancedOneDimensionalMothership->Import(*oneDimensionalMothership, *oneDimensionalMapImporter, Insert);
  oneDimensionalMothership = rebalancedOneDimensionalMothership;
  blockIDs = Teuchos::rcp((*oneDimensionalMothership)(0), false);        // block ID
  volume = Teuchos::rcp((*oneDimensionalMothership)(1), false);          // cell volume

  Teuchos::RCP<Epetra_MultiVector> rebalancedThreeDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*rebalancedThreeDimensionalMap, threeDimensionalMothership->NumVectors()));
  rebalancedThreeDimensionalMothership->Import(*threeDimensionalMothership, *threeDimensionalMapImporter, Insert);
  threeDimensionalMothership = rebalancedThreeDimensionalMothership;
  x = Teuchos::rcp((*threeDimensionalMothership)(0), false);             // initial positions
  u = Teuchos::rcp((*threeDimensionalMothership)(1), false);             // displacement
  y = Teuchos::rcp((*threeDimensionalMothership)(2), false);             // current positions
  v = Teuchos::rcp((*threeDimensionalMothership)(3), false);             // velocities
  a = Teuchos::rcp((*threeDimensionalMothership)(4), false);             // accelerations
  force = Teuchos::rcp((*threeDimensionalMothership)(5), false);         // force
  contactForce = Teuchos::rcp((*threeDimensionalMothership)(6), false);  // contact force (used only for contact simulations)
  deltaU = Teuchos::rcp((*threeDimensionalMothership)(7), false);        // increment in displacement (used only for implicit time integration)
  residual = Teuchos::rcp((*threeDimensionalMothership)(8), false);      // residual (used only for implicit time integration)
  scratch = Teuchos::rcp((*threeDimensionalMothership)(9), false);       // scratch space

  // rebalance the blocks
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->rebalance(rebalancedOneDimensionalMap,
                       rebalancedOneDimensionalOverlapMap,
                       rebalancedThreeDimensionalMap,
                       rebalancedThreeDimensionalOverlapMap,
                       rebalancedBondMap,
                       blockIDs,
                       rebalancedNeighborhoodData,
                       rebalancedContactNeighborhoodData);

  // Initialize what we can for newly-created ghosts across material boundaries.
  // These ghosts are added due to contact and will only be used by the contact evaluator.
  // There is a problem here in that history data and material-specific data are not
  // going to be available for these newly-created ghosts across material boundaries.
  // As long as the contact algorithm only uses mothership data at STATE_NONE and
  // STATE_NP1, then we can get away with the scatter operation below.  If history data
  // is needed, then some additional MPI operations are needed.  Better yet, this issue
  // could be resolved if contact were refactored to be totally separate from the material
  // models (i.e., give contact its own mothership vectors and data managers, and rebalance
  // only these objects when executing a contact search).
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->importData(*volume, Field_NS::VOLUME, Field_ENUM::STEP_NONE, Insert);

  // set all the pointers to the new maps
  oneDimensionalMap = rebalancedOneDimensionalMap;
  oneDimensionalOverlapMap = rebalancedOneDimensionalOverlapMap;
  threeDimensionalMap = rebalancedThreeDimensionalMap;
  bondMap = rebalancedBondMap;

  // update neighborhood data
  globalNeighborhoodData = rebalancedNeighborhoodData;
  globalContactNeighborhoodData = rebalancedContactNeighborhoodData;
}

QUICKGRID::Data PeridigmNS::Peridigm::currentConfigurationDecomp() {

  // Create a decomp object and fill necessary data for rebalance
  int myNumElements = oneDimensionalMap->NumMyElements();
  int dimension = 3;
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(myNumElements, dimension);

  decomp.globalNumPoints = oneDimensionalMap->NumGlobalElements();

  // fill myGlobalIDs
  Array<int> myGlobalIDs(myNumElements);
  int* myGlobalIDsPtr = myGlobalIDs.get();
  int* gIDs = oneDimensionalMap->MyGlobalElements();
  memcpy(myGlobalIDsPtr, gIDs, myNumElements*sizeof(int));
  decomp.myGlobalIDs = myGlobalIDs.get_shared_ptr();

  // fill myX
  // use current positions for x
  Array<double> myX(myNumElements*dimension);
  double* myXPtr = myX.get();
  double* yPtr;
  y->ExtractView(&yPtr);
  memcpy(myXPtr, yPtr, myNumElements*dimension*sizeof(double));
  decomp.myX = myX.get_shared_ptr();

  // fill cellVolume
  Array<double> cellVolume(myNumElements);
  double* cellVolumePtr = cellVolume.get();
  double* volumePtr;
  volume->ExtractView(&volumePtr);
  memcpy(cellVolumePtr, volumePtr, myNumElements*sizeof(double));
  decomp.cellVolume = cellVolume.get_shared_ptr();

  // call the rebalance function on the current-configuration decomp
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

  return decomp;
}

Teuchos::RCP<Epetra_BlockMap> PeridigmNS::Peridigm::createRebalancedBondMap(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                            Teuchos::RCP<const Epetra_Import> oneDimensionalMapToRebalancedOneDimensionalMapImporter) {

  // communicate the number of bonds for each point so that space for bond data can be allocated
  Teuchos::RCP<Epetra_Vector> numberOfBonds = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  for(int i=0 ; i<oneDimensionalMap->NumMyElements() ; ++i){
    int globalID = oneDimensionalMap->GID(i);
    int bondMapLocalID = bondMap->LID(globalID);
    if(bondMapLocalID != -1)
      (*numberOfBonds)[i] = (double)( bondMap->ElementSize(i) );
  }
  Teuchos::RCP<Epetra_Vector> rebalancedNumberOfBonds = Teuchos::rcp(new Epetra_Vector(*rebalancedOneDimensionalMap));
  rebalancedNumberOfBonds->Import(*numberOfBonds, *oneDimensionalMapToRebalancedOneDimensionalMapImporter, Insert);

  // create the rebalanced bond map
  // care must be taken because you cannot have an element with zero length
  int numMyElementsUpperBound = rebalancedOneDimensionalMap->NumMyElements();
  int numGlobalElements = -1; 
  int numMyElements = 0;
  int* rebalancedOneDimensionalMapGlobalElements = rebalancedOneDimensionalMap->MyGlobalElements();
  int* myGlobalElements = new int[numMyElementsUpperBound];
  int* elementSizeList = new int[numMyElementsUpperBound];
  int numPointsWithZeroNeighbors = 0;
  for(int i=0 ; i<numMyElementsUpperBound ; ++i){
    int numBonds = (int)( (*rebalancedNumberOfBonds)[i] );
    if(numBonds > 0){
      numMyElements++;
      myGlobalElements[i-numPointsWithZeroNeighbors] = rebalancedOneDimensionalMapGlobalElements[i];
      elementSizeList[i-numPointsWithZeroNeighbors] = numBonds;
    }
    else{
      numPointsWithZeroNeighbors++;
    }
  }
  int indexBase = 0;
  Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap = 
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, *peridigmComm));
  delete[] myGlobalElements;
  delete[] elementSizeList;

  return rebalancedBondMap;
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::Peridigm::createRebalancedNeighborGlobalIDList(Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                                       Teuchos::RCP<const Epetra_Import> bondMapToRebalancedBondMapImporter) {

  // construct a globalID neighbor list for the current decomposition
  Teuchos::RCP<Epetra_Vector> neighborGlobalIDs = Teuchos::rcp(new Epetra_Vector(*bondMap));
  int* neighborhoodList = globalNeighborhoodData->NeighborhoodList();
  int neighborhoodListIndex = 0;
  int neighborGlobalIDIndex = 0;
  for(int i=0 ; i<globalNeighborhoodData->NumOwnedPoints() ; ++i){
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(int j=0 ; j<numNeighbors ; ++j){
      int neighborLocalID = neighborhoodList[neighborhoodListIndex++];
      (*neighborGlobalIDs)[neighborGlobalIDIndex++] = oneDimensionalOverlapMap->GID(neighborLocalID);
    }
  }

  // redistribute the globalID neighbor list to the rebalanced configuration
  Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs = Teuchos::rcp(new Epetra_Vector(*rebalancedBondMap));
  rebalancedNeighborGlobalIDs->Import(*neighborGlobalIDs, *bondMapToRebalancedBondMapImporter, Insert);

  return rebalancedNeighborGlobalIDs;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> PeridigmNS::Peridigm::createRebalancedNeighborhoodData(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                                  Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalOverlapMap,
                                                                                                  Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                                                  Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs) {

  Teuchos::RCP<PeridigmNS::NeighborhoodData> rebalancedNeighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
  rebalancedNeighborhoodData->SetNumOwned(rebalancedOneDimensionalMap->NumMyElements());
  int* ownedIDs = rebalancedNeighborhoodData->OwnedIDs();
  for(int i=0 ; i<rebalancedOneDimensionalMap->NumMyElements() ; ++i){
    int globalID = rebalancedOneDimensionalMap->GID(i);
    int localID = rebalancedOneDimensionalOverlapMap->LID(globalID);
    TEST_FOR_EXCEPTION(localID == -1, Teuchos::RangeError, "Invalid index into rebalancedOneDimensionalOverlapMap");
    ownedIDs[i] = localID;
  }
  rebalancedNeighborhoodData->SetNeighborhoodListSize(rebalancedOneDimensionalMap->NumMyElements() + rebalancedBondMap->NumMyPoints());
  // numNeighbors1, n1LID, n2LID, n3LID, numNeighbors2, n1LID, n2LID, ...
  int* neighborhoodList = rebalancedNeighborhoodData->NeighborhoodList();
  // points into neighborhoodList, gives start of neighborhood information for each locally-owned element
  int* neighborhoodPtr = rebalancedNeighborhoodData->NeighborhoodPtr();
  // gives the offset at which the list of neighbors can be found in the rebalancedNeighborGlobalIDs vector for each locally-owned element
  int* firstPointInElementList = rebalancedBondMap->FirstPointInElementList();
  // loop over locally owned points
  int neighborhoodIndex = 0;
  for(int iLID=0 ; iLID<rebalancedOneDimensionalMap->NumMyElements() ; ++iLID){
    // location of this element's neighborhood data in the neighborhoodList
    neighborhoodPtr[iLID] = neighborhoodIndex;
    // first entry is the number of neighbors
    int globalID = rebalancedOneDimensionalMap->GID(iLID);
    int rebalancedBondMapLocalID = rebalancedBondMap->LID(globalID);
    if(rebalancedBondMapLocalID != -1){
      int numNeighbors = rebalancedBondMap->ElementSize(rebalancedBondMapLocalID);
      neighborhoodList[neighborhoodIndex++] = numNeighbors;
      // next entries record the local ID of each neighbor
      int offset = firstPointInElementList[rebalancedBondMapLocalID];
      for(int iN=0 ; iN<numNeighbors ; ++iN){
        int globalNeighborID = (int)( (*rebalancedNeighborGlobalIDs)[offset + iN] );
        int localNeighborID = rebalancedOneDimensionalOverlapMap->LID(globalNeighborID);
        TEST_FOR_EXCEPTION(localNeighborID == -1, Teuchos::RangeError, "Invalid index into rebalancedOneDimensionalOverlapMap");
        neighborhoodList[neighborhoodIndex++] = localNeighborID;
      }
    }
    else{
      neighborhoodList[neighborhoodIndex++] = 0;
    }
  }

  return rebalancedNeighborhoodData;
}

std::vector<Field_NS::FieldSpec> PeridigmNS::Peridigm::getFieldSpecs() {

  // FieldSpecs used by Peridigm class
  std::vector<Field_NS::FieldSpec> mySpecs;

  mySpecs.push_back(Field_NS::VOLUME);
  mySpecs.push_back(Field_NS::COORD3D);
  mySpecs.push_back(Field_NS::DISPL3D);
  mySpecs.push_back(Field_NS::CURCOORD3D);
  mySpecs.push_back(Field_NS::VELOC3D);
  mySpecs.push_back(Field_NS::FORCE_DENSITY3D);
  mySpecs.push_back(Field_NS::CONTACT_FORCE_DENSITY3D);

  return mySpecs;
}

void PeridigmNS::Peridigm::displayProgress(string title, double percentComplete){
  
  int numMarks = 20;

  static bool firstCall = true;

  if(firstCall)
    firstCall = false;
  else
    for(int i=0 ; i<(numMarks + 16 + (int)title.size()) ; ++i) *out << "\b";

  percentComplete = double(int(percentComplete));

  stringstream ssPercentComplete;
  ssPercentComplete << int(percentComplete);
  int stringSize = ssPercentComplete.str().size();

  stringstream ss;
  ss << title << " [";
  for(int i=0 ; i<int(percentComplete*(numMarks+2-stringSize)/100.0) ; i++)
    ss << "=";
  ss << ssPercentComplete.str() << "% Complete";
  for(int i=int(percentComplete*(numMarks+2-stringSize)/100.0) ; i<(numMarks+2-stringSize) ; i++)
    ss << "=";
  ss << "] ";

  *out << ss.str();
}

template<class T>
struct NonDeleter{
	void operator()(T* d) {}
};

void PeridigmNS::Peridigm::contactSearch(Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap, 
                                         Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap,
                                         Teuchos::RCP<const Epetra_Vector> rebalancedNeighborGlobalIDs,
                                         QUICKGRID::Data& rebalancedDecomp,
                                         Teuchos::RCP< map<int, vector<int> > > contactNeighborGlobalIDs,
                                         Teuchos::RCP< set<int> > offProcessorContactIDs)
{
// execute contact search

//	rebalancedDecomp = createAndAddNeighborhood(rebalancedDecomp, contactSearchRadius);
	shared_ptr<const Epetra_Comm> comm(peridigmComm.getRawPtr(),NonDeleter<const Epetra_Comm>());
	QUICKGRID::Data d = rebalancedDecomp;
	PDNEIGH::NeighborhoodList neighList(comm,d.zoltanPtr.get(),d.numPoints,d.myGlobalIDs,d.myX,contactSearchRadius);


//	int* searchNeighborhood = rebalancedDecomp.neighborhood.get();
	int* searchNeighborhood = neighList.get_neighborhood().get();

//	int* searchGlobalIDs = rebalancedDecomp.myGlobalIDs.get();
	int* searchGlobalIDs = neighList.get_owned_gids().get();
	int searchListIndex = 0;
	for(size_t iPt=0 ; iPt<rebalancedDecomp.numPoints ; ++iPt){

		int globalID = searchGlobalIDs[iPt];
		vector<int>& contactNeighborGlobalIDList = (*contactNeighborGlobalIDs)[globalID];

		// create a stl::list of global IDs that this point is bonded to
		list<int> bondedNeighbors;
		int tempLocalID = rebalancedBondMap->LID(globalID);
		// if there is no entry in rebalancedBondMap, then there are no bonded neighbors for this point
		if(tempLocalID != -1){
			int firstNeighbor = rebalancedBondMap->FirstPointInElementList()[tempLocalID];
			int numNeighbors = rebalancedBondMap->ElementSize(tempLocalID);
			for(int i=0 ; i<numNeighbors ; ++i){
				int neighborGlobalID = (int)( (*rebalancedNeighborGlobalIDs)[firstNeighbor + i] );
				bondedNeighbors.push_back(neighborGlobalID);
			}
		}

		// loop over the neighbors found by the contact search
		// retain only those neighbors that are not bonded
		int searchNumNeighbors = searchNeighborhood[searchListIndex++];
		for(int iNeighbor=0 ; iNeighbor<searchNumNeighbors ; ++iNeighbor){
			int globalNeighborID = searchNeighborhood[searchListIndex++];
			list<int>::iterator it = find(bondedNeighbors.begin(), bondedNeighbors.end(), globalNeighborID);
			if(it == bondedNeighbors.end()){
				contactNeighborGlobalIDList.push_back(globalNeighborID);
				if(rebalancedOneDimensionalMap->LID(globalNeighborID) == -1)
					offProcessorContactIDs->insert(globalNeighborID);
			}
		}
	}
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> PeridigmNS::Peridigm::createRebalancedContactNeighborhoodData(Teuchos::RCP<map<int, vector<int> > > contactNeighborGlobalIDs,
		Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap,
		Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalOverlapMap) {

	Teuchos::RCP<PeridigmNS::NeighborhoodData> rebalancedContactNeighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
	// record the owned IDs
	rebalancedContactNeighborhoodData->SetNumOwned(rebalancedOneDimensionalMap->NumMyElements());
	int* ownedIDs = rebalancedContactNeighborhoodData->OwnedIDs();
	for(int i=0 ; i<rebalancedOneDimensionalMap->NumMyElements() ; ++i){
		int globalID = rebalancedOneDimensionalMap->GID(i);
		int localID = rebalancedOneDimensionalOverlapMap->LID(globalID);
		TEST_FOR_EXCEPTION(localID == -1, Teuchos::RangeError, "Invalid index into rebalancedOneDimensionalOverlapMap");
		ownedIDs[i] = localID;
	}
	// determine the neighborhood list size
	int neighborhoodListSize = 0;
	for(map<int, vector<int> >::const_iterator it=contactNeighborGlobalIDs->begin() ; it!=contactNeighborGlobalIDs->end() ; it++)
		neighborhoodListSize += it->second.size() + 1;
	rebalancedContactNeighborhoodData->SetNeighborhoodListSize(neighborhoodListSize);
	// numNeighbors1, n1LID, n2LID, n3LID, numNeighbors2, n1LID, n2LID, ...
	int* neighborhoodList = rebalancedContactNeighborhoodData->NeighborhoodList();
	// points into neighborhoodList, gives start of neighborhood information for each locally-owned element
	int* neighborhoodPtr = rebalancedContactNeighborhoodData->NeighborhoodPtr();
	// loop over locally owned points
	int neighborhoodIndex = 0;
	for(int iLID=0 ; iLID<rebalancedOneDimensionalMap->NumMyElements() ; ++iLID){
		// location of this element's neighborhood data in the neighborhoodList
		neighborhoodPtr[iLID] = neighborhoodIndex;
		// get the global ID of this point and the global IDs of its neighbors
		int globalID = rebalancedOneDimensionalMap->GID(iLID);
		// require that this globalID be present as a key into contactNeighborGlobalIDs
		TEST_FOR_EXCEPTION(contactNeighborGlobalIDs->count(globalID) == 0, Teuchos::RangeError, "Invalid index into contactNeighborGlobalIDs");
		const vector<int>& neighborGlobalIDs = (*contactNeighborGlobalIDs)[globalID];
		// first entry in the neighborhoodlist is the number of neighbors
		neighborhoodList[neighborhoodIndex++] = (int) neighborGlobalIDs.size();
		// next entries record the local ID of each neighbor
		for(unsigned int iNeighbor=0 ; iNeighbor<neighborGlobalIDs.size() ; ++iNeighbor){
			neighborhoodList[neighborhoodIndex++] = rebalancedOneDimensionalOverlapMap->LID( neighborGlobalIDs[iNeighbor] );
		}
	}

	return rebalancedContactNeighborhoodData;
}


