/*! \file Peridigm.cpp
 *
 * File containing main class for Peridigm: A parallel, multi-physics,
 * peridynamics simulation code.
 */

//@HEADER
// ************************************************************************
//
//
//                             Peridigm
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met
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
//

#include <iostream>
#include <vector>
#include <map>

#include <boost/math/special_functions/fpclassify.hpp>

#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <EpetraExt_BlockMapOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Epetra_RowMatrixTransposer.h>
#include <Ifpack.h>
#include <Ifpack_IC.h>
#include <Teuchos_VerboseObject.hpp>

#include "Peridigm.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "Peridigm_PdQuickGridDiscretization.hpp"
#include "Peridigm_PartialVolumeCalculator.hpp"
#include "Peridigm_OutputManager_ExodusII.hpp"
#include "Peridigm_ComputeManager.hpp"
#include "Peridigm_BoundaryAndInitialConditionManager.hpp"
#include "Peridigm_CriticalTimeStep.hpp"
#include "Peridigm_RandomNumber.hpp"
#include "Peridigm_Timer.hpp"
#include "materials/Peridigm_MaterialFactory.hpp"
#include "damage/Peridigm_DamageModelFactory.hpp"
#include "contact/Peridigm_ContactModelFactory.hpp"
#include "mesh_input/quick_grid/QuickGrid.h"
#include "mesh_input/quick_grid/QuickGridData.h"
#include "pdneigh/PdZoltan.h"
#include "pdneigh/NeighborhoodList.h"
#include "muParser/muParser.h"
#include "muParser/muParserPeridigmFunctions.h"


using namespace std;

PeridigmNS::Peridigm::Peridigm(const Teuchos::RCP<const Epetra_Comm>& comm,
                   const Teuchos::RCP<Teuchos::ParameterList>& params)
  : numBlocks(1),
    analysisHasRebalance(false),
    rebalanceFrequency(1),
    analysisHasContact(false),
    contactRebalanceFrequency(0),
    contactSearchRadius(0.0),
    analysisHasPartialVolumes(false),
    blockIdFieldId(-1),
    volumeFieldId(-1),
    modelCoordinatesFieldId(-1),
    coordinatesFieldId(-1),
    displacementFieldId(-1),
    velocityFieldId(-1),
    accelerationFieldId(-1),
    forceDensityFieldId(-1),
    contactForceDensityFieldId(-1),
    partialVolumeFieldId(-1),
    tangentReferenceCoordinatesFieldId(-1)
{
  peridigmComm = comm;
  peridigmParams = params;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // Seed random number generator for reproducable results
  seed_rand_num( 42 );

  // Read mesh from disk or generate using geometric primatives.
  Teuchos::RCP<Teuchos::ParameterList> discParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Discretization", true) );

  // \todo When using partial volumes, the horizon should be increased by a value equal to the largest element dimension in the model (largest element diagonal for hexes).
  //       For an initial test, just double the horizon and hope for the best.
  if(analysisHasPartialVolumes)
    discParams->set("Search Horizon", 2.0*discParams->get<double>("Horizon"));
  else
    discParams->set("Search Horizon", discParams->get<double>("Horizon"));

  DiscretizationFactory discFactory(discParams);
  Teuchos::RCP<AbstractDiscretization> peridigmDisc = discFactory.create(peridigmComm);
  initializeDiscretization(peridigmDisc);

  // If the user did not provide a finite-difference probe length, use a fraction of the minimum element radius
  double minElementRadius = peridigmDisc->getMinElementRadius();
  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver");
  if(!solverParams->isParameter("Finite Difference Probe Length"))
    solverParams->set("Finite Difference Probe Length", 1.0e-6*minElementRadius);

  // Load node sets from input deck and/or input mesh file into nodeSets container
  Teuchos::RCP<Teuchos::ParameterList> bcParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Boundary Conditions") );
  initializeNodeSets(bcParams, peridigmDisc);

  boundaryAndInitialConditionManager =
    Teuchos::RCP<BoundaryAndInitialConditionManager>(new BoundaryAndInitialConditionManager(*bcParams));

  boundaryAndInitialConditionManager->initialize(peridigmDisc);

  // Instantiate material models
  instantiateMaterials();

  // Instantiate damage models
  if(params->isSublist("Damage Models"))
     instantiateDamageModels();

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  elementIdFieldId                   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id");
  blockIdFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Block_Id");
  volumeFieldId                      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  modelCoordinatesFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  coordinatesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  displacementFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement");
  velocityFieldId                    = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  accelerationFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration");
  forceDensityFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  contactForceDensityFieldId         = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density");
  tangentReferenceCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Tangent_Reference_Coordinates");
  if(analysisHasPartialVolumes)
    partialVolumeFieldId             = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::CONSTANT, "Partial_Volume");

  // Create field ids that may be required for output
  fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num");

  if(analysisHasPartialVolumes)
    contactForceDensityFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Partial_Volume");

  // Instantiate compute manager
  instantiateComputeManager();

  // Instantiate the blocks
  initializeBlocks(peridigmDisc);

  // Associate material models and damage models with blocks
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    // Set the material model
    string materialModelName = blockIt->getMaterialName();
    std::map< std::string, Teuchos::RCP<const PeridigmNS::Material> >::iterator materialModelIt;
    materialModelIt = materialModels.find(materialModelName);
    if(materialModelIt == materialModels.end()){
      string msg = "\n**** Error, invalid material model:  " + materialModelName;
      msg += "\n**** The list of defined material models is:";
      for(materialModelIt = materialModels.begin() ; materialModelIt != materialModels.end() ; ++ materialModelIt)
        msg += "  " + materialModelIt->first + ",";
      msg += "\b\n\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
    }
    blockIt->setMaterialModel(materialModelIt->second);

    // set the damage model
    string damageModelName = blockIt->getDamageModelName();
    if(damageModelName != "None"){
      std::map< std::string, Teuchos::RCP<const PeridigmNS::DamageModel> >::iterator damageModelIt;
      damageModelIt = damageModels.find(damageModelName);
      if(damageModelIt == damageModels.end()){
        string msg = "\n**** Error, invalid damage model:  " + damageModelName;
        msg += "\n**** The list of defined damage models is:";
        for(damageModelIt = damageModels.begin() ; damageModelIt != damageModels.end() ; ++ damageModelIt)
          msg += "  " + damageModelIt->first + ",";
        msg += "\b\n\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
      }
      blockIt->setDamageModel(damageModelIt->second);
    }
  }

  // Setup contact
  initializeContact();

  // If there is a contact model, assign it to all blocks
  // \todo Refactor contact!
  if(analysisHasContact){
    Teuchos::RCP<const PeridigmNS::ContactModel> contactModel = contactModels.begin()->second;
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->setContactModel(contactModel);
  }

  // Load the auxiliary field ids into the blocks (they will be
  // combined with material model and contact model ids when allocating
  // space in the Block's DataManager)
  vector<int> auxiliaryFieldIds;

  // Force the allocation of space in all blocks for the following field data
  // \todo Replace this with a query to the requested output fields, which is why we're forcing this allocation.
  auxiliaryFieldIds.push_back(modelCoordinatesFieldId);
  auxiliaryFieldIds.push_back(coordinatesFieldId);
  auxiliaryFieldIds.push_back(displacementFieldId);
  auxiliaryFieldIds.push_back(velocityFieldId);

  // Add fields from compute classes to auxiliary field vector
  vector<int> computeManagerFieldIds = computeManager->FieldIds();
  auxiliaryFieldIds.insert(auxiliaryFieldIds.end(), computeManagerFieldIds.begin(), computeManagerFieldIds.end());

  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->setAuxiliaryFieldIds(auxiliaryFieldIds);

  // Initialize the blocks (creates maps, neighborhoods, DataManager)
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->initialize(peridigmDisc->getGlobalOwnedMap(1),
                        peridigmDisc->getGlobalOverlapMap(1),
                        peridigmDisc->getGlobalOwnedMap(3),
                        peridigmDisc->getGlobalOverlapMap(3),
                        peridigmDisc->getGlobalBondMap(),
                        blockIDs,
                        globalNeighborhoodData);

  // Create a temporary vector for storing the global element ids
  Epetra_Vector elementIds(*(peridigmDisc->getCellVolume()));
  for(int i=0 ; i<elementIds.MyLength() ; ++i)
    elementIds[i] = elementIds.Map().GID(i);

  // Load initial data into the blocks
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*(peridigmDisc->getBlockID()),    blockIdFieldId,          PeridigmField::STEP_NONE, Insert);
    blockIt->importData(*(peridigmDisc->getCellVolume()), volumeFieldId,           PeridigmField::STEP_NONE, Insert);
    blockIt->importData(*(peridigmDisc->getInitialX()),   modelCoordinatesFieldId, PeridigmField::STEP_NONE, Insert);
    blockIt->importData(*(peridigmDisc->getInitialX()),   coordinatesFieldId,      PeridigmField::STEP_N,    Insert);
    blockIt->importData(*(peridigmDisc->getInitialX()),   coordinatesFieldId,      PeridigmField::STEP_NP1,  Insert);
    blockIt->importData(elementIds,                       elementIdFieldId,        PeridigmField::STEP_NONE, Insert);
  }

  // Set the density in the mothership vector
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    Teuchos::RCP<const Epetra_BlockMap> OwnedScalarPointMap = blockIt->getOwnedScalarPointMap();
    double blockDensity = blockIt->getMaterialModel()->Density();
    for(int i=0 ; i<OwnedScalarPointMap->NumMyElements() ; ++i){
      int globalID = OwnedScalarPointMap->GID(i);
      int mothershipLocalID = oneDimensionalMap->LID(globalID);
      (*density)[mothershipLocalID] = blockDensity;
    }
  }

  // compute partial volumes
  if(analysisHasPartialVolumes){
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      computePartialVolume(Teuchos::rcpFromRef(*blockIt), peridigmDisc);
  }

  // apply initial velocities
  boundaryAndInitialConditionManager->applyInitialVelocities(x, v);

  // apply initial displacements
  boundaryAndInitialConditionManager->applyInitialDisplacements(x, u, y);

  // Initialize material models and damage models
  // Initialization functions require valid initial values, e.g. velocities and displacements.
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++) {
    blockIt->initializeMaterialModel();
    blockIt->initializeDamageModel();
  }

  // Initialize the compute classes
  computeManager->initialize(blocks);

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

  // Set default value for current time;
  timeCurrent = 0.0;
}

void PeridigmNS::Peridigm::instantiateMaterials() {
  Teuchos::ParameterList& materialParams = peridigmParams->sublist("Materials", true);
  MaterialFactory materialFactory;

  // The horizon and the finite-difference probe length will be added to the material model parameters if not already present
  double horizon = peridigmParams->sublist("Discretization", true).get<double>("Horizon");
  double finiteDifferenceProbeLength = peridigmParams->sublist("Solver", true).get<double>("Finite Difference Probe Length");
  for(Teuchos::ParameterList::ConstIterator it = materialParams.begin() ; it != materialParams.end() ; it++){
    // Get the material parameters for a given material and add the horizon if necessary
    Teuchos::ParameterList& matParams = materialParams.sublist(it->first);
    if(!matParams.isParameter("Horizon"))
       matParams.set("Horizon", horizon);
    if(!matParams.isParameter("Finite Difference Probe Length"))
      matParams.set("Finite Difference Probe Length", finiteDifferenceProbeLength);
    materialModels[it->first] = materialFactory.create(matParams) ;
  }

  TEUCHOS_TEST_FOR_EXCEPT_MSG(materialModels.size() == 0, "No material models created!");
}

void PeridigmNS::Peridigm::instantiateDamageModels() {
  Teuchos::ParameterList& damageModelParams = peridigmParams->sublist("Damage Models", true);
  DamageModelFactory damageModelFactory;
  for(Teuchos::ParameterList::ConstIterator it = damageModelParams.begin() ; it != damageModelParams.end() ; it++){
    Teuchos::ParameterList& damageParams = damageModelParams.sublist(it->first);
    damageModels[it->first] = damageModelFactory.create(damageParams);
  }
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

  oneDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*oneDimensionalMap, 3));
  blockIDs = Teuchos::rcp((*oneDimensionalMothership)(0), false);        // block ID
  volume = Teuchos::rcp((*oneDimensionalMothership)(1), false);          // cell volume
  density = Teuchos::rcp((*oneDimensionalMothership)(2), false);         // density

  // \todo Do not allocate space for the contact force nor deltaU if not needed.
  threeDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*threeDimensionalMap, 9));
  x = Teuchos::rcp((*threeDimensionalMothership)(0), false);             // initial positions
  u = Teuchos::rcp((*threeDimensionalMothership)(1), false);             // displacement
  y = Teuchos::rcp((*threeDimensionalMothership)(2), false);             // current positions
  v = Teuchos::rcp((*threeDimensionalMothership)(3), false);             // velocities
  a = Teuchos::rcp((*threeDimensionalMothership)(4), false);             // accelerations
  force = Teuchos::rcp((*threeDimensionalMothership)(5), false);         // force
  contactForce = Teuchos::rcp((*threeDimensionalMothership)(6), false);  // contact force (used only for contact simulations)
  deltaU = Teuchos::rcp((*threeDimensionalMothership)(7), false);        // increment in displacement (used only for implicit time integration)
  scratch = Teuchos::rcp((*threeDimensionalMothership)(8), false);       // scratch space

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
      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
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
    TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(name) != nodeSets->end(), "**** Duplicate node set found: " + name + "\n");
    vector<int>& nodeList = it->second;
    (*nodeSets)[name] = nodeList;
  }
}

void PeridigmNS::Peridigm::initializeContact() {

  // The horizon will be added to the contact model parameter list, if needed
  double horizon = peridigmParams->sublist("Discretization", true).get<double>("Horizon");

  // Assume no contact
  analysisHasContact = false;
  contactSearchRadius = 0.0;
  contactRebalanceFrequency = 0;

  // Set up global contact parameters for rebalance and proximity search
  if(peridigmParams->isSublist("Contact")){
    Teuchos::ParameterList & contactParams = peridigmParams->sublist("Contact");
    analysisHasContact = true;
    if(!contactParams.isParameter("Search Radius"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Radius\" not specified.");
    contactSearchRadius = contactParams.get<double>("Search Radius");
    if(!contactParams.isParameter("Search Frequency"))
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Frequency\" not specified.");
    contactRebalanceFrequency = contactParams.get<int>("Search Frequency");
  }

  // Instantiate contact models
  ContactModelFactory contactModelFactory;
  if(analysisHasContact){
    Teuchos::ParameterList& contactModelParams = peridigmParams->sublist("Contact").sublist("Models");
    for(Teuchos::ParameterList::ConstIterator it = contactModelParams.begin() ; it != contactModelParams.end() ; it++){
      // Get the parameters for a given contact model and add the horizon if necessary
      Teuchos::ParameterList& modelParams = contactModelParams.sublist(it->first);
      if(!modelParams.isParameter("Friction Coefficient"))
        modelParams.set("Friction Coefficient", 0.0);
      if(!modelParams.isParameter("Horizon"))
        modelParams.set("Horizon", horizon);
      contactModels[it->first] = contactModelFactory.create(modelParams) ;
    }
  }

  // \todo Refactor contact to allow for specific contact models between specific sets of blocks
  // Print errors/warnings to let users know contact is a work in progress
  if(analysisHasContact){
    TEUCHOS_TEST_FOR_EXCEPTION(contactModels.size() > 1, Teuchos::Exceptions::InvalidParameter, "\n**** Error, the current version of Peridigm supports only a single contact model in a given analysis.\n");
    if(peridigmParams->sublist("Contact").isSublist("Interactions")){
      if(peridigmComm->MyPID() == 0)
        cout << "\n**** Warning, the current version of Peridigm does not support the specification of contact interactions, contact will be enabled for all elements.\n" << endl;
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

void PeridigmNS::Peridigm::instantiateComputeManager() {

  Teuchos::RCP<Teuchos::ParameterList> computeParams = Teuchos::rcp( new Teuchos::ParameterList("Compute Manager") );
  Teuchos::ParameterList& outputVariables  = computeParams->sublist("Output Variables");

  // If the user has provided a "Compute Class Parameters" ParameterList, add it to computeParams
  if(peridigmParams->isSublist("Compute Class Parameters"))
    computeParams->sublist("Compute Class Parameters") = peridigmParams->sublist("Compute Class Parameters");

  // Loop over high level parameter list entries to find all output lists
  for (Teuchos::ParameterList::ConstIterator it = peridigmParams->begin(); it != peridigmParams->end(); ++it) {
    // See if name of parameterlist entry contains "Output".
    const std::string output("Output");
    const std::string name(it->first);
    size_t found = name.find(output);
    Teuchos::RCP<Teuchos::ParameterList> outputParams;
    if (found!=std::string::npos) {
      // Make copy of list
      try{
        outputParams = Teuchos::rcp( new Teuchos::ParameterList( peridigmParams->sublist(name,true) ) );
      }
      catch(const std::exception &e){
        string msg = "Peridigm::instantiateComputeManager: ";
        msg+= name;
        msg+= " is not a Teuchos::ParameterList sublist.";
        TEUCHOS_TEST_FOR_EXCEPT_MSG( true, msg );
      }
      // Create union of all requested output fields
      Teuchos::ParameterList outputVariables2 = outputParams->sublist("Output Variables");
      for(Teuchos::ParameterList::ConstIterator it = outputVariables2.begin() ; it != outputVariables2.end() ; it++){
        if (!outputVariables.isParameter(it->first)) {
          outputVariables.setEntry(it->first,it->second);
        }
      }
    }
  }
  computeManager = Teuchos::rcp( new PeridigmNS::ComputeManager( computeParams, peridigmComm  ) );
}

void PeridigmNS::Peridigm::initializeBlocks(Teuchos::RCP<AbstractDiscretization> disc) {

  // Did user specify default blocks?
  bool defaultBlocks = false;
  // Parameterlist for default blocks
  Teuchos::ParameterList defaultBlockParams;

  // Create vector of blocks
  blocks = Teuchos::rcp(new std::vector<PeridigmNS::Block>());

  // Loop over each entry in "Blocks" section of input deck. 
  Teuchos::ParameterList& blockParams = peridigmParams->sublist("Blocks", true);
  for(Teuchos::ParameterList::ConstIterator it = blockParams.begin() ; it != blockParams.end() ; it++){
    const string& name = it->first;
    Teuchos::ParameterList& params = blockParams.sublist(name);
    string blockNamesString = params.get<string>("Block Names");
    // Parse space-delimited list of block names and instantiate a Block object for each
    istringstream iss(blockNamesString);
    vector<string> blockNames;
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter<vector<string> >(blockNames));
    for(vector<string>::const_iterator it=blockNames.begin() ; it!=blockNames.end() ; ++it){
      // If "default" block encountered, record parameterlist and continue on
      if ( (*it) == "Default" || (*it) == "default" ) {
        defaultBlockParams = params;
        defaultBlocks = true;
        continue;
      }
      // Assume that the block names are "block_" + the block ID
      size_t loc = it->find_last_of('_');
      TEUCHOS_TEST_FOR_EXCEPT_MSG(loc == string::npos, "\n**** Parse error, invalid block name.\n");
      stringstream blockIDSS(it->substr(loc+1, it->size()));
      int blockID;
      blockIDSS >> blockID;
      PeridigmNS::Block block(*it, blockID, params);
      blocks->push_back(block);
    }
  }

  // Add in all default blocks
  if (defaultBlocks) {
    std::vector<std::string> discretizationBlockNames = disc->getBlockNames();
    for(vector<string>::const_iterator it=discretizationBlockNames.begin() ; it!=discretizationBlockNames.end() ; ++it){
      bool blockMatch = false;
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        // if name match, break
        if ((*it) == blockIt->getName()) {
          blockMatch = true;
          break;
        }
      }
      if (!blockMatch) { // Create new block. Assume block name are "block_" + block ID
        size_t loc = it->find_last_of('_');
        TEUCHOS_TEST_FOR_EXCEPT_MSG(loc == string::npos, "\n**** Parse error, invalid block name in discretization object.\n");
        stringstream blockIDSS(it->substr(loc+1, it->size()));
        int blockID;
        blockIDSS >> blockID;
        PeridigmNS::Block block(*it, blockID, defaultBlockParams);
        blocks->push_back(block);
      }
    }
  }

  // Ensure that there is a one-to-one match between instantiated blocks and blocks defined in the discretization object
  std::vector<std::string> discreticationBlockNames = disc->getBlockNames();
  bool blockError = false;
  if(discreticationBlockNames.size() != blocks->size())
    blockError = true;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    string blockName = blockIt->getName();
    std::vector<std::string>::iterator it = find(discreticationBlockNames.begin(), discreticationBlockNames.end(), blockName);
    if(it == discreticationBlockNames.end())
      blockError = true;
  }
  if(blockError == true){
    string msg = "\n**** Error, blocks defined in mesh do not match blocks defined in input deck.";
    msg += "\n**** List of block names in mesh:";
    for(unsigned int i=0 ; i<discreticationBlockNames.size() ; ++i)
      msg += "  " + discreticationBlockNames[i] + ",";
    msg += "\b";
    msg += "\n**** List of block names in input deck:";
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      msg += "  " + blockIt->getName()  + ",";
    msg += "\b\n\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
  }

}

void PeridigmNS::Peridigm::initializeOutputManager() {

  // Create empty container for output managers
  outputManager = Teuchos::rcp(new PeridigmNS::OutputManagerContainer() );

  Teuchos::RCP<Teuchos::ParameterList> outputParams;

  // Loop over high level parameter list entries to find all output lists
  for (Teuchos::ParameterList::ConstIterator it = peridigmParams->begin(); it != peridigmParams->end(); ++it) {
    // See if name of parameterlist entry contains "Output".
    const std::string output("Output");
    const std::string name(it->first);
    size_t found = name.find(output);
    if (found!=std::string::npos) {
      // Make copy of list
      try{
        outputParams = Teuchos::rcp( new Teuchos::ParameterList( peridigmParams->sublist(name,true) ) );
      }
      catch(const std::exception &e){
        string msg = "Peridigm::initializeOutputManager: ";
        msg+= name;
        msg+= " is not a Teuchos::ParameterList sublist.";
        TEUCHOS_TEST_FOR_EXCEPT_MSG( true, msg );
      }
      // Add proc id data to copied list
      outputParams->set("NumProc", (int)(peridigmComm->NumProc()));
      outputParams->set("MyPID", (int)(peridigmComm->MyPID()));
      // Make the default format "ExodusII"
      string outputFormat = outputParams->get("Output File Type", "ExodusII");
      TEUCHOS_TEST_FOR_EXCEPTION( outputFormat != "ExodusII",
                                  std::invalid_argument,
                                  "PeridigmNS::Peridigm: \"Output File Type\" must be \"ExodusII\".");
      if (outputFormat == "ExodusII")
        outputManager->add( Teuchos::rcp(new PeridigmNS::OutputManager_ExodusII( outputParams, this, blocks ) ) );
    }
  }

}

void PeridigmNS::Peridigm::execute() {

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);

  // allowable explicit time integration schemes:  Verlet
  if(solverParams->isSublist("Verlet"))
    executeExplicit();

  // allowable implicit time integration schemes:  Implicit, QuasiStatic
  else if(solverParams->isSublist("QuasiStatic"))    
    executeQuasiStatic();
  else if(solverParams->isSublist("NOXQuasiStatic"))    
    executeNOXQuasiStatic();
  else if(solverParams->isSublist("Implicit"))    
    executeImplicit();
}

void PeridigmNS::Peridigm::executeExplicit() {

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  Teuchos::RCP<Teuchos::ParameterList> verletParams = sublist(solverParams, "Verlet", true);

  // Compute the approximate critical time step
  double criticalTimeStep = 1.0e50;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    double blockCriticalTimeStep = ComputeCriticalTimeStep(*peridigmComm, *blockIt);
    if(blockCriticalTimeStep < criticalTimeStep)
      criticalTimeStep = blockCriticalTimeStep;
  }
  double globalCriticalTimeStep;
  peridigmComm->MinAll(&criticalTimeStep, &globalCriticalTimeStep, 1);
  double dt = globalCriticalTimeStep;
  // Query for a user-supplied time step, which overrides the computed value
  double userDefinedTimeStep = 0.0;
  if(verletParams->isParameter("Fixed dt")){
    userDefinedTimeStep = verletParams->get<double>("Fixed dt");
    dt = userDefinedTimeStep;
  }
  // Multiply the time step by the user-supplied safety factor, if provided
  double safetyFactor = 1.0;
  if(verletParams->isParameter("Safety Factor")){
    safetyFactor = verletParams->get<double>("Safety Factor");
    dt *= safetyFactor;
  }
  // Write time step information to stdout
  if(peridigmComm->MyPID() == 0){
    cout << "Time step (seconds):" << endl;
    cout << "  Stable time step    " << globalCriticalTimeStep << endl;
    if(verletParams->isParameter("Fixed dt"))
      cout << "  User time step      " << dt << endl;
    else
      cout << "  User time step      not provided" << endl;
    if(verletParams->isParameter("Safety Factor"))
      cout << "  Safety factor       " << safetyFactor << endl;
    else
      cout << "  Safety factor       not provided " << endl;
    cout << "  Time step           " << dt << "\n" << endl;
  }

  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal   = solverParams->get("Final Time", 1.0);
  timeCurrent = timeInitial;
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

  // Set the prescribed displacements (allow for nonzero initial displacements).
  // Then back compute the displacement vector.  Leave the velocity as zero.
  // \todo How do we really want to handle nonzero initial displacements?
//   boundaryAndInitialConditionManager->applyKinematicBC_SetDisplacement(timeInitial, x, u);
//   for(int i=0 ; i<u->MyLength() ; ++i)
//     (*y)[i] += (*u)[i];

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // \todo The velocity copied into the DataManager is actually the midstep velocity, not the NP1 velocity; this can be fixed by creating a midstep velocity field in the DataManager and setting the NP1 value as invalid.

  // Evaluate force in initial configuration for use in first timestep
  PeridigmNS::Timer::self().startTimer("Internal Force");
  modelEvaluator->evalModel(workset);
  PeridigmNS::Timer::self().stopTimer("Internal Force");

  // Copy force from the data manager to the mothership vector
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  force->PutScalar(0.0);
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    scratch->PutScalar(0.0);
    blockIt->exportData(*scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
    force->Update(1.0, *scratch, 1.0);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  if(analysisHasContact){
    // Copy contact force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    contactForce->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(*scratch, contactForceDensityFieldId, PeridigmField::STEP_NP1, Add);
      contactForce->Update(1.0, *scratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
    // Add contact forces to forces
    force->Update(1.0, *contactForce, 1.0);
  }

  // fill the acceleration vector
  (*a) = (*force);
  for(int i=0 ; i<a->MyLength() ; ++i)
    (*a)[i] /= (*density)[i/3];

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  this->synchDataManagers();

  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  int displayTrigger = nsteps/100;
  if(displayTrigger == 0)
    displayTrigger = 1;

  for(int step=1; step<=nsteps; step++){

    double timePrevious = timeCurrent;
    timeCurrent = timeInitial + (step*dt);

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

    // Set the velocities for dof with kinematic boundary conditions.
    // This will propagate through the Verlet integrator and result in the proper
    // displacement boundary conditions on y and consistent values for v and u.
    PeridigmNS::Timer::self().startTimer("Apply kinematic B.C.");
    boundaryAndInitialConditionManager->applyKinematicBC_SetVelocity(timeCurrent, timePrevious, x, v);
    PeridigmNS::Timer::self().stopTimer("Apply kinematic B.C.");

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
      blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");

    // Copy force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    force->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(*scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
      force->Update(1.0, *scratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");    

    // Check for NaNs in force evaluation
    // We'd like to know now because a NaN will likely cause a difficult-to-unravel crash downstream.
    for(int i=0 ; i<force->MyLength() ; ++i)
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*force)[i]), "**** NaN returned by force evaluation.\n");

    if(analysisHasContact){
      // Copy contact force from the data manager to the mothership vector
      PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      contactForce->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scratch->PutScalar(0.0);
        blockIt->exportData(*scratch, contactForceDensityFieldId, PeridigmField::STEP_NP1, Add);
        contactForce->Update(1.0, *scratch, 1.0);
      }
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Check for NaNs in contact force evaluation
      for(int i=0 ; i<contactForce->MyLength() ; ++i)
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*contactForce)[i]), "**** NaN returned by contact force evaluation.\n");

      // Add contact forces to forces
      force->Update(1.0, *contactForce, 1.0);
    }

    // fill the acceleration vector
    (*a) = (*force);
    for(int i=0 ; i<a->MyLength() ; ++i)
      (*a)[i] /= (*density)[i/3];

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt2, aptr, vptr, 1, 1);

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

bool PeridigmNS::Peridigm::computeF(const Epetra_Vector& x, Epetra_Vector& FVec, NOX::Epetra::Interface::Required::FillType fillType) {
  return evaluateNOX(fillType, &x, &FVec, 0);
}

bool PeridigmNS::Peridigm::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) {
  return evaluateNOX(NOX::Epetra::Interface::Required::Jac, &x, 0, 0);
}

bool PeridigmNS::Peridigm::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& Prec, Teuchos::ParameterList* precParams) {
  cout << "ERROR: Peridigm::preconditionVector() - Use Explicit Jacobian only for NOX interface!" << endl;
  throw "Interface Error";
}

bool PeridigmNS::Peridigm::evaluateNOX(NOX::Epetra::Interface::Required::FillType flag, 
        const Epetra_Vector* soln,
        Epetra_Vector* tmp_rhs,
        Epetra_RowMatrix* tmp_matrix)
{
  //Determine what to fill (F or Jacobian)
  bool fillF = false;
  bool fillMatrix = false;
  
  // "flag" can be used to determine how accurate your fill of F should be 
  // depending on why we are calling evaluate (Could be using computeF to 
  // populate a Jacobian or Preconditioner).
  if (flag == NOX::Epetra::Interface::Required::Residual) {
    fillF = true;
  }
  else if (flag == NOX::Epetra::Interface::Required::Jac) {
    fillMatrix = true;
  }
  else if (flag == NOX::Epetra::Interface::Required::Prec) {
    // Do nothing for now
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "flag = Prec");
  }
  else if (flag == NOX::Epetra::Interface::Required::User) {
    // Do nothing for now
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "flag = User");
  }

  // copy the solution vector passed in by NOX to update the deformation 
  for(int i=0 ; i < u->MyLength() ; ++i){
    (*y)[i] = (*x)[i] + (*u)[i] + (*soln)[i];
  }

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
  } 
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  if(fillF){
    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");

    // Copy force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    force->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(*scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
      force->Update(1.0, *scratch, 1.0);
    }
    scratch->PutScalar(0.0);
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
      
    // Create residual vector
    Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));

    // copy the internal force to the residual vector
    // note that due to restrictions on CrsMatrix, these vectors have different (but equivalent) maps
    TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != force->MyLength(), "**** PeridigmNS::Peridigm::evaluateNOX() incompatible vector lengths!\n");
    for(int i=0 ; i < force->MyLength() ; ++i)
      (*residual)[i] = (*force)[i];

    // convert force density to force
    for(int i=0 ; i < residual->MyLength() ; ++i)
      (*residual)[i] *= (*volume)[i/3];

    // zero out the rows corresponding to kinematic boundary conditions
    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);
      
    // copy back to tmp_rhs 
    for(int i=0 ; i < tmp_rhs->MyLength() ; ++i)
      (*tmp_rhs)[i] = (*residual)[i];
  }

  // Compute the tangent if requested
  if( fillMatrix && m_noxJacobianUpdateCounter%m_noxTriggerJacobianUpdate == 0 ){
    tangent->PutScalar(0.0);
    PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
    modelEvaluator->evalJacobian(workset);
    int err = tangent->GlobalAssemble();
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::evaluateNOX(), GlobalAssemble() returned nonzero error code.\n");
    PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");
    boundaryAndInitialConditionManager->applyKinematicBC_InsertZerosAndSetDiagonal(tangent);
  }
  if( fillMatrix )
    m_noxJacobianUpdateCounter += 1;

  return true;
}

void PeridigmNS::Peridigm::jacobianDiagnostics(Teuchos::RCP<NOX::Epetra::Group> noxGroup){

  stringstream ss;
  ss << "  Jacobian diagonsitics:";

  NOX::Abstract::Group::ReturnType returnValue = NOX::Abstract::Group::Ok;
  if(!noxGroup->isJacobian())
    returnValue = noxGroup->computeJacobian();
  if(returnValue != NOX::Abstract::Group::Ok){
    ss << "  failed to access jacobian";
    if(peridigmComm->MyPID() == 0)
      cout << ss.str() << endl;
    return;
  }

  // Diagnostic #1, check for symmetry

  // Construct transpose
  Teuchos::RCP<Epetra_Operator> jacobianOperator = noxGroup->getLinearSystem()->getJacobianOperator();
  Epetra_CrsMatrix* jacobian = dynamic_cast<Epetra_CrsMatrix*>(jacobianOperator.get());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(jacobian == NULL, "\n****Error: jacobianDiagnostics() failed to convert jacobian to Epetra_CrsMatrix.\n");  
  Epetra_CrsMatrix jacobianTranspose(*jacobian);
  Epetra_CrsMatrix* jacobianTransposePtr = &jacobianTranspose;
  Epetra_RowMatrixTransposer jacobianTransposer(jacobian);
  bool makeDataContiguous = false;
  int returnCode = jacobianTransposer.CreateTranspose(makeDataContiguous, jacobianTransposePtr);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(returnCode != 0, "\n****Error: jacobianDiagnostics() failed to transpose jacobian.\n");  

  // Replace entries in transpose with 0.5*(J - J^T)
  int numRows = jacobian->NumMyRows();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(numRows != jacobianTranspose.NumMyRows(), "\n****Error: jacobianDiagnostics() incompatible matrices.\n");
  int numEntries, numEntriesTranspose;
  double *values, *valuesTranspose;
  int *indices, *indicesTranspose;
  for(int i=0; i<numRows; i++){
    jacobian->ExtractMyRowView(i, numEntries, values, indices);
    jacobianTranspose.ExtractMyRowView(i, numEntriesTranspose, valuesTranspose, indicesTranspose);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numEntries != numEntriesTranspose, "\n****Error: jacobianDiagnostics() incompatible matrices.\n");
    for (int j=0; j<numEntries; j++){
      TEUCHOS_TEST_FOR_EXCEPT_MSG(indices[j] != indicesTranspose[j], "\n****Error: jacobianDiagnostics() incompatible matrices.\n");
      valuesTranspose[j] = fabs(0.5*(values[j]-valuesTranspose[j]));
    }
  }
  double normFrobenius = jacobian->NormFrobenius();
  double asymmetricNormFrobenius = jacobianTranspose.NormFrobenius();
  ss << "  asymmetric Frobenius norm = " << asymmetricNormFrobenius << " (symmetric Frobenius norm = " << normFrobenius << ")";

  // Diagnostic #2, check condition number

  if(!noxGroup->isConditionNumber()){
    int conditionNumberMaxIters = 10000;
    int conditionNumberTolerance = 1.0e3;
    int conditionNumberKrylovSubspaceSize = 100;
    int conditionNumberPrintOutput = false;
    returnValue = noxGroup->computeJacobianConditionNumber(conditionNumberMaxIters,
                                                           conditionNumberTolerance,
                                                           conditionNumberKrylovSubspaceSize,
                                                           conditionNumberPrintOutput);
  }
  if(returnValue == NOX::Abstract::Group::Ok){
    double conditionNumber = noxGroup->getJacobianConditionNumber();
    ss << ", condition number = " << conditionNumber;
  }
  else{
    ss << ", condition number calculation failed";
  }

  if(peridigmComm->MyPID() == 0)
    cout << ss.str() << endl;
}

Teuchos::RCP<Epetra_CrsMatrix> PeridigmNS::Peridigm::getJacobian() {
    return tangent;
}

void PeridigmNS::Peridigm::executeNOXQuasiStatic() {

  // Allocate memory for non-zeros in global tangent and lock in the structure
  if(peridigmComm->MyPID() == 0){
    cout << "Allocating global tangent matrix...";
    cout.flush();
  }
  PeridigmNS::Timer::self().startTimer("Allocate Global Tangent");
  allocateJacobian();
  PeridigmNS::Timer::self().stopTimer("Allocate Global Tangent");
  if(peridigmComm->MyPID() == 0){
    cout << "\n  number of rows = " << tangent->NumGlobalRows() << endl;
    cout << "  number of nonzeros = " << tangent->NumGlobalNonzeros() << "\n" << endl;
  }

  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction = Teuchos::rcp(new Epetra_Vector(force->Map()));

  // Create vectors that are specific to NOX quasi-statics.
  Teuchos::RCP<Epetra_Vector> soln = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateView));
  soln->PutScalar(0.0);

  Teuchos::RCP<Epetra_Vector> initialGuess = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<NOX::Epetra::Vector> noxInitialGuess = Teuchos::rcp(new NOX::Epetra::Vector(initialGuess, NOX::Epetra::Vector::CreateView));
  initialGuess->PutScalar(0.0);
  
  // Initialize velocity to zero
  v->PutScalar(0.0);
  
  // Create a placeholder for the timestep within the workset that is passed to the model evaluator
  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;
 
  // "Solver" parameter list
  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  bool performJacobianDiagnostics = solverParams->get("Jacobian Diagnostics", false);

  // "NOXQuasiStatic" parameter list
  Teuchos::RCP<Teuchos::ParameterList> noxQuasiStaticParams = sublist(solverParams, "NOXQuasiStatic", true);
  bool noxVerbose = noxQuasiStaticParams->get("Verbose", false);
  int maxIterations = noxQuasiStaticParams->get("Max Solver Iterations", 50);
  
  // Determine tolerance, either "Relative Tolerance" or "Absolute Tolerance"
  // If none is provided, default to relative tolerance of 1.0e-6
  double tolerance = noxQuasiStaticParams->get("Relative Tolerance", 1.0e-6);
  bool useAbsoluteTolerance = false;
  if(noxQuasiStaticParams->isParameter("Absolute Tolerance")){
    useAbsoluteTolerance = true;
    tolerance = noxQuasiStaticParams->get<double>("Absolute Tolerance");
  }

  // Parameters for printing output to the screen
  Teuchos::ParameterList& printParams = noxQuasiStaticParams->sublist("Printing");
  if (noxVerbose)
    printParams.set("Output Information",
                    NOX::Utils::OuterIteration +
                    NOX::Utils::OuterIterationStatusTest +
                    NOX::Utils::InnerIteration +
                    NOX::Utils::LinearSolverDetails +
                    NOX::Utils::Parameters +
                    NOX::Utils::Details +
                    NOX::Utils::Warning +
                    NOX::Utils::Debug +
                    NOX::Utils::TestDetails +
                    NOX::Utils::Error);
  else
    printParams.set("Output Information", NOX::Utils::Error +
                    NOX::Utils::TestDetails);
  
  // Create a print class for controlling output below 
  NOX::Utils noxPrinting(printParams);

  // Create list of time steps
  // Case 1:  User provided initial time, final time, and number of load steps
  vector<double> timeSteps;
  if( solverParams->isParameter("Final Time") && noxQuasiStaticParams->isParameter("Number of Load Steps") ){
    double timeInitial = solverParams->get("Initial Time", 0.0);
    double timeFinal = solverParams->get<double>("Final Time");
    int numLoadSteps = noxQuasiStaticParams->get<int>("Number of Load Steps");
    timeSteps.push_back(timeInitial);
    for(int i=0 ; i<numLoadSteps ; ++i)
      timeSteps.push_back(timeInitial + (i+1)*(timeFinal-timeInitial)/numLoadSteps);
  }
  // Case 2:  User provided a list of time steps
  else if( noxQuasiStaticParams->isParameter("Time Steps") ){
    string timeStepString = noxQuasiStaticParams->get<string>("Time Steps");
    istringstream iss(timeStepString);
    copy(istream_iterator<double>(iss),
	 istream_iterator<double>(),
	 back_inserter<vector<double> >(timeSteps));
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "\n****Error: No valid time step data provided.\n");
  }
  timeCurrent = timeSteps[0];

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  this->synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  // Functionality for updating the Jacobian at a user-specified interval
  // This does not appear to be available for nonlinear CG in NOX, but it's important for peridynamics, so
  // we'll handle it directly within Peridigm
  m_noxTriggerJacobianUpdate = 1;
  if(noxQuasiStaticParams->sublist("Direction").sublist("Nonlinear CG").isParameter("Update Jacobian"))
    m_noxTriggerJacobianUpdate = noxQuasiStaticParams->sublist("Direction").sublist("Nonlinear CG").get<int>("Update Jacobian");

  Epetra_Time loadStepCPUTime(*peridigmComm);
  double cummulativeLoadStepCPUTime = 0.0;

  for(int step=1 ; step<(int)timeSteps.size() ; step++){

    loadStepCPUTime.ResetStartTime();

    double timePrevious = timeCurrent;
    timeCurrent = timeSteps[step];
    double timeIncrement = timeCurrent - timePrevious;
    *timeStep = timeIncrement;
 
    m_noxJacobianUpdateCounter = 0;
    
    soln->PutScalar(0.0);
    deltaU->PutScalar(0.0);

    // Use a predictor based on the velocity from the previous load step
    for(int i=0 ; i<initialGuess->MyLength() ; ++i)
      (*initialGuess)[i] = (*v)[i]*timeIncrement;
    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(initialGuess);

    // Apply the displacement increment
    // Note that the soln vector was created to be compatible with the tangent matrix, and hence needed to be
    // constructed with an Epetra_Map.  The mothership vectors were constructed with an Epetra_BlockMap, and it's
    // this map that the boundary and intial condition manager expects.  So, make sure that the boundary and initial
    // condition manager gets the right type of vector.
    boundaryAndInitialConditionManager->applyKinematicBC_SetDisplacementIncrement(timeCurrent, timePrevious, x, deltaU);
    if(step > 1){
      // Perform line search on this guess
      for(int i=0 ; i<y->MyLength() ; ++i)
        (*y)[i] = (*x)[i] + (*u)[i] + (*deltaU)[i];
      double alpha = quasiStaticsLineSearch(residual, initialGuess);
      initialGuess->Scale(alpha);
    }

    for(int i=0 ; i<u->MyLength() ; ++i)
      (*u)[i] += (*deltaU)[i];

    *soln = *initialGuess;

    double toleranceMultiplier = 1.0;
    if(!useAbsoluteTolerance){
      // compute the vector of reactions, i.e., the forces corresponding to degrees of freedom for which kinematic B.C. are applied
      for(int i=0 ; i<y->MyLength() ; ++i)
        (*y)[i] = (*x)[i] + (*u)[i];
      computeQuasiStaticResidual(residual);
      boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(force, reaction);
      // convert force density to force
      for(int i=0 ; i<reaction->MyLength() ; ++i)
        (*reaction)[i] *= (*volume)[i/3];
      double reactionNorm2;
      reaction->Norm2(&reactionNorm2);
      toleranceMultiplier = reactionNorm2;
    }
    double residualTolerance = tolerance*toleranceMultiplier;

    // Print the load step to screen
    if(peridigmComm->MyPID() == 0)
      cout << "Load step " << step << ", initial time = " << timePrevious << ", final time = " << timeCurrent <<
        ", convergence criterion = " << residualTolerance << endl;

    // Get the linear solver parameters from the proper sublist
    // \todo Handle all allowable "Direction" settings
    Teuchos::RCP<Teuchos::ParameterList> linearSystemParams = Teuchos::rcp(new Teuchos::ParameterList);
    string directionMethod = noxQuasiStaticParams->sublist("Direction").get<string>("Method");
    if(directionMethod == "Newton")
      linearSystemParams = Teuchos::rcpFromRef( noxQuasiStaticParams->sublist("Direction").sublist("Newton").sublist("Linear Solver") );
    else if(directionMethod == "NonlinearCG")
      linearSystemParams = Teuchos::rcpFromRef( noxQuasiStaticParams->sublist("Direction").sublist("Nonlinear CG").sublist("Linear Solver") );

    // Construct the NOX linear system
    Teuchos::RCP<NOX::Epetra::Interface::Required> noxInterfaceRequired = Teuchos::RCP<NOX::Epetra::Interface::Required>(this, false);
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> noxInterfaceJacobian = Teuchos::RCP<NOX::Epetra::Interface::Jacobian>(this, false);
    Teuchos::RCP<Epetra_RowMatrix> noxJacobian = getJacobian();
    const NOX::Epetra::Vector& noxCloneVector = *soln;
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linSys = 
      Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                        *linearSystemParams,
                                                        noxInterfaceRequired, 
                                                        noxInterfaceJacobian, 
                                                        noxJacobian,
                                                        noxCloneVector));

    // Create the Group
    NOX::Epetra::Vector noxInitialGuess(soln, NOX::Epetra::Vector::CreateView);
    Teuchos::RCP<NOX::Epetra::Group> noxGroup = Teuchos::rcp(new NOX::Epetra::Group(printParams, noxInterfaceRequired, noxInitialGuess, linSys));

    // Create the convergence tests
    //NOX::Abstract::Vector::NormType normType = NOX::Abstract::Vector::NormType::TwoNorm; // OneNorm, TwoNorm, MaxNorm
    //NOX::StatusTest::ToleranceType toleranceType = NOX::StatusTest::Absolute;
    NOX::StatusTest::NormF::ScaleType scaleType = NOX::StatusTest::NormF::Unscaled;
    // The following constructor defaults NormType to TwoNorm and ToleranceType to Absolute
    Teuchos::RCP<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(residualTolerance, scaleType));
    Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(maxIterations));
    Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);
    Teuchos::RCP<NOX::StatusTest::Combo> combo = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    combo->addStatusTest(fv);
    combo->addStatusTest(absresid);
    combo->addStatusTest(maxiters);

    // Create the solver and solve
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(noxGroup, combo, noxQuasiStaticParams);
    NOX::StatusTest::StatusType noxSolverStatus = NOX::StatusTest::Unevaluated;
    int solverIteration = 1;
    while(noxSolverStatus != NOX::StatusTest::Converged && noxSolverStatus != NOX::StatusTest::Failed){
      
      // carry out nonlinear iteration
      noxSolverStatus = solver->step();

      // print results to screen
      if(performJacobianDiagnostics)
        jacobianDiagnostics(noxGroup);
      double errorNorm = solver->getSolutionGroupPtr()->getFPtr()->norm();
      if(peridigmComm->MyPID() == 0)
        cout << "  iteration " << solverIteration << ": residual L2 norm = " << errorNorm << endl;

      solverIteration += 1;
    }
    TEUCHOS_TEST_FOR_EXCEPT_MSG(noxSolverStatus != NOX::StatusTest::Converged, "\n****Error:  NOX solver failed to solve system.\n");

    // Get the Epetra_Vector with the final solution from the solver
    const Epetra_Vector& finalSolution = 
      dynamic_cast<const NOX::Epetra::Vector&>(solver->getSolutionGroup().getX()).getEpetraVector();
 
    // Output the parameter list
    //if (printing.isPrintType(NOX::Utils::Parameters)) {
//     printing.out() << endl << "Final NOX Parameters" << endl << "****************" << endl;
//     solver->getList().print(printing.out());
//     printing.out() << endl;
    //}

    // Print load step timing information
    double CPUTime = loadStepCPUTime.ElapsedTime();
    cummulativeLoadStepCPUTime += CPUTime;
    if(peridigmComm->MyPID() == 0)
      cout << setprecision(2) << "  cpu time for load step = " << CPUTime << " sec., cummulative cpu time = " << cummulativeLoadStepCPUTime << " sec.\n" << endl;

    // Store the velocity
    v->PutScalar(0.0);
    boundaryAndInitialConditionManager->applyKinematicBC_SetDisplacementIncrement(timeCurrent, timePrevious, x, v);
    for(int i=0 ; i<v->MyLength() ; ++i){
      (*v)[i] += finalSolution[i];
      (*v)[i] /= timeIncrement;
    }

    // Add the converged displacement increment to the displacement
    for(int i=0 ; i<u->MyLength() ; ++i)
      (*u)[i] += finalSolution[i];

    // Write output for completed load step
    PeridigmNS::Timer::self().startTimer("Output");
    this->synchDataManagers();
    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->updateState();
  }
  if(peridigmComm->MyPID() == 0)
    cout << endl;
}

void PeridigmNS::Peridigm::executeQuasiStatic() {

  // Allocate memory for non-zeros in global tangent and lock in the structure
  if(peridigmComm->MyPID() == 0){
    cout << "Allocating global tangent matrix...";
    cout.flush();
  }
  PeridigmNS::Timer::self().startTimer("Allocate Global Tangent");
  allocateJacobian();
  PeridigmNS::Timer::self().stopTimer("Allocate Global Tangent");
  if(peridigmComm->MyPID() == 0){
    cout << "\n  number of rows = " << tangent->NumGlobalRows() << endl;
    cout << "  number of nonzeros = " << tangent->NumGlobalNonzeros() << "\n" << endl;
  }

  // Create vectors that are specific to quasi-statics.
  // These must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction = Teuchos::rcp(new Epetra_Vector(force->Map()));

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  bool solverVerbose = solverParams->get("Verbose", false);
  Teuchos::RCP<Teuchos::ParameterList> quasiStaticParams = sublist(solverParams, "QuasiStatic", true);
  int maxSolverIterations = quasiStaticParams->get("Maximum Solver Iterations", 10);
  double dampedNewtonDiagonalScaleFactor = quasiStaticParams->get("Damped Newton Diagonal Scale Factor", 1.0001);
  double dampedNewtonDiagonalShiftFactor = quasiStaticParams->get("Damped Newton Diagonal Shift Factor", 0.00001);

  // Determine tolerance
  double tolerance = quasiStaticParams->get("Relative Tolerance", 1.0e-6);
  bool useAbsoluteTolerance = false;
  if(quasiStaticParams->isParameter("Absolute Tolerance")){
    useAbsoluteTolerance = true;
    tolerance = quasiStaticParams->get<double>("Absolute Tolerance");
  }

  // Pointer index into sub-vectors for use with BLAS
  double *xptr, *uptr, *yptr, *vptr, *aptr;
  x->ExtractView( &xptr );
  u->ExtractView( &uptr );
  y->ExtractView( &yptr );
  v->ExtractView( &vptr );
  a->ExtractView( &aptr );

  // Initialize velocity to zero
  v->PutScalar(0.0);

  // Data for Belos linear solver object
  Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> linearProblem;
  string linearSolver =  quasiStaticParams->get("Belos Linear Solver", "BlockCG");
  Belos::ReturnType isConverged;
  Teuchos::ParameterList belosList;
  belosList.set( "Block Size", 1 );  // Use single-vector iteration
  belosList.set( "Maximum Iterations", quasiStaticParams->get("Belos Maximum Iterations", tangent->NumGlobalRows()) ); // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", quasiStaticParams->get("Belos Relative Tolerance", 1.0e-4) ); // Relative convergence tolerance requested
  belosList.set( "Output Frequency", -1 );
  int verbosity = Belos::Errors + Belos::Warnings;
  if( quasiStaticParams->get("Belos Print Status", false) == true ){
    verbosity += Belos::StatusTestDetails;
    belosList.set( "Output Frequency", 1 );
  }
  belosList.set( "Verbosity", verbosity );
  belosList.set( "Output Style", Belos::Brief );
  Teuchos::RCP< Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> > belosSolver;
  if (linearSolver == "BlockGMRES") {
    belosList.set( "Num Blocks", 500); // Maximum number of blocks in Krylov factorization
    belosList.set( "Maximum Restarts", 10 ); // Maximum number of restarts allowed
    belosSolver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double,Epetra_MultiVector,Epetra_Operator>(Teuchos::rcp(&linearProblem,false), Teuchos::rcp(&belosList,false)) );
  }
  else {
    linearProblem.setHermitian(); // Assume matrix is Hermitian
    belosSolver = Teuchos::rcp( new Belos::BlockCGSolMgr<double,Epetra_MultiVector,Epetra_Operator>(Teuchos::rcp(&linearProblem,false), Teuchos::rcp(&belosList,false)) );
  }

  // Create list of time steps
  
  // Case 1:  User provided initial time, final time, and number of load steps
  vector<double> timeSteps;
  if( solverParams->isParameter("Final Time") && quasiStaticParams->isParameter("Number of Load Steps") ){
    double timeInitial = solverParams->get("Initial Time", 0.0);
    double timeFinal = solverParams->get<double>("Final Time");
    int numLoadSteps = quasiStaticParams->get<int>("Number of Load Steps");
    timeSteps.push_back(timeInitial);
    for(int i=0 ; i<numLoadSteps ; ++i)
      timeSteps.push_back(timeInitial + (i+1)*(timeFinal-timeInitial)/numLoadSteps);
  }
  // Case 2:  User provided a list of time steps
  else if( quasiStaticParams->isParameter("Time Steps") ){
    string timeStepString = quasiStaticParams->get<string>("Time Steps");
    istringstream iss(timeStepString);
    copy(istream_iterator<double>(iss),
	 istream_iterator<double>(),
	 back_inserter<vector<double> >(timeSteps));
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "\n****Error: No valid time step data provided.\n");
  }

  timeCurrent = timeSteps[0];

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  this->synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  Epetra_Time loadStepCPUTime(*peridigmComm);
  double cummulativeLoadStepCPUTime = 0.0;

  for(int step=1 ; step<(int)timeSteps.size() ; step++){

    loadStepCPUTime.ResetStartTime();

    double timePrevious = timeCurrent;
    timeCurrent = timeSteps[step];
    double timeIncrement = timeCurrent - timePrevious;
    *timeStep = timeIncrement;

    // Update nodal positions for nodes with kinematic B.C.
    deltaU->PutScalar(0.0);
    boundaryAndInitialConditionManager->applyKinematicBC_SetDisplacementIncrement(timeCurrent, timePrevious, x, deltaU);

    // Set the current position
    // \todo We probably want to rework this so that the material models get valid x, u, and y values.
    // Currently the u values are from the previous load step (and if we update u here we'll be unable
    // to properly undo a time step, which we'll need to adaptive time stepping).
    for(int i=0 ; i<y->MyLength() ; ++i)
      (*y)[i] = (*x)[i] + (*u)[i] + (*deltaU)[i];

    // compute the residual
    double residualNorm = computeQuasiStaticResidual(residual);

    double toleranceMultiplier = 1.0;
    if(!useAbsoluteTolerance){
      // compute the vector of reactions, i.e., the forces corresponding to degrees of freedom for which kinematic B.C. are applied
      boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(force, reaction);
      // convert force density to force
      for(int i=0 ; i<reaction->MyLength() ; ++i)
        (*reaction)[i] *= (*volume)[i/3];
      double reactionNorm2;
      reaction->Norm2(&reactionNorm2);
      toleranceMultiplier = reactionNorm2;
    }

    if(peridigmComm->MyPID() == 0)
      cout << "Load step " << step << ", initial time = " << timePrevious << ", final time = " << timeCurrent <<
        ", convergence criterion = " << tolerance*toleranceMultiplier << endl;

    int solverIteration = 1;
    bool dampedNewton = false;
    bool usePreconditioner = true;
    int numPureNewtonSteps = 8;
    int numPreconditionerSteps = 24;
    int dampedNewtonNumStepsBetweenTangentUpdates = 8;
    double alpha = 0.0;
    while(residualNorm > tolerance*toleranceMultiplier && solverIteration <= maxSolverIterations){

      if(!solverVerbose){
        if(peridigmComm->MyPID() == 0)
          cout << "  iteration " << solverIteration << ": residual = " << residualNorm << endl;
      }
      else{
        double residualL2, residualInf;
        residual->Norm2(&residualL2);
        residual->NormInf(&residualInf);
        if(peridigmComm->MyPID() == 0)
          cout << "  iteration " << solverIteration << ": residual = " << residualNorm << ", residual L2 = " << residualL2 << ", residual inf = " << residualInf << ", alpha = " << alpha << endl;
      }

      // On the first iteration, use a predictor based on the velocity from the previous load step
      if(solverIteration == 1 && step > 1){
        for(int i=0 ; i<lhs->MyLength() ; ++i)
          (*lhs)[i] = (*v)[i]*timeIncrement;
        boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(lhs);
        isConverged = Belos::Converged;
      }
      else{
        // If we reach the specified maximum number of iterations of the nonlinear solver, switch to a damped Newton approach.
        // This should hurt the convergence rate but improve robustness.
        if(solverIteration > numPureNewtonSteps && !dampedNewton){
          if(peridigmComm->MyPID() == 0)
            cout << "  --switching nonlinear solver to damped Newton--" << endl;
          dampedNewton = true;
        }

        // If we reach the specified maximum number of iterations of the nonlinear solver, disable the preconditioner
        if(solverIteration > numPreconditionerSteps && usePreconditioner){
          if(peridigmComm->MyPID() == 0)
            cout << "  --disabling preconditioner--" << endl;
          usePreconditioner = false;
        }

        // Compute the tangent
        if( !dampedNewton || (solverIteration-numPureNewtonSteps-1)%dampedNewtonNumStepsBetweenTangentUpdates==0 ){
          tangent->PutScalar(0.0);
          // Set the coordinates at which the tangent will be evaluated
          if( solverIteration<7 || solverIteration%10==0 ){
            PeridigmNS::Timer::self().startTimer("Gather/Scatter");
            for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
              blockIt->importData(*y, tangentReferenceCoordinatesFieldId, PeridigmField::STEP_NONE, Insert);
            PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
          }
          else{
            if(peridigmComm->MyPID() == 0)
              cout << "DEBUGGING MESSAGE:  RETAINING TANGENT REFERENCE COORDINATES" << endl;
          }
          PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
          modelEvaluator->evalJacobian(workset);
          int err = tangent->GlobalAssemble();
          TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::executeQuasiStatic(), GlobalAssemble() returned nonzero error code.\n");
          PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");
          boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);
          boundaryAndInitialConditionManager->applyKinematicBC_InsertZerosAndSetDiagonal(tangent);
          tangent->Scale(-1.0);
          if(dampedNewton)
            quasiStaticsDampTangent(dampedNewtonDiagonalScaleFactor, dampedNewtonDiagonalShiftFactor);
          if(usePreconditioner)
            quasiStaticsSetPreconditioner(linearProblem);
	}
	
        // Solve linear system
        isConverged = quasiStaticsSolveSystem(residual, lhs, linearProblem, belosSolver);

        if(isConverged == Belos::Unconverged){
          // Adjust the tangent and try again
          if(peridigmComm->MyPID() == 0)
            cout << "  --switching nonlinear solver to damped Newton and deactivating preconditioner--" << endl;
          if(!dampedNewton)
            quasiStaticsDampTangent(dampedNewtonDiagonalScaleFactor, dampedNewtonDiagonalShiftFactor);
          dampedNewton = true;
          linearProblem.setLeftPrec( Teuchos::RCP<Belos::EpetraPrecOp>() );
          usePreconditioner = false;
          isConverged = quasiStaticsSolveSystem(residual, lhs, linearProblem, belosSolver);
        }
      }
      
      if(isConverged == Belos::Converged){

        // Zero out the entries corresponding to the kinematic boundary conditions
        // The solver should have returned zeros, but there may be small errors.
        boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(lhs);

	PeridigmNS::Timer::self().startTimer("Line Search");
        alpha = quasiStaticsLineSearch(residual, lhs);
	PeridigmNS::Timer::self().stopTimer("Line Search");

        // Apply increment to nodal positions
        for(int i=0 ; i<y->MyLength() ; ++i)
          (*deltaU)[i] += alpha*(*lhs)[i];
        for(int i=0 ; i<y->MyLength() ; ++i)
          (*y)[i] = (*x)[i] + (*u)[i] + (*deltaU)[i];
      
        // Compute residual
        residualNorm = computeQuasiStaticResidual(residual);

        solverIteration++;
      }
      else{
        if(peridigmComm->MyPID() == 0)
          cout << "\nError:  Belos linear solver failed to converge." << endl;
        residualNorm = 1.0e50;
        break;
      }
    } // end loop of nonlinear iterations

    if(solverIteration >= maxSolverIterations && peridigmComm->MyPID() == 0)
      cout << "\nWarning:  Nonlinear solver failed to converge in maximum allowable iterations." << endl;

    // If the maximum allowable number of load step reductions has been reached and the residual
    // is within a reasonable tolerance, then just accept the solution and forge ahead.
    // If not, abort the analysis.
    if(residualNorm > tolerance*toleranceMultiplier){
      if(residualNorm < 100.0*tolerance*toleranceMultiplier){
	if(peridigmComm->MyPID() == 0)
	  cout << "\nWarning:  Accepting current solution and progressing to next load step.\n" << endl;
      }
      else{
	if(peridigmComm->MyPID() == 0)
	  cout << "\nError:  Aborting analysis.\n" << endl;
	break;
      }
    }

    // The load step is complete
    // Update internal data and move on to the next load step

    if(!solverVerbose){
      if(peridigmComm->MyPID() == 0)
	cout << "  iteration " << solverIteration << ": residual = " << residualNorm << endl;
    }
    else{
      double residualL2, residualInf;
      residual->Norm2(&residualL2);
      residual->NormInf(&residualInf);
      if(peridigmComm->MyPID() == 0)
	cout << "  iteration " << solverIteration << ": residual = " << residualNorm << ", residual L2 = " << residualL2 << ", residual inf = " << residualInf << ", alpha = " << alpha << endl;
    }

    // Print load step timing information
    double CPUTime = loadStepCPUTime.ElapsedTime();
    cummulativeLoadStepCPUTime += CPUTime;
    if(peridigmComm->MyPID() == 0)
      cout << setprecision(2) << "  cpu time for load step = " << CPUTime << " sec., cummulative cpu time = " << cummulativeLoadStepCPUTime << " sec.\n" << endl;

    // Store the velocity
    for(int i=0 ; i<v->MyLength() ; ++i)
      (*v)[i] = (*deltaU)[i]/timeIncrement;

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

  } // end loop over load steps

  if(peridigmComm->MyPID() == 0)
    cout << endl;
}

void PeridigmNS::Peridigm::quasiStaticsSetPreconditioner(Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>& linearProblem) {
  Ifpack IFPFactory;
  Teuchos::ParameterList ifpackList;

  if (linearProblem.isHermitian()) { // assume matrix Hermitian; construct IC preconditioner
    Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( IFPFactory.Create("IC", &(*tangent), 0) );
    Teuchos::ParameterList ifpackList;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(Prec->SetParameters(ifpackList), 
  		      "**** PeridigmNS::Peridigm::executeQuasiStatic(), Prec->SetParameters() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(Prec->Initialize(), 
		      "**** PeridigmNS::Peridigm::executeQuasiStatic(), Prec->Initialize() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(Prec->Compute(), 
		      "**** PeridigmNS::Peridigm::executeQuasiStatic(), Prec->Compute() returned nonzero error code.\n");
    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( Prec ) );
    linearProblem.setLeftPrec( belosPrec );
  }
  else { // assume matrix non-Hermitian; construct ILU preconditioner
    std::string PrecType = "ILU"; // incomplete LU
    int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1, param is ignored.
    Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( IFPFactory.Create(PrecType, &(*tangent), OverlapLevel) );
    // specify parameters for ILU
    ifpackList.set("fact: drop tolerance", 1e-9);
    ifpackList.set("fact: ilut level-of-fill", 1);
    // the combine mode is on the following: "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
    ifpackList.set("schwarz: combine mode", "Add");
    // sets the parameters
    TEUCHOS_TEST_FOR_EXCEPT_MSG(Prec->SetParameters(ifpackList), 
  		      "**** PeridigmNS::Peridigm::executeQuasiStatic(), Prec->SetParameters() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(Prec->Initialize(), 
		      "**** PeridigmNS::Peridigm::executeQuasiStatic(), Prec->Initialize() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(Prec->Compute(), 
		      "**** PeridigmNS::Peridigm::executeQuasiStatic(), Prec->Compute() returned nonzero error code.\n");
    // Create the Belos preconditioned operator from the Ifpack preconditioner.
    // NOTE:  This is necessary because Belos expects an operator to apply the
    //        preconditioner with Apply() NOT ApplyInverse().
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp( new Belos::EpetraPrecOp( Prec ) );
    linearProblem.setLeftPrec( belosPrec );
  }
}

void PeridigmNS::Peridigm::quasiStaticsDampTangent(double dampedNewtonDiagonalScaleFactor,
                                                   double dampedNewtonDiagonalShiftFactor) {
  // Create a vector to store the diagonal
  static Teuchos::RCP<Epetra_Vector> diagonal;
  if(diagonal.is_null() || !diagonal->Map().SameAs(tangent->Map()))
    diagonal = Teuchos::rcp(new Epetra_Vector(tangent->Map()));

  // Extract the diagonal, modify it, and re-insert it into the tangent
  tangent->ExtractDiagonalCopy(*diagonal);
  diagonal->Scale(dampedNewtonDiagonalScaleFactor);
  double diagonalNormInf;
  diagonal->NormInf(&diagonalNormInf);
  double* diagonalPtr;
  diagonal->ExtractView(&diagonalPtr);
  for(int i=0 ; i<diagonal->MyLength() ; ++i)
    diagonalPtr[i] += dampedNewtonDiagonalShiftFactor*diagonalNormInf;
  tangent->ReplaceDiagonalValues(*diagonal);
}

Belos::ReturnType PeridigmNS::Peridigm::quasiStaticsSolveSystem(Teuchos::RCP<Epetra_Vector> residual,
								Teuchos::RCP<Epetra_Vector> lhs,
								Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>& linearProblem,
								Teuchos::RCP< Belos::SolverManager<double,Epetra_MultiVector,Epetra_Operator> >& belosSolver)
{
  PeridigmNS::Timer::self().startTimer("Solve Linear System");

  Belos::ReturnType isConverged(Belos::Unconverged);

  lhs->PutScalar(0.0);
  linearProblem.setOperator(tangent);
  bool isSet = linearProblem.setProblem(lhs, residual);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!isSet, "**** Belos::LinearProblem::setProblem() returned nonzero error code.\n");
  try{
    isConverged = belosSolver->solve();
  }
  catch(const std::exception &e){
    if(peridigmComm->MyPID() == 0)
      cout << "Warning:  Belos linear solver aborted with " << Teuchos::typeName(e) << "." << endl; // can dump details by printing e.what()
    isConverged = Belos::Unconverged;
  }

  PeridigmNS::Timer::self().stopTimer("Solve Linear System");

  // Debugging code: Debug linear system to disk
  bool writeMatrixNow = false;
  static int solverCount = 1;
  solverCount++;
  if (writeMatrixNow) {
     char matFilename[50];
     char LHSFilename[50];
     char RHSFilename[50];
     sprintf(matFilename,"A_%03i.mat",solverCount);
     sprintf(LHSFilename,"x_%03i.mat",solverCount);
     sprintf(RHSFilename,"b_%03i.mat",solverCount);
     EpetraExt::RowMatrixToMatrixMarketFile(matFilename, *tangent, "Matrix", "Matrix");
     EpetraExt::MultiVectorToMatrixMarketFile(LHSFilename, *(linearProblem.getLHS()), "LHS", "LHS", true );
     EpetraExt::MultiVectorToMatrixMarketFile(RHSFilename, *(linearProblem.getRHS()), "RHS", "RHS", true );
  }

  return isConverged;
}

double PeridigmNS::Peridigm::quasiStaticsLineSearch(Teuchos::RCP<Epetra_Vector> residual,
						    Teuchos::RCP<Epetra_Vector> lhs)
{
  Teuchos::RCP<Epetra_Vector> tempVector = Teuchos::rcp(new Epetra_Vector(*deltaU));

  double *lhsPtr, *residualPtr, *deltaUPtr, *xPtr, *yPtr, *uPtr;
  lhs->ExtractView(&lhsPtr);
  residual->ExtractView(&residualPtr);
  deltaU->ExtractView(&deltaUPtr);
  x->ExtractView(&xPtr);
  y->ExtractView(&yPtr);
  u->ExtractView(&uPtr);

  // compute the current residual
  double unperturbedResidualNorm = computeQuasiStaticResidual(residual);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite(unperturbedResidualNorm), "**** NaN detected in residual calculation in quasiStaticsLineSearch().\n");
  if(unperturbedResidualNorm == 0.0)
    return 0.0;

  double bestAlpha = 1.0;
  double bestResidual = 1.0e50;

  vector<double> candidateAlphas;

  // a systematic guess for alpha
  double epsilon = 1.0e-4;
  for(int i=0 ; i<y->MyLength() ; ++i)
    deltaUPtr[i] += epsilon*lhsPtr[i];
  for(int i=0 ; i<y->MyLength() ; ++i)
    yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
  Teuchos::RCP<Epetra_Vector> perturbedResidual = Teuchos::rcp(new Epetra_Vector(*residual));
  computeQuasiStaticResidual(perturbedResidual);
  double SR, SPerturbedR;
  lhs->Dot(*residual, &SR);
  lhs->Dot(*perturbedResidual, &SPerturbedR);
  double tempAlpha = -1.0*epsilon*SR/(SPerturbedR - SR);
  if(tempAlpha > -0.1 && tempAlpha < 10.0)
    candidateAlphas.push_back(tempAlpha);

  // include a few brute-force gueses
  candidateAlphas.push_back(0.1);
  candidateAlphas.push_back(0.2);
  candidateAlphas.push_back(0.3);
  candidateAlphas.push_back(0.5);
  candidateAlphas.push_back(0.75);
  candidateAlphas.push_back(1.0);
  
  // compute the residual for each candidate alpha
  for(unsigned int i=0 ; i<candidateAlphas.size(); ++i){
    double alpha = candidateAlphas[i];
    for(int i=0 ; i<y->MyLength() ; ++i)
      deltaUPtr[i] += alpha*lhsPtr[i];
    for(int i=0 ; i<y->MyLength() ; ++i)
      yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
    double residualNorm = computeQuasiStaticResidual(residual);
    *deltaU = *tempVector;
    if(boost::math::isfinite(residualNorm) && residualNorm < bestResidual){
      bestAlpha = alpha;
      bestResidual = residualNorm;
    }
  }

  // if the residual is reduced by 2%, then call it quits
  double percentReduced = (unperturbedResidualNorm - bestResidual)/unperturbedResidualNorm;
  if(percentReduced > 0.02)
    return bestAlpha;

  // if the residual has not been reduced by 2%, try additional brute-force searches

  // numerous brute-force guesses for alpha
  candidateAlphas.clear();
  for(int i=0 ; i<25 ; ++i)
    candidateAlphas.push_back(0.002*(i+1));
  for(int i=0 ; i<15 ; ++i)
    candidateAlphas.push_back(0.01*(i+1)+0.05);
  for(int i=0 ; i<6 ; ++i)
    candidateAlphas.push_back(0.05*(i+1)+0.2);
  candidateAlphas.push_back(0.75);
  candidateAlphas.push_back(1.0);
  candidateAlphas.push_back(1.5);
  candidateAlphas.push_back(2.0);
  candidateAlphas.push_back(5.0);
  candidateAlphas.push_back(10.0);

  // for alphas that increase the residual, track the percent increases,
  // if no alpha can reduce the residual one of these alphas will be
  // chosen in attempt to bump the solver out of a local minimum

  vector< pair<double,double> > residualData;

  // compute the residual for each candidate alpha
  for(unsigned int i=0 ; i<candidateAlphas.size(); ++i){
    double alpha = candidateAlphas[i];
    for(int i=0 ; i<y->MyLength() ; ++i)
      deltaUPtr[i] += alpha*lhsPtr[i];
    for(int i=0 ; i<y->MyLength() ; ++i)
      yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
    double residualNorm = computeQuasiStaticResidual(residual);
    *deltaU = *tempVector;
    if(boost::math::isfinite(residualNorm)){
      if(residualNorm < bestResidual){
	bestAlpha = alpha;
	bestResidual = residualNorm;
      }
      else if(residualNorm > unperturbedResidualNorm){
	double percentIncreased = (residualNorm - unperturbedResidualNorm)/unperturbedResidualNorm;
	residualData.push_back( pair<double,double>(alpha, percentIncreased) );
      }
    }
  }

  // if the residual can be reduced by 0.01%, call it quits
  percentReduced = (unperturbedResidualNorm - bestResidual)/unperturbedResidualNorm;
  if(percentReduced > 0.0001)
    return bestAlpha;  

  // if the residual cannot be effectively reduced, try to bounce the solver out of a local
  // minimum by increasing the residual by 10%
  double targetIncrease = 0.1;
  if(peridigmComm->MyPID() == 0)
    cout << "  --line search unable to reduce residual, attempting to increase residual by 10% in hopes of escaping local minimum--" << endl;
  double bestResult = 1.0e50;
  double bestIncrease = 1.0e50;
  for(unsigned int i=0 ; i<residualData.size() ; ++i){
    double alpha = residualData[i].first;
    double percentIncrease = residualData[i].second;
    if( fabs(percentIncrease - targetIncrease) < bestResult ){
      bestAlpha = alpha;
      bestIncrease = percentIncrease;
      bestResult = fabs(percentIncrease - targetIncrease);
    }
  }
  if(peridigmComm->MyPID() == 0)
    cout << "  --selecting alpha = " << bestAlpha << ", increasing residual by " << bestIncrease << "--" << endl;

  return bestAlpha;
}

void PeridigmNS::Peridigm::executeImplicit() {

  // Allocate memory for non-zeros in global Jacobain and lock in the structure
  if(peridigmComm->MyPID() == 0)
    cout << "Allocating global tangent matrix...";
  PeridigmNS::Timer::self().startTimer("Allocate Global Tangent");
  allocateJacobian();
  PeridigmNS::Timer::self().stopTimer("Allocate Global Tangent");
  if(peridigmComm->MyPID() == 0){
    cout << "\n  number of rows = " << tangent->NumGlobalRows() << endl;
    cout << "  number of nonzeros = " << tangent->NumGlobalNonzeros() << "\n" << endl;
  }

  // Create vectors that are specific to implicit dynamics
  // The residual must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> un = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> vn = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> an = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver", true);
  Teuchos::RCP<Teuchos::ParameterList> quasiStaticParams = sublist(solverParams, "Implicit", true);
  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal = solverParams->get("Final Time", 1.0);
  timeCurrent = timeInitial;
  double absoluteTolerance       = quasiStaticParams->get("Absolute Tolerance", 1.0e-6);
  int maxSolverIterations        = quasiStaticParams->get("Maximum Solver Iterations", 10);
  double dt                      = quasiStaticParams->get("Fixed dt", 1.0);
  double beta                    = quasiStaticParams->get("Beta", 0.25);
  double gamma                   = quasiStaticParams->get("Gamma", 0.50);
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

  // Create temporary owned (e.g., "mothership") vectors for data at timestep n
  // to be used in Newmark integration
  Teuchos::RCP<Epetra_Vector> u2 = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> v2 = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

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
      blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*a, accelerationFieldId, PeridigmField::STEP_NP1, Insert);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");

    // Copy force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->exportData(*force, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Compute the residual
    // residual = beta*dt*dt*(M*a - force)
    // Note that due to restrictions to CrsMatrix, the residual has a different (but equivalent) map
    // than the force and acceleration
    for(int i=0 ; i<residual->MyLength() ; ++i)
      (*residual)[i] = beta*dt2*( (*density)[i/3] * (*a)[i] - (*force)[i] );

    // Modify residual for kinematic BC
    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);

    double residualNorm;
    residual->Norm2(&residualNorm);

    int NLSolverIteration = 0;
    while(residualNorm > absoluteTolerance && NLSolverIteration <= maxSolverIterations){

      if(peridigmComm->MyPID() == 0)
        cout << "Time step " << step << ", Newton iteration = " << NLSolverIteration << ", norm(residual) = " << residualNorm << endl;

      // Fill the Jacobian
      computeImplicitJacobian(beta);

      // Modify Jacobian for kinematic BC
      boundaryAndInitialConditionManager->applyKinematicBC_InsertZerosAndSetDiagonal(tangent);

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
      Belos::ReturnType isConverged = belosSolver->solve();
      if(isConverged != Belos::Converged && peridigmComm->MyPID() == 0)
        cout << "Warning:  Belos linear solver failed to converge!  Proceeding with nonconverged solution..." << endl;
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
        blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(*a, accelerationFieldId, PeridigmField::STEP_NP1, Insert);
      }
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Update forces based on new positions
      PeridigmNS::Timer::self().startTimer("Internal Force");
      modelEvaluator->evalModel(workset);
      PeridigmNS::Timer::self().stopTimer("Internal Force");

      // Copy force from the data manager to the mothership vector
      PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->exportData(*force, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Compute residual vector and its norm
      // residual = beta*dt*dt*(M*a - force)
      // Note that due to restrictions to CrsMatrix, the residual has a different (but equivalent) map
      // than the force and acceleration
      for(int i=0 ; i<residual->MyLength() ; ++i)
        (*residual)[i] = beta*dt2*( (*density)[i/3] * (*a)[i] - (*force)[i] );

      // Modify residual for kinematic BC
      boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);

      residual->Norm2(&residualNorm);

      NLSolverIteration++;
    }

    if(peridigmComm->MyPID() == 0)
      cout << "Time step " << step << ", Newton iteration = " << NLSolverIteration << ", norm(residual) = " << residualNorm << endl;

    timeCurrent = timeInitial + (step+1)*dt;

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
  vector<int> myGlobalElements(numMyElements);
  int* oneDimensionalMapGlobalElements = oneDimensionalMap->MyGlobalElements();
  for(int iElem=0 ; iElem<oneDimensionalMap->NumMyElements() ; ++iElem){
    myGlobalElements[3*iElem]     = 3*oneDimensionalMapGlobalElements[iElem];
    myGlobalElements[3*iElem + 1] = 3*oneDimensionalMapGlobalElements[iElem] + 1;
    myGlobalElements[3*iElem + 2] = 3*oneDimensionalMapGlobalElements[iElem] + 2;
  }
  int indexBase = 0;
  tangentMap = Teuchos::rcp(new Epetra_Map(numGlobalElements, numMyElements, &myGlobalElements[0], indexBase, *peridigmComm));
  myGlobalElements.clear();

  // Create the global tangent matrix
  Epetra_DataAccess CV = Copy;
  int numEntriesPerRow = 0;  // Indicates allocation will take place during the insertion phase
  bool ignoreNonLocalEntries = false;
  tangent = Teuchos::rcp(new Epetra_FECrsMatrix(CV, *tangentMap, numEntriesPerRow, ignoreNonLocalEntries));

  // Store nonzero columns for each row, with everything in global indices
  map<int, set<int> > rowEntries;

  // Loop over the neighborhood for each locally-owned point and record non-zero entries in the matrix.
  // Entries will exist for any two points that are bonded, and any two points that are bonded to a common third point.
  int* neighborhoodList = globalNeighborhoodData->NeighborhoodList();
  int neighborhoodListIndex = 0;
  vector<int> globalIndices;
  int numOwnedPoints = globalNeighborhoodData->NumOwnedPoints();
  for(int LID=0 ; LID<numOwnedPoints ; ++LID){
    int GID = oneDimensionalOverlapMap->GID(LID);
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    unsigned int numEntries = 3*(numNeighbors+1);
    if(globalIndices.size() < numEntries)
      globalIndices.resize(numEntries);
    globalIndices[0] = 3*GID;
    globalIndices[1] = 3*GID + 1;
    globalIndices[2] = 3*GID + 2;
    for(int j=0 ; j<numNeighbors ; ++j){
      int neighborLocalID = neighborhoodList[neighborhoodListIndex++];
      int neighborGlobalID = oneDimensionalOverlapMap->GID(neighborLocalID);
      globalIndices[3*j+3] = 3*neighborGlobalID;
      globalIndices[3*j+4] = 3*neighborGlobalID + 1;
      globalIndices[3*j+5] = 3*neighborGlobalID + 2;
    }

    // The entries going into the tangent are a dense matrix of size numEntries by numEntries.
    // Each global ID in the list interacts with all other global IDs in the list.
    for(unsigned int i=0 ; i<numEntries ; ++i){
      for(unsigned int j=0 ; j<numEntries ; ++j)
	rowEntries[globalIndices[i]].insert(globalIndices[j]);
    }
  }

  // Allocate space in the global Epetra_FECrsMatrix
  vector<int> indices;
  vector<double> zeros;
  for(map<int, set<int> >::iterator rowEntry=rowEntries.begin(); rowEntry!=rowEntries.end() ; ++rowEntry){
    unsigned int numRowNonzeros = rowEntry->second.size();
    if(zeros.size() < numRowNonzeros)
      zeros.resize(numRowNonzeros, 0.0);

    // Load indices into a sorted vector
    indices.resize(numRowNonzeros);
    int i=0;
    for(set<int>::const_iterator globalIndex=rowEntry->second.begin() ; globalIndex!=rowEntry->second.end() ; ++globalIndex)
      indices[i++] = *globalIndex;
    sort(indices.begin(), indices.end());

    // Allocate space in the global matrix
    int err = tangent->InsertGlobalValues(rowEntry->first, numRowNonzeros, (const double*)&zeros[0], (const int*)&indices[0]);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err < 0, "**** PeridigmNS::Peridigm::allocateJacobian(), InsertGlobalValues() returned negative error code.\n");

    rowEntry->second.clear();
  }
  int err = tangent->GlobalAssemble();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::allocateJacobian(), GlobalAssemble() returned nonzero error code.\n");

  // create the serial Jacobian
  overlapJacobian = Teuchos::rcp(new PeridigmNS::SerialMatrix(tangent));
  workset->jacobian = overlapJacobian;
}

double PeridigmNS::Peridigm::computeQuasiStaticResidual(Teuchos::RCP<Epetra_Vector> residual) {

  PeridigmNS::Timer::self().startTimer("Compute Residual");

  // The residual is computed as the L2 norm of the internal force vector with the
  // entries corresponding to kinematic BC zeroed out.

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // Update forces based on new positions
  PeridigmNS::Timer::self().startTimer("Internal Force");
  modelEvaluator->evalModel(workset);
  PeridigmNS::Timer::self().stopTimer("Internal Force");

  // Copy force from the data manager to the mothership vector
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  force->PutScalar(0.0);
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    scratch->PutScalar(0.0);
    blockIt->exportData(*scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
    force->Update(1.0, *scratch, 1.0);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
  scratch->PutScalar(0.0);

  // copy the internal force to the residual vector
  // note that due to restrictions on CrsMatrix, these vectors have different (but equivalent) maps
  TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != force->MyLength(), "**** PeridigmNS::Peridigm::computeQuasiStaticResidual() incompatible vector lengths!\n");
  for(int i=0 ; i<force->MyLength() ; ++i)
    (*residual)[i] = (*force)[i];

  // convert force density to force
  for(int i=0 ; i<residual->MyLength() ; ++i)
    (*residual)[i] *= (*volume)[i/3];

  // zero out the rows corresponding to kinematic boundary conditions and compute the residual
  boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);
  double residualNorm2;
  residual->Norm2(&residualNorm2);
  double residualNormInf;
  residual->NormInf(&residualNormInf);

  PeridigmNS::Timer::self().stopTimer("Compute Residual");

  return residualNorm2 + 20.0*residualNormInf;
}

void PeridigmNS::Peridigm::computeImplicitJacobian(double beta) {

  // MLP: Pass this in from integrator routine
  double dt = *(workset->timeStep);
  double dt2 = dt*dt;

  // Compute the tangent
  tangent->PutScalar(0.0);
  PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
  modelEvaluator->evalJacobian(workset);
  int err = tangent->GlobalAssemble();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::computeImplicitJacobian(), GlobalAssemble() returned nonzero error code.\n");
  PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");

  // Code to symmeterize Jacobian

  // First, construct transpose
  bool makeDataContiguous = true;
  EpetraExt::RowMatrix_Transpose transposer( makeDataContiguous );
  Epetra_CrsMatrix& transTangent = dynamic_cast<Epetra_CrsMatrix&>(transposer(*tangent));

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

  // tangent = M - beta*dt*dt*K
  tangent->Scale(-beta*dt2);

  Epetra_Vector diagonal1(tangent->RowMap());
  Epetra_Vector diagonal2(tangent->RowMap());
  tangent->ExtractDiagonalCopy(diagonal1);
  for(int i=0 ; i<diagonal2.MyLength() ; ++i)
    diagonal2[i] = (*density)[i/3];
  diagonal1.Update(1.0, diagonal2, 1.0);
  tangent->ReplaceDiagonalValues(diagonal1);

}

void PeridigmNS::Peridigm::synchDataManagers() {
  // Need to ensure these primal fields are synchronized: VOLUME, BLOCK_ID, COORD3D, DISPL3D, CURCOORD3D, VELOC3D, FORCE_DENSITY3D, CONTACT_FORCE_DENSITY_3D

  // Copy data from mothership vectors to overlap vectors in blocks
  // VOLUME and BLOCK_ID are synched during creation and rebalance, and otherwise never changes
  // COORD3D is synched during creation and rebalance, and otherwise never changes
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*force, forceDensityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*contactForce, contactForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
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
  density = Teuchos::rcp((*oneDimensionalMothership)(2), false);         // density

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
  scratch = Teuchos::rcp((*threeDimensionalMothership)(8), false);       // scratch space

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
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*volume, volumeFieldId, PeridigmField::STEP_NONE, Insert);
    blockIt->importData(*blockIDs, blockIdFieldId, PeridigmField::STEP_NONE, Insert);
  }

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
  UTILITIES::Array<int> myGlobalIDs(myNumElements);
  int* myGlobalIDsPtr = myGlobalIDs.get();
  int* gIDs = oneDimensionalMap->MyGlobalElements();
  memcpy(myGlobalIDsPtr, gIDs, myNumElements*sizeof(int));
  decomp.myGlobalIDs = myGlobalIDs.get_shared_ptr();

  // fill myX
  // use current positions for x
  UTILITIES::Array<double> myX(myNumElements*dimension);
  double* myXPtr = myX.get();
  double* yPtr;
  y->ExtractView(&yPtr);
  memcpy(myXPtr, yPtr, myNumElements*dimension*sizeof(double));
  decomp.myX = myX.get_shared_ptr();

  // fill cellVolume
  UTILITIES::Array<double> cellVolume(myNumElements);
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
    TEUCHOS_TEST_FOR_EXCEPTION(localID == -1, Teuchos::RangeError, "Invalid index into rebalancedOneDimensionalOverlapMap");
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
        TEUCHOS_TEST_FOR_EXCEPTION(localNeighborID == -1, Teuchos::RangeError, "Invalid index into rebalancedOneDimensionalOverlapMap");
        neighborhoodList[neighborhoodIndex++] = localNeighborID;
      }
    }
    else{
      neighborhoodList[neighborhoodIndex++] = 0;
    }
  }

  return rebalancedNeighborhoodData;
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
		TEUCHOS_TEST_FOR_EXCEPTION(localID == -1, Teuchos::RangeError, "Invalid index into rebalancedOneDimensionalOverlapMap");
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
		TEUCHOS_TEST_FOR_EXCEPTION(contactNeighborGlobalIDs->count(globalID) == 0, Teuchos::RangeError, "Invalid index into contactNeighborGlobalIDs");
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
