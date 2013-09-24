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

#include <boost/unordered_set.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "Peridigm_Field.hpp"
#include "Peridigm_InfluenceFunction.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "Peridigm_PartialVolumeCalculator.hpp"
#include "Peridigm_OutputManager_ExodusII.hpp"
#include "Peridigm_SolverManager.hpp"
#include "Peridigm_SolverManagerContainer.hpp"
#include "Peridigm_ComputeManager.hpp"
#include "Peridigm_BoundaryAndInitialConditionManager.hpp"
#include "Peridigm_CriticalTimeStep.hpp"
#include "Peridigm_RandomNumber.hpp"
#include "Peridigm_Timer.hpp"
#include "materials/Peridigm_MaterialFactory.hpp"
#include "damage/Peridigm_DamageModelFactory.hpp"
#include "muParser/muParser.h"
#include "muParser/muParserPeridigmFunctions.h"
#include "Peridigm.hpp"

#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Epetra_RowMatrixTransposer.h>
#include <Ifpack.h>
#include <Ifpack_IC.h>
#include <Teuchos_VerboseObject.hpp>

using namespace std;

PeridigmNS::Peridigm::Peridigm(Teuchos::RCP<const Epetra_Comm> comm,
                               Teuchos::RCP<Teuchos::ParameterList> params)
  : analysisHasContact(false),
    analysisHasPartialVolumes(false),
    blockIdFieldId(-1),
    volumeFieldId(-1),
    modelCoordinatesFieldId(-1),
    coordinatesFieldId(-1),
    displacementFieldId(-1),
    velocityFieldId(-1),
    accelerationFieldId(-1),
    deltaTemperatureFieldId(-1),
    forceDensityFieldId(-1),
    contactForceDensityFieldId(-1),
    externalForceDensityFieldId(-1),
    partialVolumeFieldId(-1)
{
  peridigmComm = comm;
  peridigmParams = params;
  // set the comm for memory use statistics
  Memstat * memstat = Memstat::Instance();
  memstat->setComm(comm);

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // Seed random number generator for reproducable results
  seed_rand_num( 42 );

  // Initialize the influence function
  string influenceFunctionString = peridigmParams->sublist("Discretization").get<string>("Influence Function", "One");
  PeridigmNS::InfluenceFunction::self().setInfluenceFunction( influenceFunctionString );

  // Read mesh from disk or generate using geometric primatives.
  Teuchos::RCP<Teuchos::ParameterList> discParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Discretization", true) );

  // \todo When using partial volumes, the horizon should be increased by a value equal to the largest element dimension in the model (largest element diagonal for hexes).
  //       For an initial test, just double the horizon and hope for the best.
  // if(analysisHasPartialVolumes)
  //   discParams->set("Search Horizon", 2.0*discParams->get<double>("Horizon"));
  // else
  //   discParams->set("Search Horizon", discParams->get<double>("Horizon"));

  // The horizon may no longer be specified in the discretization block
  string msg = "\n**** Error, \"Horizon\" is no longer an allowable Discretization parameter.\n";
  msg +=         "****        A horizon for each block must be specified in the Blocks section.\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(discParams->isParameter("Horizon"), msg);

  // Create a list of horizon values for each block and add it to the discretization parameter list
  Teuchos::ParameterList& blockParams = peridigmParams->sublist("Blocks", true);
  map<string, double> blockHorizonValues = parseHorizonValuesFromBlockParameters(blockParams);
  for(map<string, double>::const_iterator it = blockHorizonValues.begin() ; it != blockHorizonValues.end() ; it++){
    string name = "Horizon " + it->first;
    discParams->set(name, it->second);
  }

  DiscretizationFactory discFactory(discParams);
  Teuchos::RCP<Discretization> peridigmDisc = discFactory.create(peridigmComm);
  checkHorizon(peridigmDisc, blockHorizonValues);
  initializeDiscretization(peridigmDisc);

  // If the user did not provide a finite-difference probe length, use a fraction of the minimum element radius
  double minElementRadius = peridigmDisc->getMinElementRadius();
  Teuchos::RCP<Teuchos::ParameterList> solverParams = sublist(peridigmParams, "Solver");
  if(!solverParams->isParameter("Finite Difference Probe Length"))
    solverParams->set("Finite Difference Probe Length", 1.0e-6*minElementRadius);

  // Instantiate and initialize the boundary and initial condition manager
  Teuchos::RCP<Teuchos::ParameterList> bcParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Boundary Conditions") );
  boundaryAndInitialConditionManager =
    Teuchos::RCP<BoundaryAndInitialConditionManager>(new BoundaryAndInitialConditionManager(*bcParams, this));

  boundaryAndInitialConditionManager->initialize(peridigmDisc);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  elementIdFieldId                   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id");
  blockIdFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Block_Id");
  volumeFieldId                      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");
  modelCoordinatesFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  coordinatesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  displacementFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement");
  velocityFieldId                    = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  accelerationFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration");
  deltaTemperatureFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
  forceDensityFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  contactForceDensityFieldId         = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density");
  externalForceDensityFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "External_Force_Density");
  if(analysisHasPartialVolumes)
    partialVolumeFieldId             = fieldManager.getFieldId(PeridigmField::BOND,    PeridigmField::SCALAR, PeridigmField::CONSTANT, "Partial_Volume");

  // Create field ids that may be required for output
  fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num");

  // Instantiate the contact manager
  Teuchos::ParameterList contactParams;
  if(peridigmParams->isSublist("Contact")){
    analysisHasContact = true;
    contactParams = peridigmParams->sublist("Contact");
    checkContactSearchRadius(contactParams,peridigmDisc);
    contactManager =
      Teuchos::RCP<ContactManager>(new ContactManager(contactParams, peridigmDisc, peridigmParams));
    contactManager->initialize(oneDimensionalMap,
                               threeDimensionalMap,
                               oneDimensionalOverlapMap,
                               bondMap,
                               globalNeighborhoodData,
                               peridigmDisc->getBlockID(),
                               blockHorizonValues);
    // contactManager->loadNeighborhoodData(globalNeighborhoodData,
    //                                      oneDimensionalMap,
    //                                      oneDimensionalOverlapMap);
    contactManager->loadAllMothershipData(blockIDs,
                                          volume,
                                          y, 
                                          v);
    contactManager->initializeContactBlocks();

    const std::string statTag = "Contact Initialized";
    memstat->addStat(statTag);
  }

  // Instantiate the blocks
  initializeBlocks(peridigmDisc);

  // Obtain parameter lists and factories for material models ane damage models
  // Material models
  Teuchos::ParameterList materialParams = peridigmParams->sublist("Materials", true);
  MaterialFactory materialFactory;
  // Damage models
  Teuchos::ParameterList damageModelParams;
  if(peridigmParams->isSublist("Damage Models"))
    damageModelParams = peridigmParams->sublist("Damage Models");
  DamageModelFactory damageModelFactory;

  // Associate material models and damage models with blocks
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){

    // Obtain the horizon for this block
    string blockName = blockIt->getName();
    double blockHorizon(0.0);
    if(blockHorizonValues.find(blockName) != blockHorizonValues.end())
      blockHorizon = blockHorizonValues[blockName];
    else if(blockHorizonValues.find("default") != blockHorizonValues.end())
      blockHorizon = blockHorizonValues["default"];
    else{
      string msg = "\n**** Error, no Horizon parameter supplied for block " + blockName + " and no default block parameter list provided.\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
    }

    // Set the material model
    string materialName = blockIt->getMaterialName();
    Teuchos::ParameterList matParams = materialParams.sublist(materialName);

    // Assign the horizon to the material model
    // Make sure the user did not try to set the horizon in the material block
    TEUCHOS_TEST_FOR_EXCEPT_MSG(matParams.isParameter("Horizon") , "\n**** Error, Horizon is an invalid material parameter.\n");
    matParams.set("Horizon", blockHorizon);

    // Assign the finite difference probe length
    double finiteDifferenceProbeLength = peridigmParams->sublist("Solver", true).get<double>("Finite Difference Probe Length");
    if(!matParams.isParameter("Finite Difference Probe Length"))
      matParams.set("Finite Difference Probe Length", finiteDifferenceProbeLength);

    // Instantiate the material model for this block
    Teuchos::RCP<const PeridigmNS::Material> materialModel = materialFactory.create(matParams);
    blockIt->setMaterialModel(materialModel);

    // Set the damage model (if any)
    string damageModelName = blockIt->getDamageModelName();
    if(damageModelName != "None"){
      Teuchos::ParameterList damageParams = damageModelParams.sublist(damageModelName, true);
      Teuchos::RCP<const PeridigmNS::DamageModel> damageModel = damageModelFactory.create(damageParams);
      blockIt->setDamageModel(damageModel);
    }
  }

  // Instantiate compute manager
  instantiateComputeManager();

  // Load the auxiliary field ids into the blocks (they will be
  // combined with material model and damage model ids when allocating
  // space in the Block's DataManager)
  vector<int> auxiliaryFieldIds;

  // Force the allocation of space in all blocks for the following field data
  // \todo Replace this with a query to the requested output fields, which is why we're forcing this allocation.
  auxiliaryFieldIds.push_back(modelCoordinatesFieldId);
  auxiliaryFieldIds.push_back(coordinatesFieldId);
  auxiliaryFieldIds.push_back(displacementFieldId);
  auxiliaryFieldIds.push_back(velocityFieldId);
  auxiliaryFieldIds.push_back(externalForceDensityFieldId);
  if(analysisHasContact)
    auxiliaryFieldIds.push_back(contactForceDensityFieldId);

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

  // apply initial conditions
  PeridigmNS::Timer::self().startTimer("Apply Initial Conditions");
  boundaryAndInitialConditionManager->applyInitialConditions();
  PeridigmNS::Timer::self().stopTimer("Apply Initial Conditions");

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
  
  // Initialize solver manager
  initializeSolverManager();

  // Initialize output manager
  initializeOutputManager();

  // Call rebalance function if analysis has contact
  // this is required to set up proper contact neighbor list
  if(analysisHasContact)
    contactManager->rebalance(0);

  // Create service manager
  serviceManager = Teuchos::rcp(new PeridigmNS::ServiceManager());
  serviceManager->requestService(computeManager->Services());

  // Check if request for allocation request of tangent matrix
  bool allocate_tangent = false;
  for (Teuchos::ParameterList::ConstIterator it = peridigmParams->begin(); it != peridigmParams->end(); ++it) {
    // See if name of parameterlist entry contains "Solver"
    const std::string output("Solver");
    const std::string name(it->first);
    size_t found = name.find(output);
    if (found!=std::string::npos) {
      Teuchos::RCP<Teuchos::ParameterList> slvrParams = sublist(peridigmParams, name);
      if(slvrParams->isSublist("QuasiStatic") ||
          slvrParams->isSublist("NOXQuasiStatic") ||
          slvrParams->isSublist("Implicit"))
        allocate_tangent = true;
    }
  }

  // Perform requested services
  if (serviceManager->isRequested(PeridigmNS::PeridigmService::ALLOCATE_TANGENT) || allocate_tangent) {
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
  }

  // Check if request for allocation of block diagonal tangent stiffness matrix
  if (serviceManager->isRequested(PeridigmNS::PeridigmService::ALLOCATE_BLOCK_DIAGONAL_TANGENT)) {
    // Allocate memory for non-zeros in global tangent and lock in the structure
    if(peridigmComm->MyPID() == 0 && !allocate_tangent){
      cout << "Allocating global block diagonal tangent matrix...";
      cout.flush();
    }
    PeridigmNS::Timer::self().startTimer("Allocate Global Block Diagonal Tangent");
    allocateBlockDiagonalJacobian();
    PeridigmNS::Timer::self().stopTimer("Allocate Global Block Diagonal Tangent");
    if(peridigmComm->MyPID() == 0 && !allocate_tangent){
      cout << "\n  number of rows = " << blockDiagonalTangent->NumGlobalRows() << endl;
      cout << "  number of nonzeros = " << blockDiagonalTangent->NumGlobalNonzeros() << "\n" << endl;
    }
  }

  // Set default value for current time;
  timeCurrent = 0.0;
}

void PeridigmNS::Peridigm::checkHorizon(Teuchos::RCP<Discretization> peridigmDisc, map<string, double> & blockHorizonValues){
  // Warn the user if it appears the horizon is too large and may cause memory to fill
  const double warningPerc = 0.80; // TODO: This might need to be adjusted
  const int mesh_size = peridigmDisc->getNumElem();
  const bool mesh_size_large = mesh_size > 10000; // TODO This might need to be adjusted
  const double maxPerc = (double)peridigmDisc->getMaxNumBondsPerElem() / (double)mesh_size;
  if(maxPerc >= warningPerc && mesh_size_large){
    if(peridigmComm->MyPID() == 0){
      cout << "** Warning: elements were detected with large neighborhood sizes relative to the number of local elements.\n"
           << "** The largest element neighborhood contains " << maxPerc * 100 << "% of the elements on this processor.\n"
           << "** This may indicate the horizon was selected too large and could lead to memory capacity being exceeded.\n\n";
    }
    return;
  }
  // Warn the user if a block's horizon is too large compared to the max element diameter
  const double warningSizeRatio = 25.0; // TODO: This might need to be adjusted
  string blockName = "";
  double horizonValue = 0.0;
  const double maxRad = peridigmDisc->getMaxElementRadius();
  for(map<string, double>::const_iterator it = blockHorizonValues.begin() ; it != blockHorizonValues.end() ; it++){
    if(it->second > warningSizeRatio * maxRad && mesh_size_large){
      blockName = it->first;
      horizonValue = it->second;
      if(peridigmComm->MyPID() == 0){
        cout << "** Warning: The horizon for " << blockName << " is " << horizonValue << ", which is\n"
             << "** more than " << warningSizeRatio <<  " times the max element radius (" << maxRad << ").\n"
             << "** This may indicate the horizon was selected too large\n"
             << "** and could lead to memory capacity being exceeded.\n\n";
      }
    }
  }
}

void PeridigmNS::Peridigm::checkContactSearchRadius(const Teuchos::ParameterList& contactParams, Teuchos::RCP<Discretization> peridigmDisc){
  if(!contactParams.isParameter("Search Radius"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Radius\" not specified.");

  const double maxRatio = 10.0;  // TODO: this might need to be adjusted
  const double contactRad = contactParams.get<double>("Search Radius");
  const double maxRad = peridigmDisc->getMaxElementRadius();

  if(contactRad/maxRad >= maxRatio){
    if(peridigmComm->MyPID() == 0){
      cout << "** Warning: the selected contact search radius, " << contactRad << ", is large\n"
           << "** relative to the maximum element diameter (" << maxRad << ").\n"
           << "** This may lead to the memory capacity being exceeded.\n\n";
    }
  }
}

map<string, double> PeridigmNS::Peridigm::parseHorizonValuesFromBlockParameters(Teuchos::ParameterList& blockParams) {

  map<string, double> blockHorizonValues;

  // Find the horizon value for each block and record the default horizon value (if any)
  for(Teuchos::ParameterList::ConstIterator it = blockParams.begin() ; it != blockParams.end() ; it++){
    const string& name = it->first;
    Teuchos::ParameterList& params = blockParams.sublist(name);
    string blockNamesString = params.get<string>("Block Names");
    // Parse space-delimited list of block names
    istringstream iss(blockNamesString);
    vector<string> blockNames;
    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter<vector<string> >(blockNames));
    for(vector<string>::const_iterator it = blockNames.begin() ; it != blockNames.end() ; ++it){
      double horizonValue = params.get<double>("Horizon");
      if( *it == "Default" || *it == "default")
        blockHorizonValues["default"] = horizonValue;
      else
        blockHorizonValues[*it] = horizonValue;
    }
  }
  
  return blockHorizonValues;
}

void PeridigmNS::Peridigm::initializeDiscretization(Teuchos::RCP<Discretization> peridigmDisc) {

  // oneDimensionalMap
  // used for cell volumes and scalar constitutive data
  oneDimensionalMap = peridigmDisc->getGlobalOwnedMap(1); 

  // oneDimensionalOverlapMap
  // used for initializing tangent structure
  // includes ghosts
  oneDimensionalOverlapMap = peridigmDisc->getGlobalOverlapMap(1);

  // threeDimensionalMap
  // used for positions, displacements, velocities and vector constitutive data
  threeDimensionalMap = peridigmDisc->getGlobalOwnedMap(3);

  // bondConstitutiveDataMap
  // a non-overlapping map used for storing constitutive data on bonds
  bondMap = peridigmDisc->getGlobalBondMap();

  // Create mothership vectors

  // \todo Do not allocate space for deltaTemperature if not needed.
  oneDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*oneDimensionalMap, 4));
  blockIDs = Teuchos::rcp((*oneDimensionalMothership)(0), false);         // block ID
  volume = Teuchos::rcp((*oneDimensionalMothership)(1), false);           // cell volume
  density = Teuchos::rcp((*oneDimensionalMothership)(2), false);          // density
  deltaTemperature = Teuchos::rcp((*oneDimensionalMothership)(3), false); // change in temperature

  // \todo Do not allocate space for the contact force nor deltaU if not needed.
  threeDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*threeDimensionalMap, 10));
  x = Teuchos::rcp((*threeDimensionalMothership)(0), false);             // initial positions
  u = Teuchos::rcp((*threeDimensionalMothership)(1), false);             // displacement
  y = Teuchos::rcp((*threeDimensionalMothership)(2), false);             // current positions
  v = Teuchos::rcp((*threeDimensionalMothership)(3), false);             // velocities
  a = Teuchos::rcp((*threeDimensionalMothership)(4), false);             // accelerations
  force = Teuchos::rcp((*threeDimensionalMothership)(5), false);         // force
  contactForce = Teuchos::rcp((*threeDimensionalMothership)(6), false);  // contact force (used only for contact simulations)
  externalForce = Teuchos::rcp((*threeDimensionalMothership)(7), false); // external force
  deltaU = Teuchos::rcp((*threeDimensionalMothership)(8), false);        // increment in displacement (used only for implicit time integration)
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

void PeridigmNS::Peridigm::initializeWorkset() {
  workset = Teuchos::rcp(new PHPD::Workset);
  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  *timeStep = 0.0;
  workset->timeStep = timeStep;
  workset->jacobian = overlapJacobian;
  workset->blocks = blocks;
  if(!contactManager.is_null())
    workset->contactManager = contactManager;
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

  // Initialize the parameterlist containing global Peridigm data (not stored in blocks)
  computeClassGlobalData = Teuchos::rcp(new Teuchos::ParameterList());
  Teuchos::RCP<Epetra_FECrsMatrix> *tmp1 = &( tangent );
  Teuchos::RCP<Epetra_FECrsMatrix> *tmp2 = &( blockDiagonalTangent );
  Teuchos::RCP<PeridigmNS::SerialMatrix> *tmp3 = &( overlapJacobian );
  Teuchos::RCP<Epetra_Map> *tmp4 = &( blockDiagonalTangentMap );
  computeClassGlobalData->set("tangent",tmp1);
  computeClassGlobalData->set("blockDiagonalTangent",tmp2);
  computeClassGlobalData->set("overlapJacobian",tmp3);
  computeClassGlobalData->set("blockDiagonalTangentMap",tmp4);

  computeManager = Teuchos::rcp( new PeridigmNS::ComputeManager( computeParams, peridigmComm, computeClassGlobalData ) );
}

void PeridigmNS::Peridigm::initializeBlocks(Teuchos::RCP<Discretization> disc) {

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

void PeridigmNS::Peridigm::initializeSolverManager() {

  // Create empty container for solver managers
  solverManager = Teuchos::rcp(new PeridigmNS::SolverManagerContainer() );

  Teuchos::RCP<Teuchos::ParameterList> solverParams;

  // Loop over high level parameter list entries to find all solver lists
  for (Teuchos::ParameterList::ConstIterator it = peridigmParams->begin(); it != peridigmParams->end(); ++it) {
    // See if name of parameterlist entry contains "Solver".
    const std::string output("Solver");
    const std::string name(it->first);
    size_t found = name.find(output);
    if (found!=std::string::npos) {
      // Make copy of list
      try{
        solverParams = Teuchos::rcp( new Teuchos::ParameterList( peridigmParams->sublist(name,true) ) );
      }
      catch(const std::exception &e){
        string msg = "Peridigm::initializeSolverManager: ";
        msg+= name;
        msg+= " is not a Teuchos::ParameterList sublist.";
        TEUCHOS_TEST_FOR_EXCEPT_MSG( true, msg );
      }

      solverManager->add( Teuchos::rcp(new PeridigmNS::SolverManager( solverParams, this) ) );
    }
  }

}

void PeridigmNS::Peridigm::execute(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  TEUCHOS_TEST_FOR_EXCEPT_MSG(solverParams.is_null(), "Error in Peridigm::execute, solverParams is null.\n");
 
  // allowable explicit time integration schemes:  Verlet
  if(solverParams->isSublist("Verlet")){
    executeExplicit(solverParams);}

  // allowable implicit time integration schemes:  Implicit, QuasiStatic
  else if(solverParams->isSublist("QuasiStatic"))    
    executeQuasiStatic(solverParams);
  else if(solverParams->isSublist("NOXQuasiStatic"))    
    executeNOXQuasiStatic(solverParams);
  else if(solverParams->isSublist("Implicit"))    
    executeImplicit(solverParams);

  PeridigmNS::Memstat * memstat = PeridigmNS::Memstat::Instance();
  const std::string statTag = "Post Execute";
  memstat->addStat(statTag);
}

void PeridigmNS::Peridigm::executeSolvers() {
  // Call the solver manager to execute all solvers in sequence
  solverManager->executeSolvers();

}

void PeridigmNS::Peridigm::executeExplicit(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

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

  // Pointer index into sub-vectors for use with BLAS
  double *xPtr, *uPtr, *yPtr, *vPtr, *aPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  a->ExtractView( &aPtr );
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
    blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
  }
  if(analysisHasContact)
    contactManager->importData(volume, y, v);
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // \todo The velocity copied into the DataManager is actually the midstep velocity, not the NP1 velocity; this can be fixed by creating a midstep velocity field in the DataManager and setting the NP1 value as invalid.

  // Evaluate internal force and contact force in initial configuration for use in first timestep
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
  if(analysisHasContact){
    contactManager->exportData(contactForce);
    force->Update(1.0, *contactForce, 1.0);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // evaluate the external (body) forces:
  boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,0.0); // external forces are dirichlet BCs so the previous time is defaulted to 0.0

  // fill the acceleration vector
  (*a) = (*force);
  for(int i=0 ; i<a->MyLength() ; ++i)
  {
    (*a)[i] += (*externalForce)[i];
    (*a)[i] /= (*density)[i/3];
  }
  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  synchDataManagers();
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
    PeridigmNS::Timer::self().startTimer("Rebalance");
    // \todo Should we load updated information first?  If so, only do this if we're really going to rebalance.
    if(analysisHasContact)
      contactManager->rebalance(step);
    PeridigmNS::Timer::self().stopTimer("Rebalance");

    // Do one step of velocity-Verlet

    // V^{n+1/2} = V^{n} + (dt/2)*A^{n}
    // blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt2, aPtr, vPtr, 1, 1);

    // Set the velocities for dof with kinematic boundary conditions.
    // This will propagate through the Verlet integrator and result in the proper
    // displacement boundary conditions on y and consistent values for v and u.
    PeridigmNS::Timer::self().startTimer("Apply kinematic B.C.");
    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply kinematic B.C.");

    // evaluate the external (body) forces:
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,0.0); // external forces are dirichlet BCs so the previous time is defaulted to 0.0

    // Y^{n+1} = X_{o} + U^{n} + (dt)*V^{n+1/2}
    // \todo Replace with blas call
    for(int i=0 ; i<y->MyLength() ; ++i)
      yPtr[i] = xPtr[i] + uPtr[i] + dt*vPtr[i];

    // U^{n+1} = U^{n} + (dt)*V^{n+1/2}
    // blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt, vPtr, uPtr, 1, 1);

    // \todo The velocity copied into the DataManager is actually the midstep velocity, not the NP1 velocity; this can be fixed by creating a midstep velocity field in the DataManager and setting the NP1 value as invalid.

    // Copy data from mothership vectors to overlap vectors in data manager
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
    }
    if(analysisHasContact)
      contactManager->importData(volume, y, v);

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

    // Check for NaNs in force evaluation
    // We'd like to know now because a NaN will likely cause a difficult-to-unravel crash downstream.
    for(int i=0 ; i<externalForce->MyLength() ; ++i)
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*externalForce)[i]), "**** NaN returned by external force evaluation.\n");

    if(analysisHasContact){
      contactManager->exportData(contactForce);
      // Check for NaNs in contact force evaluation
      for(int i=0 ; i<contactForce->MyLength() ; ++i)
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*contactForce)[i]), "**** NaN returned by contact force evaluation.\n");
      // Add contact forces to forces
      force->Update(1.0, *contactForce, 1.0);
    }

    // fill the acceleration vector
    (*a) = (*force);
    for(int i=0 ; i<a->MyLength() ; ++i)
    {
      (*a)[i] += (*externalForce)[i];
      (*a)[i] /= (*density)[i/3];
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt2, aPtr, vPtr, 1, 1);

    PeridigmNS::Timer::self().startTimer("Output");
    synchDataManagers();
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
    (*v)[i] = (*soln)[i] / (*(workset->timeStep));
  }
  v->Update(1.0, *noxVelocityAtDOFWithKinematicBC, 1.0); // Necessary because soln is zero at dof with kinematic bc

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
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
    double conditionNumberTolerance = 1.0e3;
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

void PeridigmNS::Peridigm::executeNOXQuasiStatic(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction = Teuchos::rcp(new Epetra_Vector(force->Map()));

  // Create vectors that are specific to NOX quasi-statics.
  Teuchos::RCP<Epetra_Vector> soln = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateView));
  soln->PutScalar(0.0);

  Teuchos::RCP<Epetra_Vector> initialGuess = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<NOX::Epetra::Vector> noxInitialGuess = Teuchos::rcp(new NOX::Epetra::Vector(initialGuess, NOX::Epetra::Vector::CreateView));
  initialGuess->PutScalar(0.0);
  
  noxVelocityAtDOFWithKinematicBC = Teuchos::rcp(new Epetra_Vector(v->Map()));

  // Initialize velocity to zero
  v->PutScalar(0.0);
  
  // Create a placeholder for the timestep within the workset that is passed to the model evaluator
  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;
 
  // Pointers into mothership vectors
  double *xPtr, *uPtr, *yPtr, *vPtr, *deltaUPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  deltaU->ExtractView( &deltaUPtr );

  // "Solver" parameter list
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
  synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  // Functionality for updating the Jacobian at a user-specified interval
  // This does not appear to be available for nonlinear CG in NOX, but it's important for peridynamics, so
  // we'll handle it directly within Peridigm
  m_noxTriggerJacobianUpdate = 1;
  if(noxQuasiStaticParams->sublist("Direction").sublist("Nonlinear CG").isParameter("Update Jacobian"))
    m_noxTriggerJacobianUpdate = noxQuasiStaticParams->sublist("Direction").sublist("Nonlinear CG").get<int>("Update Jacobian");

  Epetra_Time loadStepCPUTime(*peridigmComm);
  double cumulativeLoadStepCPUTime = 0.0;

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

    v->PutScalar(0.0); 

    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(initialGuess);

    // Apply the displacement increment and update the temperature field
    // Note that the soln vector was created to be compatible with the tangent matrix, and hence needed to be
    // constructed with an Epetra_Map.  The mothership vectors were constructed with an Epetra_BlockMap, and it's
    // this map that the boundary and intial condition manager expects.  So, make sure that the boundary and initial
    // condition manager gets the right type of vector.

    PeridigmNS::Timer::self().startTimer("Apply kinematic B.C.");
    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent, timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply kinematic B.C.");

    // For NOX, add the increment in displacement BC directly into the displacement vector
    for(int i=0 ; i<u->MyLength() ; ++i)
      uPtr[i] += deltaUPtr[i];

    // Note:  applyBoundaryConditions() sets the velocity as well as the displacement.
    //        This needs to be stored because NOX returns a deltaU of zero for these dof
    *noxVelocityAtDOFWithKinematicBC = *v;

    // evaluate the external (body) forces:
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent, timePrevious);

    *soln = *initialGuess;

    double toleranceMultiplier = 1.0;
    if(!useAbsoluteTolerance){
      // compute the vector of reactions, i.e., the forces corresponding to degrees of freedom for which kinematic B.C. are applied
      for(int i=0 ; i<y->MyLength() ; ++i)
        yPtr[i] = xPtr[i] + uPtr[i];
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
 
    // Print load step timing information
    double CPUTime = loadStepCPUTime.ElapsedTime();
    cumulativeLoadStepCPUTime += CPUTime;
    if(peridigmComm->MyPID() == 0)
      cout << setprecision(2) << "  cpu time for load step = " << CPUTime << " sec., cumulative cpu time = " << cumulativeLoadStepCPUTime << " sec.\n" << endl;

    // Store the velocity
    for(int i=0 ; i<v->MyLength() ; ++i)
      (*v)[i] = finalSolution[i]/timeIncrement;
    v->Update(1.0, *noxVelocityAtDOFWithKinematicBC, 1.0);

    // Add the converged displacement increment to the displacement
    for(int i=0 ; i<u->MyLength() ; ++i)
      (*u)[i] += finalSolution[i];

    // Write output for completed load step
    PeridigmNS::Timer::self().startTimer("Output");
    synchDataManagers();
    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->updateState();
  }
  if(peridigmComm->MyPID() == 0)
    cout << endl;
}

void PeridigmNS::Peridigm::executeQuasiStatic(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  // Create vectors that are specific to quasi-statics.
  // These must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction = Teuchos::rcp(new Epetra_Vector(force->Map()));

  // Vector for predictor
  Epetra_Vector predictor(v->Map());

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

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
  double *xPtr, *uPtr, *deltaUPtr, *yPtr, *vPtr, *aPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  deltaU->ExtractView( &deltaUPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  a->ExtractView( &aPtr );

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
  synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  Epetra_Time loadStepCPUTime(*peridigmComm);
  double cumulativeLoadStepCPUTime = 0.0;

  for(int step=1 ; step<(int)timeSteps.size() ; step++){

    loadStepCPUTime.ResetStartTime();

    double timePrevious = timeCurrent;
    timeCurrent = timeSteps[step];
    double timeIncrement = timeCurrent - timePrevious;
    *timeStep = timeIncrement;

    // Update nodal positions for nodes with kinematic B.C.
    deltaU->PutScalar(0.0);
    PeridigmNS::Timer::self().startTimer("Apply kinematic B.C.");
    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply kinematic B.C.");

    // evaluate the external (body) forces:
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);

    // Set the current position and velocity
    // \todo We probably want to rework this so that the material models get valid x, u, and y values.
    // Currently the u values are from the previous load step (and if we update u here we'll be unable
    // to properly undo a time step, which we'll need to adaptive time stepping).
    for(int i=0 ; i<y->MyLength() ; ++i){
      yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
      vPtr[i] = deltaUPtr[i]/timeIncrement;
    }

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
    bool usePreconditioner = false; // \todo Determine why ifpack preconditioners started exhibiting problems with Trilinos 11.2.5 (Jul-11-2013).
                                    //       For the record, Trilinos 11.2.4 (Jun-20-2013) works.
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
          (*lhs)[i] = predictor[i]*timeIncrement;
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
        alpha = quasiStaticsLineSearch(residual, lhs, timeIncrement);
        PeridigmNS::Timer::self().stopTimer("Line Search");

        // Apply increment to nodal positions
        for(int i=0 ; i<y->MyLength() ; ++i){
          deltaUPtr[i] += alpha*(*lhs)[i];
          yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
          vPtr[i] = deltaUPtr[i]/timeIncrement;
        }

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
    cumulativeLoadStepCPUTime += CPUTime;
    if(peridigmComm->MyPID() == 0)
      cout << setprecision(2) << "  cpu time for load step = " << CPUTime << " sec., cumulative cpu time = " << cumulativeLoadStepCPUTime << " sec.\n" << endl;

    // Add the converged displacement increment to the displacement
    for(int i=0 ; i<u->MyLength() ; ++i)
      (*u)[i] += (*deltaU)[i];

    // Store the velocity for use as a predictor in the next load step
    predictor = *v;

    // Write output for completed load step
    PeridigmNS::Timer::self().startTimer("Output");
    synchDataManagers();
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
                                                    Teuchos::RCP<Epetra_Vector> lhs,
                                                    double dt)
{
  Teuchos::RCP<Epetra_Vector> tempVector = Teuchos::rcp(new Epetra_Vector(*deltaU));

  double *lhsPtr, *residualPtr, *deltaUPtr, *xPtr, *yPtr, *uPtr, *vPtr;
  lhs->ExtractView(&lhsPtr);
  residual->ExtractView(&residualPtr);
  deltaU->ExtractView(&deltaUPtr);
  x->ExtractView(&xPtr);
  y->ExtractView(&yPtr);
  u->ExtractView(&uPtr);
  v->ExtractView(&vPtr);

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
  for(int i=0 ; i<y->MyLength() ; ++i){
    deltaUPtr[i] += epsilon*lhsPtr[i];
    yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
    vPtr[i] = deltaUPtr[i]/dt;
  }

  Teuchos::RCP<Epetra_Vector> perturbedResidual = Teuchos::rcp(new Epetra_Vector(*residual));
  computeQuasiStaticResidual(perturbedResidual);
  double SR, SPerturbedR;
  lhs->Dot(*residual, &SR);
  lhs->Dot(*perturbedResidual, &SPerturbedR);
  double tempAlpha = -1.0*epsilon*SR/(SPerturbedR - SR);
  if(tempAlpha > -0.1 && tempAlpha < 10.0)
    candidateAlphas.push_back(tempAlpha);
  *deltaU = *tempVector;

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
    for(int i=0 ; i<y->MyLength() ; ++i){
      deltaUPtr[i] += alpha*lhsPtr[i];
      yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
      vPtr[i] = deltaUPtr[i]/dt;
    }
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
    for(int i=0 ; i<y->MyLength() ; ++i){
      deltaUPtr[i] += alpha*lhsPtr[i];
      yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
      vPtr[i] = deltaUPtr[i]/dt;
    }
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

void PeridigmNS::Peridigm::executeImplicit(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  // Create vectors that are specific to implicit dynamics
  // The residual must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> un = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> vn = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> an = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  Teuchos::RCP<double> timeStep = Teuchos::rcp(new double);
  workset->timeStep = timeStep;

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
  double *xPtr, *uPtr, *yPtr, *vPtr, *aPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  a->ExtractView( &aPtr );

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
  synchDataManagers();
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
      blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
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

    // evaluate the external (body) forces:
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,0.0);

    // Compute the residual
    // residual = beta*dt*dt*(M*a - force)
    // Note that due to restrictions to CrsMatrix, the residual has a different (but equivalent) map
    // than the force and acceleration
    for(int i=0 ; i<residual->MyLength() ; ++i)
      (*residual)[i] = beta*dt2*( (*density)[i/3] * (*a)[i] - (*force)[i] - (*externalForce)[i]);

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
          cout << std::endl << "ERROR: Belos::LinearProblem failed to set up correctly!" << endl;
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
        blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
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
        (*residual)[i] = beta*dt2*( (*density)[i/3] * (*a)[i] - (*force)[i]  - (*externalForce)[i]);

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
    synchDataManagers();
    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->updateState();

    cout << endl;
  }
}

void PeridigmNS::Peridigm::allocateJacobian() {

  // do not re-allocate if already allocated
  if (tangent != Teuchos::null) return;

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
  map<int, boost::unordered_set<int> > rowEntries;

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
  for(map<int, boost::unordered_set<int> >::iterator rowEntry=rowEntries.begin(); rowEntry!=rowEntries.end() ; ++rowEntry){
    unsigned int numRowNonzeros = rowEntry->second.size();
    if(zeros.size() < numRowNonzeros)
      zeros.resize(numRowNonzeros, 0.0);

    // Load indices into a sorted vector
    indices.resize(numRowNonzeros);
    int i=0;
    for(boost::unordered_set<int>::const_iterator globalIndex=rowEntry->second.begin() ; globalIndex!=rowEntry->second.end() ; ++globalIndex)
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

  PeridigmNS::Memstat * memstat = PeridigmNS::Memstat::Instance();
  const std::string statTag = "Allocated Jacobian";
  memstat->addStat(statTag);
}

void PeridigmNS::Peridigm::allocateBlockDiagonalJacobian() {

  // do not re-allocate if already allocated
  if (blockDiagonalTangent != Teuchos::null) return;

  if(tangent != Teuchos::null) // the tangent matrix is alread allocated (i.e. this is an implicit or QS simulation) so use that one instead
  {
    blockDiagonalTangent = tangent;
    blockDiagonalTangentMap = tangentMap;
    return;
  }


  // err code
  int err;

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
  blockDiagonalTangentMap = Teuchos::rcp(new Epetra_Map(numGlobalElements, numMyElements, &myGlobalElements[0], indexBase, *peridigmComm));
  myGlobalElements.clear();

  // Create the global tangent matrix
  Epetra_DataAccess CV = Copy;
  int numEntriesPerRow = 3;  // Indicates allocation will take place during the insertion phase
  bool ignoreNonLocalEntries = false;
  blockDiagonalTangent = Teuchos::rcp(new Epetra_FECrsMatrix(CV, *blockDiagonalTangentMap, numEntriesPerRow, ignoreNonLocalEntries));

  // Store nonzero columns for each row, with everything in global indices
  map<int, boost::unordered_set<int> > rowEntries;

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
  vector<int> blockDiagonalIndices;
  blockDiagonalIndices.resize(3);
  vector<double> zeros;
  for(map<int, boost::unordered_set<int> >::iterator rowEntry=rowEntries.begin(); rowEntry!=rowEntries.end() ; ++rowEntry){
    unsigned int numRowNonzeros = rowEntry->second.size();
    if(zeros.size() < numRowNonzeros)
      zeros.resize(numRowNonzeros, 0.0);

    // Load indices into a sorted vector
    indices.resize(numRowNonzeros);
    int i=0;
    for(boost::unordered_set<int>::const_iterator globalIndex=rowEntry->second.begin() ; globalIndex!=rowEntry->second.end() ; ++globalIndex)
      indices[i++] = *globalIndex;
    sort(indices.begin(), indices.end());

    // Figure out node ID from global DOF ID
    int myID = rowEntry->first/3;
    blockDiagonalIndices[0] = 3*myID + 0;
    blockDiagonalIndices[1] = 3*myID + 1;
    blockDiagonalIndices[2] = 3*myID + 2;

    // Allocate space in the global matrix
    //err = blockDiagonalTangent->InsertGlobalValues(rowEntry->first, numRowNonzeros, (const double*)&zeros[0], (const int*)&indices[0]);
    err = blockDiagonalTangent->InsertGlobalValues(rowEntry->first, numEntriesPerRow, (const double*)&zeros[0], (const int*)&blockDiagonalIndices[0]);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err < 0, "**** PeridigmNS::Peridigm::allocateblockDiagonalJacobian(), InsertGlobalValues() returned negative error code.\n");

    rowEntry->second.clear();
  }
  err = blockDiagonalTangent->GlobalAssemble();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::allocateBlockDiagonalJacobian(), GlobalAssemble() returned nonzero error code.\n");

  // create the serial Jacobian
  overlapJacobian = Teuchos::rcp(new PeridigmNS::SerialMatrix(blockDiagonalTangent));
  workset->jacobian = overlapJacobian;

  PeridigmNS::Memstat * memstat = PeridigmNS::Memstat::Instance();
  const std::string statTag = "Alloc Blk Diag Jacobian";
  memstat->addStat(statTag);
}

double PeridigmNS::Peridigm::computeQuasiStaticResidual(Teuchos::RCP<Epetra_Vector> residual) {

  PeridigmNS::Timer::self().startTimer("Compute Residual");

  // The residual is computed as the norm of the internal force vector with the
  // entries corresponding to kinematic BC zeroed out.
  // The specific residual measure is the L2 norm plus twenty times the infinity norm

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
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

  // Check for NaNs in force evaluation
  // We'd like to know now because a NaN will likely cause a difficult-to-unravel crash downstream.
  for(int i=0 ; i<force->MyLength() ; ++i)
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*force)[i]), "**** NaN returned by force evaluation.\n");
  for(int i=0 ; i<externalForce->MyLength() ; ++i)
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!boost::math::isfinite((*externalForce)[i]), "**** NaN returned by external force evaluation.\n");

  // copy the internal force to the residual vector
  // note that due to restrictions on CrsMatrix, these vectors have different (but equivalent) maps
  TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != force->MyLength(), "**** PeridigmNS::Peridigm::computeQuasiStaticResidual() incompatible vector lengths!\n");
  for(int i=0 ; i<force->MyLength() ; ++i)
    (*residual)[i] = (*force)[i];

  TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != externalForce->MyLength(), "**** PeridigmNS::Peridigm::computeQuasiStaticResidual() incompatible vector lengths!\n");
  for(int i=0 ; i<externalForce->MyLength() ; ++i)
    (*residual)[i] += (*externalForce)[i];

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
  EpetraExt::RowMatrix_Transpose transposer;
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

  // Copy data from mothership vectors to overlap vectors in blocks
  // Volume and Block_Id are synched during creation and rebalance, and otherwise never changes
  // Model_Coordinates is synched during creation and rebalance, and otherwise never changes

  PeridigmNS::Timer::self().startTimer("Gather/Scatter");

  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(*u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*force, forceDensityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*contactForce, contactForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*externalForce, externalForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(*deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
  }

  // The hourglass force density is a special case.  It needs to be parallel assembled
  // prior to output.

  static Teuchos::RCP<Epetra_Vector> tempVector;  

  if(PeridigmNS::FieldManager::self().hasField("Hourglass_Force_Density")){
    int hourglassForceDensityFieldId = PeridigmNS::FieldManager::self().getFieldId("Hourglass_Force_Density");
    if(tempVector.is_null())
      tempVector = Teuchos::rcp(new Epetra_Vector(scratch->Map()));
    tempVector->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(*scratch, hourglassForceDensityFieldId, PeridigmField::STEP_NP1, Add);
      tempVector->Update(1.0, *scratch, 1.0);
    }
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->importData(*tempVector, hourglassForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
  }

  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
}

Teuchos::RCP< map< string, vector<int> > > PeridigmNS::Peridigm::getExodusNodeSets(){
  Teuchos::RCP< map< string, vector<int> > > nodeSets = boundaryAndInitialConditionManager->getNodeSets();
  Teuchos::RCP< map< string, vector<int> > > exodusNodeSets = Teuchos::rcp(new map< string, vector<int> >() );
  map< string, vector<int> >::iterator it;
  for(it=nodeSets->begin() ; it!=nodeSets->end() ; it++){
    const string& nodeSetName = it->first;
    const vector<int>& nodeSet = it->second;
    (*exodusNodeSets)[nodeSetName] = vector<int>(); // \todo Preallocate space, once we're sure the node sets are the right size
    vector<int>& exodusNodeSet = (*exodusNodeSets)[nodeSetName];
    for(unsigned int i=0 ; i<nodeSet.size() ; ++i){
      int localId = oneDimensionalMap->LID(nodeSet[i]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(localId == -1, "**** Error, Peridigm::getExodusNodeSets() encountered off-processor node in node set.\n");
      exodusNodeSet.push_back(localId + 1);
    }
  }
  return exodusNodeSets;
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
