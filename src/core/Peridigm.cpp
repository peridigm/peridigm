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
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <unordered_set>
#include <iterator>
#include <cmath>

#include "Peridigm_Field.hpp"
#include "Peridigm_HorizonManager.hpp"
#include "Peridigm_InfluenceFunction.hpp"
#include "Peridigm_DiscretizationFactory.hpp"
#include "Peridigm_OutputManager_ExodusII.hpp"
#include "Peridigm_ComputeManager.hpp"
#include "Peridigm_ContactModelFactory.hpp"
#include "Peridigm_BoundaryAndInitialConditionManager.hpp"
#include "Peridigm_DegreesOfFreedomManager.hpp"
#include "Peridigm_CriticalTimeStep.hpp"
#include "Peridigm_Timer.hpp"
#include "Peridigm_MaterialFactory.hpp"
#include "Peridigm_DamageModelFactory.hpp"
#include "Peridigm_InterfaceAwareDamageModel.hpp"
#include "Peridigm_UserDefinedTimeDependentCriticalStretchDamageModel.hpp"
#include "Peridigm_ShortRangeForceContactModel.hpp"
#include "Peridigm_UserDefinedTimeDependentShortRangeForceContactModel.hpp"
#include "Peridigm.hpp"
#include "correspondence.h" // For Invert3by3Matrix // TODO this should go
#include "Peridigm_DataManager.hpp" //For readBlocktoDisk & writeBlocktoDisk
#ifdef PERIDIGM_PV
  #include "Peridigm_PartialVolumeCalculator.hpp"
#endif

#include <Epetra_Import.h>
#include <Epetra_LinearProblem.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_Transpose_RowMatrix.h>
#include <Epetra_RowMatrixTransposer.h>
#include <Ifpack.h>
#include <Ifpack_IC.h>
#include <Teuchos_VerboseObject.hpp>

// required for restart
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"
#include <sys/stat.h>

using namespace std;

PeridigmNS::Peridigm::Peridigm(const MPI_Comm& comm,
                               Teuchos::RCP<Teuchos::ParameterList> params,
                               Teuchos::RCP<Discretization> inputPeridigmDiscretization)
  : agePeridigmPreconditioner(0),
    maxAgePeridigmPreconditioner(0),
    analysisHasContact(false),
    analysisHasDataLoader(false),
    analysisHasMultiphysics(false),
    computeIntersections(false),
    constructInterfaces(false),
    blockIdFieldId(-1),
    horizonFieldId(-1),
    volumeFieldId(-1),
    modelCoordinatesFieldId(-1),
    coordinatesFieldId(-1),
    displacementFieldId(-1),
    velocityFieldId(-1),
    accelerationFieldId(-1),
    temperatureFieldId(-1),
    concentrationFieldId(-1),
    deltaTemperatureFieldId(-1),
    fluxDivergenceFieldId(-1),
    concentrationFluxDivergenceFieldId(-1),
    forceDensityFieldId(-1),
    contactForceDensityFieldId(-1),
    externalForceDensityFieldId(-1),
    partialVolumeFieldId(-1),
    fluidPressureYFieldId(-1),
    fluidPressureUFieldId(-1),
    fluidPressureVFieldId(-1),
    fluidFlowDensityFieldId(-1),
    numMultiphysDoFs(0),
    analysisHasBondAssociatedHypoelasticModel(false),
    damageFieldId(-1),
    jacobianDeterminantFieldId(-1),
    weightedVolumeFieldId(-1),
    velocityGradientXFieldId(-1),
    velocityGradientYFieldId(-1),
    velocityGradientZFieldId(-1)
{
#ifdef HAVE_MPI
  peridigmComm = Teuchos::rcp(new Epetra_MpiComm(comm));
#else
  peridigmComm = Teuchos::rcp(new Epetra_SerialComm);
#endif
  if(peridigmComm->MyPID() == 0)
	  if(params->isParameter("Multiphysics") && params->isParameter("Restart") ){
		  TEUCHOS_TEST_FOR_EXCEPT_MSG((params->isParameter("Multiphysics") && params->isParameter("Restart") ), "Error: Restart for Multiphysics is not implemented yet.\n");
		  MPI_Finalize();
		  exit(0);
	  }
  peridigmParams = params;
  // set the comm for memory use statistics
  Memstat * memstat = Memstat::Instance();
  memstat->setComm(peridigmComm);

  // Tracker for recording the total number of iterations taken by the nonlinear solver
  nonlinearSolverIterations = Teuchos::rcp(new int);
  *nonlinearSolverIterations = 0;

  out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // Process and validate requests for multiphysics
  // Can the number of multiphysics DoFs be accomodated?
  string multiphysError;
  if(peridigmParams->isParameter("Multiphysics")){
    std::cout<<"\n**** Multiphysics is selected.\n"<< std::endl;
    if(peridigmParams->get<int>("Multiphysics") == 1){
      analysisHasMultiphysics = true;
      numMultiphysDoFs = peridigmParams->get<int>("Multiphysics");
      std::cout<<"\n**** Multiphysics is enabled, pending material model screening.\n"<< std::endl;
    }
    else{
      analysisHasMultiphysics = false;
      numMultiphysDoFs = 0;
      multiphysError = "\n**** Error, number of requested Multiphysics DoFs cannot be accomodated.\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, multiphysError);
    }
  }
  else{
    analysisHasMultiphysics = false;
    numMultiphysDoFs = 0;
  }

  // Initialize the influence function
  string influenceFunctionString = peridigmParams->sublist("Discretization").get<string>("Influence Function", "One");
  PeridigmNS::InfluenceFunction::self().setInfluenceFunction( influenceFunctionString );

  // Read mesh from disk or generate using geometric primatives.
  Teuchos::RCP<Teuchos::ParameterList> discParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Discretization", true) );

  // The horizon may no longer be specified in the discretization block
  // Throw an exception if the user is running an old input deck with the horizon in the discretization parameter list
  string msg = "\n**** Error, \"Horizon\" is no longer an allowable Discretization parameter.\n";
  msg +=         "****        A horizon for each block must be specified in the Blocks section.\n";
  TEUCHOS_TEST_FOR_EXCEPT_MSG(discParams->isParameter("Horizon"), msg);

  // Check for command to compute horizon-element intersections
  if(discParams->isParameter("Compute Element-Horizon Intersections"))
    computeIntersections = discParams->get<bool>("Compute Element-Horizon Intersections");
#ifndef PERIDIGM_PV
  TEUCHOS_TEST_FOR_EXCEPT_MSG(computeIntersections, "\n**** Error:  Horizon-Element intersections not enabled, recompile with -DUSE_PV.\n");
#endif

  // Pass the blockParams to the HorizonManager
  Teuchos::ParameterList& blockParams = peridigmParams->sublist("Blocks", true);
  PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();
  horizonManager.loadHorizonInformationFromBlockParameters(blockParams);

  // Create a list containing parameters for each solver
  for (Teuchos::ParameterList::ConstIterator it = peridigmParams->begin(); it != peridigmParams->end(); ++it) {
    // Check for string "Solver" in parameter list entry
    const std::string name(it->first);
    size_t found = name.find("Solver");
    if (found!=std::string::npos)
      solverParameters.push_back( sublist(peridigmParams, name) );
  }

  // For the case where multiple solvers are used, assume that the degrees of freedom (multiphysics) are the same for all solvers.
  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  Teuchos::ParameterList solverParamsForDofManager;
  if (solverParameters.size() > 0) {
    solverParamsForDofManager = *(solverParameters[0]);
  }
  dofManager.initialize(solverParamsForDofManager);
  if(peridigmComm->MyPID() == 0) {
    dofManager.print();
  }

  // Check solver parameters for request to allocate tangent matrix
  // Note that Peridigm can be run with multiple solvers, applied in sequence
  bool implicitTimeIntegration(false), userSpecifiedFullTangent(false), userSpecifiedBlockDiagonalTangent(false);
  for(unsigned int i=0 ; i<solverParameters.size() ; ++i){
    if(solverParameters[i]->isSublist("QuasiStatic") || solverParameters[i]->isSublist("NOXQuasiStatic") || solverParameters[i]->isSublist("Implicit")){
      implicitTimeIntegration = true;
    }
    if(solverParameters[i]->isSublist("ImplicitDiffusion")){
      implicitTimeIntegration = true;
    }
    if(solverParameters[i]->isParameter("Peridigm Preconditioner")){
      std::string peridigmPreconditionerType = solverParameters[i]->get<string>("Peridigm Preconditioner");
      if(peridigmPreconditionerType == "Full Tangent")
        userSpecifiedFullTangent = true;
      // Note:  Currently, Peridigm must have some sort of tangent to avoid null pointer errors (\todo:  Fix this!)
      //        For the time being, if the users requests NOX with no precondioner, go ahead and allocate the 3x3
      if(peridigmPreconditionerType == "Block 3x3" || peridigmPreconditionerType == "None")
        userSpecifiedBlockDiagonalTangent = true;
    }
  }
  bool allocateTangent(false), allocateBlockDiagonalTangent(false);
  if(userSpecifiedFullTangent)
    allocateTangent = true;
  if(userSpecifiedBlockDiagonalTangent)
    allocateBlockDiagonalTangent = true;
  if(implicitTimeIntegration && (!userSpecifiedFullTangent && !userSpecifiedBlockDiagonalTangent))
    allocateTangent = true;
  if(peridigmParams->isParameter("Optimization Based Coupling"))
    allocateTangent = true;

  // For the sake of the bond-associated hypoelastic model
  if(peridigmParams->isParameter("Enable Bond-Associated Hypoelastic Model"))
    analysisHasBondAssociatedHypoelasticModel = peridigmParams->get<bool>("Enable Bond-Associated Hypoelastic Model") == true;

  // If a discretization was passed into the constructor, use it.  This is done for code coupling with Albany.
  // If not, create one based on the Discretization ParameterList in the input deck.
  Teuchos::RCP<Discretization> peridigmDiscretization = inputPeridigmDiscretization;
  if(peridigmDiscretization.is_null()){
    DiscretizationFactory discFactory(discParams);
    peridigmDiscretization = discFactory.create(peridigmComm);
  }
  initializeDiscretization(peridigmDiscretization);

  // Instantiate and initialize the boundary and initial condition manager
  Teuchos::RCP<Teuchos::ParameterList> bcParams =
    Teuchos::rcpFromRef( peridigmParams->sublist("Boundary Conditions") );

  // Set a flag for creation of the RANK_DEFICIENT_NODES node set if the simulation
  // uses implicit time integration and has bond failure
  bool hasDamage(false);
  for(Teuchos::ParameterList::ConstIterator it = peridigmParams->sublist("Blocks").begin() ; it != peridigmParams->sublist("Blocks").end() ; it++){
    if(blockParams.sublist(it->first).isParameter("Damage Model"))
      hasDamage = true;
  }
  // Note that allocateTangent is true only iff it's an implicit solve
  if(allocateTangent && hasDamage){
    if(!bcParams->isParameter("Create Node Set For Rank Deficient Nodes"))
      bcParams->set<bool>("Create Node Set For Rank Deficient Nodes", true);
  }

  boundaryAndInitialConditionManager =
    Teuchos::RCP<BoundaryAndInitialConditionManager>(new BoundaryAndInitialConditionManager(*bcParams, this));

  boundaryAndInitialConditionManager->initialize(peridigmDiscretization);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  elementIdFieldId                   = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Element_Id");
  blockIdFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Block_Id");
  horizonFieldId                     = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Horizon");
  volumeFieldId                      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Volume");

  if(analysisHasMultiphysics){
    fluidPressureYFieldId            = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Fluid_Pressure_Y");
    fluidPressureUFieldId            = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Fluid_Pressure_U");
    fluidPressureVFieldId            = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Fluid_Pressure_V");
    fluidFlowDensityFieldId          = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flux_Density");
  }

  modelCoordinatesFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  coordinatesFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Coordinates");
  displacementFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Displacement");
  velocityFieldId                    = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Velocity");
  accelerationFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Acceleration");
  temperatureFieldId                 = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature");
  concentrationFieldId               = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Concentration");
  deltaTemperatureFieldId            = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Temperature_Change");
  fluxDivergenceFieldId              = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Flux_Divergence");
  concentrationFluxDivergenceFieldId = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Concentration_Flux_Divergence");
  forceDensityFieldId                = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Force_Density");
  contactForceDensityFieldId         = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "Contact_Force_Density");
  externalForceDensityFieldId        = fieldManager.getFieldId(PeridigmField::NODE,    PeridigmField::VECTOR, PeridigmField::TWO_STEP, "External_Force_Density");
  damageFieldId                      = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  jacobianDeterminantFieldId         = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Jacobian_Determinant");
  weightedVolumeFieldId              = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Undamaged_Weighted_Volume");
  velocityGradientXFieldId           = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Velocity_Gradient_X");
  velocityGradientYFieldId           = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Velocity_Gradient_Y");
  velocityGradientZFieldId           = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Velocity_Gradient_Z");

  // Create field ids that may be required for output
  fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Proc_Num");

  // Instantiate the contact manager
  Teuchos::ParameterList contactParams;
  if(peridigmParams->isSublist("Contact")){
    analysisHasContact = true;
    contactParams = peridigmParams->sublist("Contact");
    checkContactSearchRadius(contactParams,peridigmDiscretization);
    contactManager =
      Teuchos::RCP<ContactManager>(new ContactManager(contactParams, peridigmDiscretization, peridigmParams));
    contactManager->initialize(oneDimensionalMap,
                               threeDimensionalMap,
                               oneDimensionalOverlapMap,
                               bondMap,
                               globalNeighborhoodData,
                               peridigmDiscretization->getBlockID());
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

    ContactModelFactory contactModelFactory;

    contactBlocks = contactManager->getContactBlocks();

    double currentValue = 0.0;
    double previousValue = 0.0;
    double timeCurrent = 0.0;
    double timePrevious = 0.0;

    for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
    	contactModel = contactBlockIt->getContactModel();
        if(contactModel->Name() == "Time-Dependent Short-Range Force"){
            New_contactModel = Teuchos::rcp_const_cast<PeridigmNS::ContactModel> (contactModel);
            New_contactModel->evaluateParserFriction(currentValue, previousValue, timeCurrent, timePrevious);
        }
    }
  }

  // Instantiate the data loader, if requested
  if(peridigmParams->isSublist("Data Loader")){
    analysisHasDataLoader = true;
    dataLoader = Teuchos::RCP<DataLoader>(new DataLoader(peridigmParams->sublist("Data Loader"),
                                                         oneDimensionalMap));
  }

  // Instantiate the blocks
  initializeBlocks(peridigmDiscretization);

  // Determine a default finite-difference probe length
  double minElementRadius = peridigmDiscretization->getMinElementRadius();
  double defaultFiniteDifferenceProbeLength = 1.0e-6*minElementRadius;

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
    bool constantHorizon = horizonManager.blockHasConstantHorizon(blockName);
    double blockHorizon(0.0);
    if(constantHorizon)
      blockHorizon = horizonManager.getBlockConstantHorizonValue(blockName);

    // Set the material model
    string materialName = blockIt->getMaterialName();

    // Generate a relevant error message and then check whether is should be heard
    // Material names tagged with "MP" somewhere in their name are considered by
    // this code to be multiphysics compatible.
    multiphysError = "\n**** Error, for material block ";
    multiphysError += blockName;
    multiphysError += ", material ";
    multiphysError += materialName;
    multiphysError += ", is not multiphysics compatible.\n";
    //The following: If we tried to enable multiphysics, but aren't using the right material model in each material block, raise an exception.
    TEUCHOS_TEST_FOR_EXCEPT_MSG((analysisHasMultiphysics && (materialName.find("Multiphysics") == std::string::npos)), "\n**** Error, material model is not multiphysics compatible.\n");
    //The following: If we have not tried to enable multiphysics, yet are attempting to use a multiphysics material model, raise an exception.
    TEUCHOS_TEST_FOR_EXCEPT_MSG((!analysisHasMultiphysics && (materialName.find("Multiphysics") != std::string::npos)), "\n**** Error, multiphysics must be enabled at the top level of the input deck.\n");

    Teuchos::ParameterList matParams = materialParams.sublist(materialName);

    // Is the material name that of one designed for multiphysics when multiphysics is enabled?

    // If the horizon is a constant value, assign it to the material model
    // Make sure the user did not try to set the horizon in the material block
    TEUCHOS_TEST_FOR_EXCEPT_MSG(matParams.isParameter("Horizon") , "\n**** Error, Horizon is an invalid material parameter.\n");
    if(constantHorizon)
      matParams.set("Horizon", blockHorizon);

    // Assign the finite difference probe length
    if(!matParams.isParameter("Finite Difference Probe Length"))
      matParams.set("Finite Difference Probe Length", defaultFiniteDifferenceProbeLength);

    // Instantiate the material model for this block
    Teuchos::RCP<PeridigmNS::Material> materialModel = materialFactory.create(matParams);
    materialModel->setBCManager(boundaryAndInitialConditionManager);
    blockIt->setMaterialModel(materialModel);


    // Set the damage model (if any)
    double currentValue = 0.0;
    double previousValue = 0.0;
    double timeCurrent = 0.0;
    double timePrevious = 0.0;
    string damageModelName = blockIt->getDamageModelName();
    if(damageModelName != "None"){
      Teuchos::ParameterList damageParams = damageModelParams.sublist(damageModelName, true);
      Teuchos::RCP<PeridigmNS::DamageModel> damageModel = damageModelFactory.create(damageParams);
      blockIt->setDamageModel(damageModel);
      if(damageModel->Name() =="Interface Aware"){
        Teuchos::RCP< PeridigmNS::InterfaceAwareDamageModel > IADamageModel = Teuchos::rcp_dynamic_cast< PeridigmNS::InterfaceAwareDamageModel >(damageModel);
        IADamageModel->setBCManager(boundaryAndInitialConditionManager);
      }
      else if(damageModel->Name() =="Time Dependent Critical Stretch"){
        CSDamageModel = Teuchos::rcp_dynamic_cast< PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel >(damageModel);
        CSDamageModel->evaluateParserDmg(currentValue, previousValue, timeCurrent, timePrevious);
      }
    }
  }

  // Instantiate compute manager
  instantiateComputeManager(peridigmDiscretization);

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
  if(analysisHasMultiphysics) {
    auxiliaryFieldIds.push_back(fluidPressureYFieldId);
    auxiliaryFieldIds.push_back(fluidPressureUFieldId);
		auxiliaryFieldIds.push_back(fluidPressureVFieldId);
  }
  if(analysisHasBondAssociatedHypoelasticModel) {
    auxiliaryFieldIds.push_back(damageFieldId);
    auxiliaryFieldIds.push_back(jacobianDeterminantFieldId);
    auxiliaryFieldIds.push_back(weightedVolumeFieldId);
    auxiliaryFieldIds.push_back(velocityGradientXFieldId);
    auxiliaryFieldIds.push_back(velocityGradientYFieldId);
    auxiliaryFieldIds.push_back(velocityGradientZFieldId);
  }
  if(computeIntersections){
    int tempFieldId;
    auxiliaryFieldIds.push_back(blockIdFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Neighbor_Volume");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Neighbor_Centroid_X");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Neighbor_Centroid_Y");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Neighbor_Centroid_Z");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Self_Volume");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Self_Centroid_X");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Self_Centroid_Y");
    auxiliaryFieldIds.push_back(tempFieldId);
    tempFieldId = fieldManager.getFieldId(PeridigmField::BOND, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Self_Centroid_Z");
    auxiliaryFieldIds.push_back(tempFieldId);
  }

  // Add fields from compute classes to auxiliary field vector
  vector<int> computeManagerFieldIds = computeManager->FieldIds();
  auxiliaryFieldIds.insert(auxiliaryFieldIds.end(), computeManagerFieldIds.begin(), computeManagerFieldIds.end());

  // Add fields from data loader, if any
  if(analysisHasDataLoader){
    vector<int> dataLoaderFieldIds = dataLoader->getFieldIds();
    auxiliaryFieldIds.insert(auxiliaryFieldIds.end(), dataLoaderFieldIds.begin(), dataLoaderFieldIds.end());
  }

  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->setAuxiliaryFieldIds(auxiliaryFieldIds);

  // Initialize the blocks (creates maps, neighborhoods, DataManager)
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
    blockIt->initialize(peridigmDiscretization->getGlobalOwnedMap(1),
                        peridigmDiscretization->getGlobalOverlapMap(1),
                        peridigmDiscretization->getGlobalOwnedMap(3),
                        peridigmDiscretization->getGlobalOverlapMap(3),
                        peridigmDiscretization->getGlobalBondMap(),
                        blockIDs,
                        globalNeighborhoodData);

  // Create a temporary vector for storing the global element ids
  Teuchos::RCP<Epetra_Vector> elementIds = Teuchos::rcp(new Epetra_Vector(*(peridigmDiscretization->getCellVolume())));
  for(int i=0 ; i<elementIds->MyLength() ; ++i) {
    (*elementIds)[i] = elementIds->Map().GID(i);
  }

  // Load initial data into the blocks
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(peridigmDiscretization->getBlockID(),    blockIdFieldId,          PeridigmField::STEP_NONE, Insert);
    blockIt->importData(peridigmDiscretization->getHorizon(),    horizonFieldId,          PeridigmField::STEP_NONE, Insert);
    blockIt->importData(peridigmDiscretization->getCellVolume(), volumeFieldId,           PeridigmField::STEP_NONE, Insert);
    blockIt->importData(peridigmDiscretization->getInitialX(),   modelCoordinatesFieldId, PeridigmField::STEP_NONE, Insert);
    blockIt->importData(peridigmDiscretization->getInitialX(),   coordinatesFieldId,      PeridigmField::STEP_N,    Insert);
    blockIt->importData(peridigmDiscretization->getInitialX(),   coordinatesFieldId,      PeridigmField::STEP_NP1,  Insert);
    blockIt->importData(elementIds,                                 elementIdFieldId,        PeridigmField::STEP_NONE, Insert);

		if(analysisHasMultiphysics){
			scalarScratch->PutScalar(0.0);
			blockIt->importData(scalarScratch, fluidPressureYFieldId, PeridigmField::STEP_N, Insert);
			blockIt->importData(scalarScratch, fluidPressureYFieldId, PeridigmField::STEP_NP1, Insert);
		}
  }


  // Compute element-horizon intersections
#ifdef PERIDIGM_PV
  if(computeIntersections){

    // Grab parameters from discretization section of input deck
    PartialVolumeScheme partialVolumeScheme = PV;
    string partialVolumeSchemeString = "PV (default)";
    if(discParams->isParameter("Element-Horizon Intersection Partial Volume Scheme")){
      partialVolumeSchemeString = discParams->get<string>("Element-Horizon Intersection Partial Volume Scheme");
      if(partialVolumeSchemeString == "FV")
        partialVolumeScheme = FV;
      else if(partialVolumeSchemeString == "PDLAMMPS")
        partialVolumeScheme = PDLAMMPS;
      else if(partialVolumeSchemeString == "HHB")
        partialVolumeScheme = HHB;
      else if(partialVolumeSchemeString == "PV")
        partialVolumeScheme = PV;
      else
        TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Invalid option for \"Element-Horizon Intersection Partial Volume Scheme\".");
    }
    int computeIntersectionsNumRecursion = 2;
    if(discParams->isParameter("Element-Horizon Intersection Recursion Level"))
      computeIntersectionsNumRecursion = discParams->get<int>("Element-Horizon Intersection Recursion Level");
    int computeIntersectionsNumSamples = 2;
    if(discParams->isParameter("Element-Horizon Intersection Number Of Samples"))
      computeIntersectionsNumSamples = discParams->get<int>("Element-Horizon Intersection Number Of Samples");
    double computeIntersectionsCharacteristicElementLength = 0.0;
    if(discParams->isParameter("Element-Horizon Intersection Characteristic Element Length"))
      computeIntersectionsCharacteristicElementLength = discParams->get<double>("Element-Horizon Intersection Characteristic Element Length");
    bool useLookupTable = true;
    if(discParams->isParameter("Element-Horizon Intersection Use Lookup Table"))
      useLookupTable = discParams->get<bool>("Element-Horizon Intersection Use Lookup Table");

    PeridigmNS::Timer::self().startTimer("Element-Horizon Intersections");
    if(peridigmComm->MyPID() == 0){
      cout << "Computing element-horizon intersections, scheme = " << partialVolumeSchemeString << endl;
      cout.flush();
    }
    computePartialVolume(blocks,
                         peridigmDiscretization,
                         computeIntersectionsNumRecursion,
                         computeIntersectionsNumSamples,
                         partialVolumeScheme,
                         computeIntersectionsCharacteristicElementLength,
                         useLookupTable);
    if(peridigmComm->MyPID() == 0){
      cout << "\n  Intersection calculations complete.\n" << endl;
      cout.flush();
    }
    PeridigmNS::Timer::self().stopTimer("Element-Horizon Intersections");
  }
#endif

  // Store the locations of the original Exodus nodes for each element.
  // This is only done if the fields "Exodus_Node_1", "Exodus_Node_2", etc., have been registered.
  // The point of storing the node positions as element variables is to make them available in
  // compute classes or for output (which enables various types of post-processing).
  if( fieldManager.hasField("Exodus_Node_1") ){

    int m_exodusNode1FieldId = fieldManager.getFieldId("Exodus_Node_1");
    int m_exodusNode2FieldId = fieldManager.getFieldId("Exodus_Node_2");
    int m_exodusNode3FieldId = fieldManager.getFieldId("Exodus_Node_3");
    int m_exodusNode4FieldId = fieldManager.getFieldId("Exodus_Node_4");
    int m_exodusNode5FieldId = fieldManager.getFieldId("Exodus_Node_5");
    int m_exodusNode6FieldId = fieldManager.getFieldId("Exodus_Node_6");
    int m_exodusNode7FieldId = fieldManager.getFieldId("Exodus_Node_7");
    int m_exodusNode8FieldId = fieldManager.getFieldId("Exodus_Node_8");

    int globalId;
    unsigned int numExodusNodes;
    vector<double> nodePositions;
    Teuchos::RCP<const Epetra_BlockMap> blockScalarPointMap;
    double *node1, *node2, *node3, *node4, *node5, *node6, *node7, *node8;
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++) {
      blockScalarPointMap = blockIt->getOwnedScalarPointMap();
      blockIt->getData(m_exodusNode1FieldId, PeridigmField::STEP_NONE)->ExtractView(&node1);
      blockIt->getData(m_exodusNode2FieldId, PeridigmField::STEP_NONE)->ExtractView(&node2);
      blockIt->getData(m_exodusNode3FieldId, PeridigmField::STEP_NONE)->ExtractView(&node3);
      blockIt->getData(m_exodusNode4FieldId, PeridigmField::STEP_NONE)->ExtractView(&node4);
      blockIt->getData(m_exodusNode5FieldId, PeridigmField::STEP_NONE)->ExtractView(&node5);
      blockIt->getData(m_exodusNode6FieldId, PeridigmField::STEP_NONE)->ExtractView(&node6);
      blockIt->getData(m_exodusNode7FieldId, PeridigmField::STEP_NONE)->ExtractView(&node7);
      blockIt->getData(m_exodusNode8FieldId, PeridigmField::STEP_NONE)->ExtractView(&node8);

      for(int i=0 ; i<blockScalarPointMap->NumMyElements() ; ++i){
        globalId = blockScalarPointMap->GID(i);
        peridigmDiscretization->getExodusMeshNodePositions(globalId, nodePositions);
        numExodusNodes = nodePositions.size()/3;
        if(numExodusNodes >= 1){
          for(int j=0 ; j<3 ; ++j)
            node1[3*i+j] = nodePositions[j];
        }
        if(numExodusNodes >= 2){
          for(int j=0 ; j<3 ; ++j)
            node2[3*i+j] = nodePositions[3+j];
        }
        if(numExodusNodes >= 3){
          for(int j=0 ; j<3 ; ++j)
            node3[3*i+j] = nodePositions[6+j];
        }
        if(numExodusNodes >= 4){
          for(int j=0 ; j<3 ; ++j)
            node4[3*i+j] = nodePositions[9+j];
        }
        if(numExodusNodes >= 5){
          for(int j=0 ; j<3 ; ++j)
            node5[3*i+j] = nodePositions[12+j];
        }
        if(numExodusNodes >= 6){
          for(int j=0 ; j<3 ; ++j)
            node6[3*i+j] = nodePositions[15+j];
        }
        if(numExodusNodes >= 7){
          for(int j=0 ; j<3 ; ++j)
            node7[3*i+j] = nodePositions[18+j];
        }
        if(numExodusNodes >= 8){
          for(int j=0 ; j<3 ; ++j)
            node8[3*i+j] = nodePositions[21+j];
        }
      }
    }
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

  	if(analysisHasMultiphysics){
    	double blockFluidDensity = blockIt->getMaterialModel()->lookupMaterialProperty("Fluid density");
    	double blockFluidCompressibility = blockIt->getMaterialModel()->lookupMaterialProperty("Fluid compressibility");
			for(int i=0 ; i<OwnedScalarPointMap->NumMyElements() ; ++i){
				int globalID = OwnedScalarPointMap->GID(i);
				int mothershipLocalID = oneDimensionalMap->LID(globalID);
				(*fluidDensity)[mothershipLocalID] = blockFluidDensity;
				(*fluidCompressibility)[mothershipLocalID] = blockFluidCompressibility;
			}
		}
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

  // Manage data synchronization across block boundaries and MPI partitions, as needed by material models
  PeridigmNS::DataManagerSynchronizer& dataManagerSynchronizer = PeridigmNS::DataManagerSynchronizer::self();
  dataManagerSynchronizer.initialize(oneDimensionalMap, threeDimensionalMap);
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++) {
    std::vector<int> fieldIdsForSynchronizationAfterInitialize = blockIt->getMaterialModel()->FieldIdsForSynchronizationAfterInitialize();
    dataManagerSynchronizer.setFieldIdsToSynchronizeAfterInitialize(fieldIdsForSynchronizationAfterInitialize);
    std::vector<int> fieldIdsForSynchronizationAfterPrecompute = blockIt->getMaterialModel()->FieldIdsForSynchronizationAfterPrecompute();
    dataManagerSynchronizer.setFieldIdsToSynchronizeAfterPrecompute(fieldIdsForSynchronizationAfterPrecompute);
  }
  dataManagerSynchronizer.checkFieldValidity(blocks);
  dataManagerSynchronizer.synchronizeDataAfterInitialize(blocks);

  // Create the model evaluator
  modelEvaluator = Teuchos::rcp(new PeridigmNS::ModelEvaluator());

  // Initialize output manager
  initializeOutputManager();

  // Call rebalance function if analysis has contact
  // this is required to set up proper contact neighbor list
  if(analysisHasContact)
    contactManager->rebalance(0);

  // Create service manager
  serviceManager = Teuchos::rcp(new PeridigmNS::ServiceManager());
  serviceManager->requestService(computeManager->Services());

  jacobianType = PeridigmNS::Material::UNDEFINED;

  // Perform requested services
  if (serviceManager->isRequested(PeridigmNS::PeridigmService::ALLOCATE_TANGENT) || allocateTangent) {
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
      if(numMultiphysDoFs > 0)
        cout << "  of those rows, " << numMultiphysDoFs << " are interspersed multiphysics terms." << endl;
      cout << "  number of nonzeros = " << tangent->NumGlobalNonzeros() << "\n" << endl;
    }
    jacobianType = PeridigmNS::Material::FULL_MATRIX;
  }

  // Check if request for allocation of block diagonal tangent stiffness matrix
  if (serviceManager->isRequested(PeridigmNS::PeridigmService::ALLOCATE_BLOCK_DIAGONAL_TANGENT) || allocateBlockDiagonalTangent) {
    // Allocate memory for non-zeros in global tangent and lock in the structure
    if(peridigmComm->MyPID() == 0 && !allocateTangent){
      cout << "Allocating global block diagonal tangent matrix...";
      cout.flush();
    }
    PeridigmNS::Timer::self().startTimer("Allocate Global Block Diagonal Tangent");
    allocateBlockDiagonalJacobian();
    // If both the full tangent and the block diagonal tangent are flagged for allocation,
    // only the full tangent is allocated and the block diagonal just points to the full tangent.
    // If only the block diagonal is allocated, the the pointers for the tangent should
    // be set to the block diagonal tangent so that the block diagonal tangent gets filled.
    tangentMap = blockDiagonalTangentMap;
    tangent = blockDiagonalTangent;
    PeridigmNS::Timer::self().stopTimer("Allocate Global Block Diagonal Tangent");
    if(peridigmComm->MyPID() == 0 && !allocateTangent){
      cout << "\n  number of rows = " << blockDiagonalTangent->NumGlobalRows() << endl;
      cout << "  number of nonzeros = " << blockDiagonalTangent->NumGlobalNonzeros() << "\n" << endl;
    }
    if(jacobianType == PeridigmNS::Material::UNDEFINED)
      jacobianType = PeridigmNS::Material::BLOCK_DIAGONAL;
  }
  //Initialize restart if requested in the input file
  if(peridigmParams->isParameter("Restart")){
	 InitializeRestart();
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

void PeridigmNS::Peridigm::initializeDiscretization(Teuchos::RCP<Discretization> peridigmDisc) {

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();

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

  // threeDimensionalOverlapMap
  // used for NeumannBC
  threeDimensionalOverlapMap = peridigmDisc->getGlobalOverlapMap(3);

  // unknonwnsMap
  // used for time integrators / solvers
  unknownsMap = Teuchos::rcp(new Epetra_BlockMap(-1,
                                                 oneDimensionalMap->NumMyElements(),
                                                 oneDimensionalMap->MyGlobalElements(),
                                                 dofManager.totalNumberOfDegreesOfFreedom(),
                                                 0,
                                                 oneDimensionalMap->Comm()));

  // bondConstitutiveDataMap
  // a non-overlapping map used for storing constitutive data on bonds
  bondMap = peridigmDisc->getGlobalBondMap();

  // Create mothership vectors
  int numOneDimensionalMothershipVectors = 9;
  bool initializeToZero = true;
  if(analysisHasMultiphysics)
    numOneDimensionalMothershipVectors += 7;
  if(analysisHasBondAssociatedHypoelasticModel)
    numOneDimensionalMothershipVectors += 3;

  oneDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*oneDimensionalMap, numOneDimensionalMothershipVectors, initializeToZero));
  blockIDs = Teuchos::rcp((*oneDimensionalMothership)(0), false);                    // block ID
  horizon = Teuchos::rcp((*oneDimensionalMothership)(1), false);                     // horizon for each point
  volume = Teuchos::rcp((*oneDimensionalMothership)(2), false);                      // cell volume
  density = Teuchos::rcp((*oneDimensionalMothership)(3), false);                     // solid density
  temperature = Teuchos::rcp((*oneDimensionalMothership)(4), false);                 // temperature
  deltaTemperature = Teuchos::rcp((*oneDimensionalMothership)(5), false);            // change in temperature
  fluxDivergence = Teuchos::rcp((*oneDimensionalMothership)(6), false);              // divergence of the flux (e.g., heat flux)
  concentrationFluxDivergence = Teuchos::rcp((*oneDimensionalMothership)(7), false); // divergence of the flux of chemical concentration
  scalarScratch = Teuchos::rcp((*oneDimensionalMothership)(8), false);               // scratch vector corresponding to oneDimensionalMap
  if (analysisHasMultiphysics) {
    fluidPressureU = Teuchos::rcp((*oneDimensionalMothership)(9), false);        // fluid pressure displacement
    fluidPressureY = Teuchos::rcp((*oneDimensionalMothership)(10), false);       // fluid pressure current coordinates at anode
    fluidPressureV = Teuchos::rcp((*oneDimensionalMothership)(11), false);       // fluid pressure first time derv at a node
    fluidFlow = Teuchos::rcp((*oneDimensionalMothership)(12), false);            // flux through a node
    fluidPressureDeltaU = Teuchos::rcp((*oneDimensionalMothership)(13), false);  // fluid pressure displacement analogue increment
    fluidDensity = Teuchos::rcp((*oneDimensionalMothership)(14), false); 		     // fluid density at a node
    fluidCompressibility = Teuchos::rcp((*oneDimensionalMothership)(15), false); // fluid compressibility at a node
    if(analysisHasBondAssociatedHypoelasticModel){
      damage = Teuchos::rcp((*oneDimensionalMothership)(16), false);              // damage
      jacobianDeterminant = Teuchos::rcp((*oneDimensionalMothership)(17), false); // jacobian determinant (J)
      weightedVolume = Teuchos::rcp((*oneDimensionalMothership)(18), false);      // weighted volume
      damage->PutScalar(0.0);
      jacobianDeterminant->PutScalar(1.0);
    }
  }
  if(analysisHasBondAssociatedHypoelasticModel){
    damage = Teuchos::rcp((*oneDimensionalMothership)(9), false);              // damage
    jacobianDeterminant = Teuchos::rcp((*oneDimensionalMothership)(10), false); // jacobian determinant (J)
    weightedVolume = Teuchos::rcp((*oneDimensionalMothership)(11), false);      // weighted volume
    damage->PutScalar(0.0);
    jacobianDeterminant->PutScalar(1.0);
  }

  int numThreeDimensionalMothershipVectors = 10;
  if(analysisHasBondAssociatedHypoelasticModel)
    numThreeDimensionalMothershipVectors += 3;

  threeDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*threeDimensionalMap, numThreeDimensionalMothershipVectors, initializeToZero));
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
  if(analysisHasBondAssociatedHypoelasticModel){
    velocityGradientX = Teuchos::rcp((*threeDimensionalMothership)(10), false);  // velocity gradient XX XY XZ (L)
    velocityGradientY = Teuchos::rcp((*threeDimensionalMothership)(11), false); // velocity gradient YX YY YZ (L)
    velocityGradientZ = Teuchos::rcp((*threeDimensionalMothership)(12), false); // velocity gradient ZX ZY ZZ (L)
  }

  unknownsMothership = Teuchos::rcp(new Epetra_MultiVector(*unknownsMap, 5, initializeToZero));
  unknownsU = Teuchos::rcp((*unknownsMothership)(0), false);             // abstract displacement
  unknownsY = Teuchos::rcp((*unknownsMothership)(1), false);             // abstract current positions
  unknownsV = Teuchos::rcp((*unknownsMothership)(2), false);             // abstract velocities
  unknownsForce = Teuchos::rcp((*unknownsMothership)(3), false);         // abstract force
  unknownsDeltaU = Teuchos::rcp((*unknownsMothership)(4), false);        // abstract increment in displacement (used only for implicit time integration)

  // Set the block IDs
  double* bID;
  peridigmDisc->getBlockID()->ExtractView(&bID);
  double* blockIDsPtr;
  blockIDs->ExtractView(&blockIDsPtr);
  blas.COPY(blockIDs->MyLength(), bID, blockIDsPtr);

  // Set the horizon values
  double* discHorizonPtr;
  peridigmDisc->getHorizon()->ExtractView(&discHorizonPtr);
  double* horizonPtr;
  horizon->ExtractView(&horizonPtr);
  blas.COPY(horizon->MyLength(), discHorizonPtr, horizonPtr);

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

  if(peridigmDisc->InterfacesAreConstructed()){
    interfaceData = peridigmDisc->getInterfaceData();
    constructInterfaces = true;
  }
}

void PeridigmNS::Peridigm::initializeWorkset() {
  workset = Teuchos::rcp(new Workset);
  workset->timeStep = 0.0;
  workset->blocks = blocks;
  if(!contactManager.is_null())
    workset->contactManager = contactManager;
  workset->jacobianType = Teuchos::rcpFromRef(jacobianType);
  workset->jacobian = overlapJacobian;
}

std::string getCmdOutput(const std::string& mStr)
{
    std::string result, file;
    FILE* pipe;
    char buffer[256];
    pipe=popen(mStr.c_str(), "r");
    while(fgets(buffer, sizeof(buffer), pipe) != NULL)
    {
        file = buffer;
        result += file.substr(0, file.size() - 1);
    }

    pclose(pipe);
    return result;
}
std::string firstNumbersSring(std::string const & str)
{
  std::size_t const n = str.find_first_of("0123456789");
  if (n != std::string::npos)
  {
    std::size_t const m = str.find_first_not_of("0123456789", n);
    return str.substr(n, m != std::string::npos ? m-n : m);
  }
  return std::string();
}
void PeridigmNS::Peridigm::InitializeRestart() {
	std::string str;
	struct stat sb;
	char const * restart_directory_namePtr;
	if (stat("restart-000001", &sb) == 0 && S_ISDIR(sb.st_mode)){
	    str=getCmdOutput("ls -td -- ./restart*/ | head -n1 | cut -d'/' -f2");
	    if (str != ""){
	        if(peridigmComm->MyPID() == 0){
	        	cout <<"Restart folder exists, will attempt to read the restart files. \n"<< endl;
	        	cout.flush();
	        }
	    	std::vector<char> writable(str.begin(), str.end());
	    	writable.push_back('\0');
	    	restart_directory_namePtr=&*writable.begin();
	    	setRestartNames(restart_directory_namePtr);
	    	readRestart();
	        if(peridigmComm->MyPID() == 0){
			    	cout <<"Restart is initialized." << endl;
		        	cout.flush();
		    }
	    }else{
	    	if(peridigmComm->MyPID() == 0){
	    		TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "Error: Initial restart folder exists, but it is not suitable for a restart. \n");
	    		MPI_Finalize();
	    		exit(0);
	    	}
	    }
	}else{
	    if(peridigmComm->MyPID() == 0){
	    	cout <<"Initial restart folder does not exist." << endl;
			cout.flush();
		}
		restart_directory_namePtr ="restart-000000";
		setRestartNames(restart_directory_namePtr);
	}
}

void PeridigmNS::Peridigm::setRestartNames(	char const * restart_directory_namePtr) {
char pathname[100];
//path to current restart folder
restartFiles["path"] = restart_directory_namePtr;
//Current time restart file
sprintf(pathname,"%s/currentTime.txt",restart_directory_namePtr);
restartFiles["currentTime"] = pathname;
//blockIDs restart file
sprintf(pathname,"%s/blockIDs.mat",restart_directory_namePtr);
restartFiles["blockIDs"] = pathname;
//horizon restart file
sprintf(pathname,"%s/horizon.mat",restart_directory_namePtr);
restartFiles["horizon"] = pathname;

//volume restart file
sprintf(pathname,"%s/volume.mat",restart_directory_namePtr);
restartFiles["volume"] = pathname;

//density restart file
sprintf(pathname,"%s/density.mat",restart_directory_namePtr);
restartFiles["density"] = pathname;

//deltaTemperature restart file
sprintf(pathname,"%s/deltaTemperature.mat",restart_directory_namePtr);
restartFiles["deltaTemperature"] = pathname;

//x restart file
sprintf(pathname,"%s/x.mat",restart_directory_namePtr);
restartFiles["x"] = pathname;

//u restart file
sprintf(pathname,"%s/u.mat",restart_directory_namePtr);
restartFiles["u"] = pathname;

//y restart file
sprintf(pathname,"%s/y.mat",restart_directory_namePtr);
restartFiles["y"] = pathname;

//v restart file
sprintf(pathname,"%s/v.mat",restart_directory_namePtr);
restartFiles["v"] = pathname;

//a restart file
sprintf(pathname,"%s/a.mat",restart_directory_namePtr);
restartFiles["a"] = pathname;

//force restart file
sprintf(pathname,"%s/force.mat",restart_directory_namePtr);
restartFiles["force"] = pathname;

//contactForce restart file
sprintf(pathname,"%s/contactForce.mat",restart_directory_namePtr);
restartFiles["contactForce"] = pathname;

//externalForce restart file
sprintf(pathname,"%s/externalForce.mat",restart_directory_namePtr);
restartFiles["externalForce"] = pathname;

//deltaU restart file
sprintf(pathname,"%s/deltaU.mat",restart_directory_namePtr);
restartFiles["deltaU"] = pathname;

//scratch restart file
sprintf(pathname,"%s/scratch.mat",restart_directory_namePtr);
restartFiles["scratch"] = pathname;
}
void PeridigmNS::Peridigm::instantiateComputeManager(Teuchos::RCP<Discretization> peridigmDiscretization) {

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
  Teuchos::RCP<Discretization> *tmp5 = &( peridigmDiscretization );
  Teuchos::RCP<int> *tmp6 = &( nonlinearSolverIterations );
  computeClassGlobalData->set("tangent",tmp1);
  computeClassGlobalData->set("blockDiagonalTangent",tmp2);
  computeClassGlobalData->set("overlapJacobian",tmp3);
  computeClassGlobalData->set("blockDiagonalTangentMap",tmp4);
  computeClassGlobalData->set("discretization",tmp5);
  computeClassGlobalData->set("nonlinearSolverIterations",tmp6);

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

void PeridigmNS::Peridigm::execute(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  TEUCHOS_TEST_FOR_EXCEPT_MSG(solverParams.is_null(), "Error in Peridigm::execute, solverParams is null.\n");

  if(solverParams->isSublist("Verlet")){
    executeExplicit(solverParams);}
  else if(solverParams->isSublist("QuasiStatic"))
    executeQuasiStatic(solverParams);
  else if(solverParams->isSublist("NOXQuasiStatic"))
    executeNOXQuasiStatic(solverParams);
  else if(solverParams->isSublist("Implicit"))
    executeImplicit(solverParams);
  else if(solverParams->isSublist("ImplicitDiffusion"))
    executeImplicitDiffusion(solverParams);
  else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error: Unrecognized time integration scheme.\n");
  }

  PeridigmNS::Memstat * memstat = PeridigmNS::Memstat::Instance();
  const std::string statTag = "Post Execute";
  memstat->addStat(statTag);
}

void PeridigmNS::Peridigm::executeSolvers() {
  for(unsigned int i=0 ; i<solverParameters.size() ; ++i){
    execute(solverParameters[i]);
    if(peridigmParams->isParameter("Restart")){
    	writeRestart(solverParameters[i]);
    }
  }
}

void PeridigmNS::Peridigm::executeExplicit(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

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
  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal   = solverParams->get("Final Time", 1.0);
  double timeCurrent = timeInitial;
  double timePrevious = timeCurrent;
  workset->timeStep = dt;
  double dt2 = dt/2.0;
  int nsteps = static_cast<int>( floor((timeFinal-timeInitial)/dt) );

  // Check to make sure the number of time steps is sane
  if(floor((timeFinal-timeInitial)/dt) > static_cast<double>(INT_MAX)){
    if(peridigmComm->MyPID() == 0){
      cout << "WARNING:  The number of time steps exceed the maximum allowable value for an integer." << endl;
      cout << "          The number of steps will be reduced to " << INT_MAX << "." << endl;
      cout << "          Any chance you botched the units in your input deck?\n" << endl;
    }
    nsteps = INT_MAX;
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
    cout << "Total number of time steps " << nsteps << "\n" << endl;
  }

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
    blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
    if(analysisHasBondAssociatedHypoelasticModel){
      blockIt->importData(damage, damageFieldId, PeridigmField::STEP_N, Insert); // Note that damage lags one step in the model evaluation
      blockIt->importData(jacobianDeterminant, jacobianDeterminantFieldId, PeridigmField::STEP_N, Insert); // Note that J lags one step in the model evaluation
    }
  }
  if(analysisHasContact)
    contactManager->importData(volume, y, v);
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  if(analysisHasBondAssociatedHypoelasticModel){
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->computeVelocityGradient(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");
    
    // Copy data from mothership vectors to overlap vectors in data manager
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    jacobianDeterminant->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scalarScratch->PutScalar(0.0);
      blockIt->exportData(scalarScratch, jacobianDeterminantFieldId, PeridigmField::STEP_NP1, Add);
      jacobianDeterminant->Update(1.0, *scalarScratch, 1.0);
    }
    weightedVolume->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scalarScratch->PutScalar(0.0);
      blockIt->exportData(scalarScratch, weightedVolumeFieldId, PeridigmField::STEP_NONE, Add);
      weightedVolume->Update(1.0, *scalarScratch, 1.0);
    }
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(weightedVolume, weightedVolumeFieldId, PeridigmField::STEP_NONE, Insert); 
    }
    velocityGradientX->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(scratch, velocityGradientXFieldId, PeridigmField::STEP_NONE, Add);
      velocityGradientX->Update(1.0, *scratch, 1.0);
    }
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(velocityGradientX, velocityGradientXFieldId, PeridigmField::STEP_NONE, Insert); 
    }
    velocityGradientY->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(scratch, velocityGradientYFieldId, PeridigmField::STEP_NONE, Add);
      velocityGradientY->Update(1.0, *scratch, 1.0);
    }
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(velocityGradientY, velocityGradientYFieldId, PeridigmField::STEP_NONE, Insert); 
    }
    velocityGradientZ->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(scratch, velocityGradientZFieldId, PeridigmField::STEP_NONE, Add);
      velocityGradientZ->Update(1.0, *scratch, 1.0);
    }
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(velocityGradientZ, velocityGradientZFieldId, PeridigmField::STEP_NONE, Insert); 
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Compute bond-level velocity gradient
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->computeBondVelocityGradient(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");
  }

  // Load the data manager with data from disk, if requested
  if(analysisHasDataLoader){
    PeridigmNS::Timer::self().startTimer("Data Loader");
    dataLoader->loadData(timeCurrent, blocks);
    PeridigmNS::Timer::self().stopTimer("Data Loader");
  }

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
    blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
    force->Update(1.0, *scratch, 1.0);
  }
  if(analysisHasContact){
    contactManager->exportData(contactForce);
    force->Update(1.0, *contactForce, 1.0);
  }
  if(analysisHasBondAssociatedHypoelasticModel){
    damage->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scalarScratch->PutScalar(0.0);
      blockIt->exportData(scalarScratch, damageFieldId, PeridigmField::STEP_NP1, Add);
      damage->Update(1.0, *scalarScratch, 1.0);
    }
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  // Apply BC at time zero
  PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
  boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");
  PeridigmNS::Timer::self().startTimer("Apply Body Forces");
  boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

  // fill the acceleration vector
  (*a) = (*force);
  for(int i=0 ; i<a->MyLength() ; ++i){
    (*a)[i] += (*externalForce)[i];
    (*a)[i] /= (*density)[i/3];
  }
  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  synchDataManagers();
  if(analysisHasDataLoader){
    dataLoader->loadData(timeCurrent, blocks);
  }
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  int displayTrigger = nsteps/100;
  if(displayTrigger == 0)
    displayTrigger = 1;

  Teuchos::ParameterList damageModelParams;
  if(peridigmParams->isSublist("Damage Models"))
    damageModelParams = peridigmParams->sublist("Damage Models");
  DamageModelFactory damageModelFactory;

  double currentValue = 0.0;
  double previousValue = 0.0;

  for(int step=1; step<=nsteps; step++){

    timePrevious = timeCurrent;
    timeCurrent = timeInitial + (step*dt);

    // TODO this should not be here
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      string damageModelName = blockIt->getDamageModelName();
      if(damageModelName != "None"){
         Teuchos::ParameterList damageParams = damageModelParams.sublist(damageModelName, true);
         Teuchos::RCP<PeridigmNS::DamageModel> damageModel = damageModelFactory.create(damageParams);
         blockIt->setDamageModel(damageModel);
         if(damageModel->Name() == "Time Dependent Critical Stretch"){
           CSDamageModel = Teuchos::rcp_dynamic_cast< PeridigmNS::UserDefinedTimeDependentCriticalStretchDamageModel >(damageModel);
           CSDamageModel->evaluateParserDmg(currentValue, previousValue, timeCurrent, timePrevious);
         }
      }
    }

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
    PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent, timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");

    // evaluate the external (body) forces:
    PeridigmNS::Timer::self().startTimer("Apply Body Forces");
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent, timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

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
      blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
      if(analysisHasBondAssociatedHypoelasticModel){
        blockIt->importData(damage, damageFieldId, PeridigmField::STEP_N, Insert); // Note that damage lags one step in the model evaluation
        blockIt->importData(jacobianDeterminant, jacobianDeterminantFieldId, PeridigmField::STEP_N, Insert); // Note that J lags one step in the model evaluation
      }
    }
    if(analysisHasContact){
      for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
        contactModel = contactBlockIt->getContactModel();
        if(contactModel->Name() == "Time-Dependent Short-Range Force"){
          New_contactModel = Teuchos::rcp_const_cast<PeridigmNS::ContactModel> (contactModel);
          New_contactModel->evaluateParserFriction(currentValue, previousValue, timeCurrent, timePrevious);
        }
      }
      contactManager->importData(volume, y, v);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    if(analysisHasBondAssociatedHypoelasticModel){
      PeridigmNS::Timer::self().startTimer("Internal Force");
      modelEvaluator->computeVelocityGradient(workset);
      PeridigmNS::Timer::self().stopTimer("Internal Force");
      
      // Copy data from mothership vectors to overlap vectors in data manager
      PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      jacobianDeterminant->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scalarScratch->PutScalar(0.0);
        blockIt->exportData(scalarScratch, jacobianDeterminantFieldId, PeridigmField::STEP_NP1, Add);
        jacobianDeterminant->Update(1.0, *scalarScratch, 1.0);
      }
      weightedVolume->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scalarScratch->PutScalar(0.0);
        blockIt->exportData(scalarScratch, weightedVolumeFieldId, PeridigmField::STEP_NONE, Add);
        weightedVolume->Update(1.0, *scalarScratch, 1.0);
      }
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(weightedVolume, weightedVolumeFieldId, PeridigmField::STEP_NONE, Insert); 
      }
      velocityGradientX->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scratch->PutScalar(0.0);
        blockIt->exportData(scratch, velocityGradientXFieldId, PeridigmField::STEP_NONE, Add);
        velocityGradientX->Update(1.0, *scratch, 1.0);
      }
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(velocityGradientX, velocityGradientXFieldId, PeridigmField::STEP_NONE, Insert); 
      }
      velocityGradientY->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scratch->PutScalar(0.0);
        blockIt->exportData(scratch, velocityGradientYFieldId, PeridigmField::STEP_NONE, Add);
        velocityGradientY->Update(1.0, *scratch, 1.0);
      }
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(velocityGradientY, velocityGradientYFieldId, PeridigmField::STEP_NONE, Insert); 
      }
      velocityGradientZ->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scratch->PutScalar(0.0);
        blockIt->exportData(scratch, velocityGradientZFieldId, PeridigmField::STEP_NONE, Add);
        velocityGradientZ->Update(1.0, *scratch, 1.0);
      }
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(velocityGradientZ, velocityGradientZFieldId, PeridigmField::STEP_NONE, Insert); 
      }
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Compute bond-level velocity gradient
      PeridigmNS::Timer::self().startTimer("Internal Force");
      modelEvaluator->computeBondVelocityGradient(workset);
      PeridigmNS::Timer::self().stopTimer("Internal Force");
    }
    // Load the data manager with data from disk, if requested
    if(analysisHasDataLoader){
      PeridigmNS::Timer::self().startTimer("Data Loader");
      dataLoader->loadData(timeCurrent, blocks);
      PeridigmNS::Timer::self().stopTimer("Data Loader");
    }

    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");

    // Copy force from the data manager to the mothership vector
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    force->PutScalar(0.0);
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
      force->Update(1.0, *scratch, 1.0);
    }
    if(analysisHasBondAssociatedHypoelasticModel){
      damage->PutScalar(0.0);
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        scalarScratch->PutScalar(0.0);
        blockIt->exportData(scalarScratch, damageFieldId, PeridigmField::STEP_NP1, Add);
        damage->Update(1.0, *scalarScratch, 1.0);
      }
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Check for NaNs in force evaluation
    // We'd like to know now because a NaN will likely cause a difficult-to-unravel crash downstream.
    for(int i=0 ; i<force->MyLength() ; ++i)
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*force)[i]), "**** NaN returned by force evaluation.\n");

    // Check for NaNs in force evaluation
    // We'd like to know now because a NaN will likely cause a difficult-to-unravel crash downstream.
    for(int i=0 ; i<externalForce->MyLength() ; ++i)
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*externalForce)[i]), "**** NaN returned by external force evaluation.\n");

    if(analysisHasContact){
      contactManager->exportData(contactForce);
      // Check for NaNs in contact force evaluation
      for(int i=0 ; i<contactForce->MyLength() ; ++i)
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*contactForce)[i]), "**** NaN returned by contact force evaluation.\n");
      // Add contact forces to forces
      force->Update(1.0, *contactForce, 1.0);
    }

    // fill the acceleration vector
    (*a) = (*force);
    for(int i=0 ; i<a->MyLength() ; ++i){
      (*a)[i] += (*externalForce)[i];
      (*a)[i] /= (*density)[i/3];
    }

    // V^{n+1}   = V^{n+1/2} + (dt/2)*A^{n+1}
    //blas.AXPY(const int N, const double ALPHA, const double *X, double *Y, const int INCX=1, const int INCY=1) const
    blas.AXPY(length, dt2, aPtr, vPtr, 1, 1);

    PeridigmNS::Timer::self().startTimer("Output");
    synchDataManagers();
    if(analysisHasDataLoader){
      dataLoader->loadData(timeCurrent, blocks);
    }
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
  return evaluateNOX(fillType, &x, &FVec);
}

bool PeridigmNS::Peridigm::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) {
  return evaluateNOX(NOX::Epetra::Interface::Required::Jac, &x, 0);
}

bool PeridigmNS::Peridigm::computePreconditioner(const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams) {

  if(agePeridigmPreconditioner > 0 && agePeridigmPreconditioner < maxAgePeridigmPreconditioner)
    return true;
  if(agePeridigmPreconditioner >= maxAgePeridigmPreconditioner)
    agePeridigmPreconditioner = 0;
  agePeridigmPreconditioner += 1;

  // Call evaluateNOX() with the Jac flag to evaluate the tangent (or 3x3 sub-tangent)
  evaluateNOX(NOX::Epetra::Interface::Required::Jac, &x, NULL);

  // Invert the 3x3 block tangent
  PeridigmNS::Timer::self().startTimer("Invert 3x3 Block Tangent");
  TEUCHOS_TEST_FOR_EXCEPT_MSG(tangent->NumMyRows()%3 != 0, "****Error in Peridigm::computePreconditioner(), invalid number of rows.\n");
  int numEntries, err;
  double *valuesRow1, *valuesRow2, *valuesRow3;
  double matrix[9], determinant, inverse[9];
  for(int iBlock=0 ; iBlock<tangent->NumMyRows() ; iBlock+=3){
    err = tangent->ExtractMyRowView(iBlock, numEntries, valuesRow1);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::computePreconditioner(), tangent->ExtractMyRowView() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numEntries != 3, "**** PeridigmNS::Peridigm::computePreconditioner(), number of row entries not equal to three (block 3x3 matrix required).\n");
    for(int i=0 ; i<3 ; ++i)
      matrix[i] = valuesRow1[i];
    err = tangent->ExtractMyRowView(iBlock+1, numEntries, valuesRow2);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::computePreconditioner(), tangent->ExtractMyRowView() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numEntries != 3, "**** PeridigmNS::Peridigm::computePreconditioner(), number of row entries not equal to three (block 3x3 matrix required).\n");
    for(int i=0 ; i<3 ; ++i)
      matrix[3+i] = valuesRow2[i];
    err = tangent->ExtractMyRowView(iBlock+2, numEntries, valuesRow3);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::computePreconditioner(), tangent->ExtractMyRowView() returned nonzero error code.\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numEntries != 3, "**** PeridigmNS::Peridigm::computePreconditioner(), number of row entries not equal to three (block 3x3 matrix required).\n");
    for(int i=0 ; i<3 ; ++i)
      matrix[6+i] = valuesRow3[i];
    err = CORRESPONDENCE::Invert3by3Matrix(matrix, determinant, inverse);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::computePreconditioner(), Invert3by3Matrix() returned nonzero error code.\n");
    for(int i=0 ; i<3 ; ++i){
      valuesRow1[i] = inverse[i];
      valuesRow2[i] = inverse[3+i];
      valuesRow3[i] = inverse[6+i];
    }
  }
  PeridigmNS::Timer::self().stopTimer("Invert 3x3 Block Tangent");

  return true;
}

bool PeridigmNS::Peridigm::evaluateNOX(NOX::Epetra::Interface::Required::FillType flag, 
                                       const Epetra_Vector* soln,
                                       Epetra_Vector* tmp_rhs)
{
  //Determine what to fill (F or Jacobian)
  bool fillF = false;
  bool fillMatrix = false;
  
  // "flag" can be used to determine how accurate your fill of F should be 
  // depending on why we are calling evaluate (Could be using computeF to 
  // populate a Jacobian or Preconditioner).
  if (flag == NOX::Epetra::Interface::Required::Residual ||
      flag == NOX::Epetra::Interface::Required::FD_Res   ||
      flag == NOX::Epetra::Interface::Required::MF_Res   ||
      flag == NOX::Epetra::Interface::Required::MF_Jac) {
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

  // Multiphysics: copy the solution vector passed in by NOX to update the deformation 
  if(analysisHasMultiphysics){
    for(int i=0 ; i < soln->MyLength() ; i+=(3+numMultiphysDoFs)){
        for(int j=0 ; j < 3 ; ++j){
            (*y)[i*3/(3+numMultiphysDoFs) + j] = (*unknownsU)[i+j] + (*x)[i*3/(3+numMultiphysDoFs) +j] + (*soln)[i+j];
            (*v)[i*3/(3+numMultiphysDoFs) + j] = (*soln)[i+j] / (workset->timeStep);
				}
	      (*fluidPressureY)[i/(3+numMultiphysDoFs)] = (*unknownsU)[i+3] + (*soln)[i+3];
	      (*fluidPressureV)[i/(3+numMultiphysDoFs)] = (*soln)[i+3] / (workset->timeStep);
    }
		//Necessary because soln is zero at dof with kinematic bc
		v->Update(1.0, *noxVelocityAtDOFWithKinematicBC, 1.0);
		fluidPressureV->Update(1.0, *noxPressureVAtDOFWithKinematicBC, 1.0);
		for(int i=0 ; i<v->MyLength() ; i+=3){
			for(int j=0 ; j<3 ; ++j){
				(*unknownsY)[i/3*(3+numMultiphysDoFs) + j] = (*y)[i+j];
				(*unknownsV)[i/3*(3+numMultiphysDoFs) + j] = (*v)[i+j];
			}
			(*unknownsY)[i/3*(3+numMultiphysDoFs) + 3] = (*fluidPressureY)[i/3];
			(*unknownsV)[i/3*(3+numMultiphysDoFs) + 3] = (*fluidPressureV)[i/3];
		}
  }
  else{
    // copy the solution vector passed in by NOX to update the deformation 
    for(int i=0 ; i < u->MyLength() ; ++i){
        (*y)[i] = (*x)[i] + (*u)[i] + (*soln)[i];
        (*v)[i] = (*soln)[i] / (workset->timeStep);
    }
  	v->Update(1.0, *noxVelocityAtDOFWithKinematicBC, 1.0); // Necessary because soln is zero at dof with kinematic bc
  }

  // Copy data from mothership vectors to overlap vectors in data manager
  if(analysisHasMultiphysics){
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(fluidPressureU, fluidPressureUFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(fluidPressureY, fluidPressureYFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(fluidPressureV, fluidPressureVFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
  }
  else{
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
        blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
  }

  if(fillF){
    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");

    if(analysisHasMultiphysics){
    	PeridigmNS::Timer::self().startTimer("Gather/Scatter");
      unknownsForce->PutScalar(0.0);
    	for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      		scratch->PutScalar(0.0);
      		scalarScratch->PutScalar(0.0);
      		//scratchCombined->PutScalar(0.0);

					//Assign values to scratch vectors from overlap data
      		blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
      		blockIt->exportData(scalarScratch, fluidFlowDensityFieldId, PeridigmField::STEP_NP1, Add);

					//Add the values from the scratch vectors to the uncombined mothership vectors
      		force->Update(1.0, *scratch, 1.0);
					fluidFlow->Update(1.0, *scalarScratch, 1.0);

					//Copy the values from the uncombined mothership vectors into the combined mothership vector 
					for(int i=0 ; i < force->MyLength() ; i+=3){
						for(int j=0 ; j < 3 ; ++j){
							(*unknownsForce)[i/3*(3+numMultiphysDoFs) + j] = (*force)[i+j];
						}
						(*unknownsForce)[i/3*(3+numMultiphysDoFs) + 3] = (*fluidFlow)[i/3];
					}
			}
    }
    else {
	    // Copy force from the data manager to the mothership vector
	    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
	    force->PutScalar(0.0);

	    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
	      scratch->PutScalar(0.0);
	      blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
	      force->Update(1.0, *scratch, 1.0);
	    }
    }
    scratch->PutScalar(0.0);
    if(analysisHasMultiphysics){
	    scalarScratch->PutScalar(0.0);
    }

    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Create residual vector
    Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));

    // copy the internal force to the residual vector
    // note that due to restrictions on CrsMatrix, these vectors have different (but equivalent) maps
    if(not analysisHasMultiphysics){
        TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != force->MyLength(), "**** PeridigmNS::Peridigm::evaluateNOX() incompatible vector lengths!\n");
    }
    else{
        TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != unknownsForce->MyLength(), "**** PeridigmNS::Peridigm::evaluateNOX() incompatible vector lengths! (residual with unknownsForce)\n");
    }

    if(analysisHasMultiphysics){
	//Store abstract force density immediately converted to force
        for(int i=0 ; i < unknownsForce->MyLength() ; ++i)
            (*residual)[i] = (*unknownsForce)[i]*(*volume)[i/(3+numMultiphysDoFs)];
    }
    else{
        for(int i=0 ; i < force->MyLength() ; ++i)
	  (*residual)[i] = (*force)[i] + (*externalForce)[i];
	// convert force density to force
    	for(int i=0 ; i < residual->MyLength() ; ++i)
      	    (*residual)[i] *= (*volume)[i/3];
    }

    // zero out the rows corresponding to kinematic boundary conditions
    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);

    // copy back to tmp_rhs
    TEUCHOS_TEST_FOR_EXCEPT_MSG(residual->MyLength() != tmp_rhs->MyLength(), "**** PeridigmNS::Peridigm::evaluateNOX() incompatible vector lengths! (tmp_rhs with residual)\n");
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

void PeridigmNS::Peridigm::computeInternalForce()
{

   TEUCHOS_TEST_FOR_EXCEPT_MSG(analysisHasMultiphysics, "**** PeridigmNS::Peridigm::computeInternalForce() is not multiphysics compatible.\n");
  // This function is intended for use when Peridigm is called as an external library (e.g., code coupling)
  // It is assumed that the global vectors x, u, y, and v have already been set by the driver application

  // Run some checks to make sure things haven't gone haywire
  for(int i=0 ; i<u->MyLength() ; ++i){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*x)[i]), "**** NaN detetected in vector x in Peridigm::computeInternalForce().\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*u)[i]), "**** NaN detetected in vector u in Peridigm::computeInternalForce().\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*y)[i]), "**** NaN detetected in vector y in Peridigm::computeInternalForce().\n");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*v)[i]), "**** NaN detetected in vector v in Peridigm::computeInternalForce().\n");
  }

  // Copy data from mothership vectors to overlap vectors in data manager
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
  }

  // Call the model evaluator
  modelEvaluator->evalModel(workset);

  // Copy force from the data manager to the mothership vector
  force->PutScalar(0.0);
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    scratch->PutScalar(0.0);
    blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
    force->Update(1.0, *scratch, 1.0);
  }
  scratch->PutScalar(0.0);

  // Run some checks to make sure things haven't gone haywire
  for(int i=0 ; i<force->MyLength() ; ++i){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite((*force)[i]), "**** NaN detetected in force vector in Peridigm::computeInternalForce().\n");
  }

  // convert force density to force
  double *f, *v;
  force->ExtractView(&f);
  volume->ExtractView(&v);
  for(int i=0 ; i < force->MyLength() ; ++i)
    f[i] *= v[i/3];
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

void PeridigmNS::Peridigm::executeNOXQuasiStatic(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  // The tangent map was made with multiphysics compatibility already.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction;

  // The reaction vector is created here and will be longer if multiphysics is enabled
  if(analysisHasMultiphysics)
    reaction = Teuchos::rcp(new Epetra_Vector(unknownsForce->Map()));
  else
    reaction = Teuchos::rcp(new Epetra_Vector(force->Map()));

  // Create vectors that are specific to NOX quasi-statics.
  // These are already sized right for multiphysics simulations because they are created with the tangent map
  Teuchos::RCP<Epetra_Vector> soln = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<NOX::Epetra::Vector> noxSoln = Teuchos::rcp(new NOX::Epetra::Vector(soln, NOX::Epetra::Vector::CreateView));
  soln->PutScalar(0.0);

  Teuchos::RCP<Epetra_Vector> initialGuess = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<NOX::Epetra::Vector> noxInitialGuess = Teuchos::rcp(new NOX::Epetra::Vector(initialGuess, NOX::Epetra::Vector::CreateView));
  initialGuess->PutScalar(0.0);

  noxVelocityAtDOFWithKinematicBC = Teuchos::rcp(new Epetra_Vector(v->Map()));
  if(analysisHasMultiphysics)
    noxPressureVAtDOFWithKinematicBC = Teuchos::rcp(new Epetra_Vector(fluidPressureV->Map()));

  // Initialize velocity to zero
  v->PutScalar(0.0);
  if(analysisHasMultiphysics){
      unknownsV->PutScalar(0.0);
			fluidPressureV->PutScalar(0.0);
	}

  // Pointers into mothership vectors
  double *xPtr, *uPtr, *yPtr, *vPtr, *deltaUPtr;
  double *unknownsUPtr, *unknownsYPtr, *unknownsVPtr, *unknownsDeltaUPtr;
  double *fluidPressureUPtr, *fluidPressureDeltaUPtr, *fluidPressureYPtr, *fluidPressureVPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  deltaU->ExtractView( &deltaUPtr );
  if(analysisHasMultiphysics){
      fluidPressureU->ExtractView( &fluidPressureUPtr );
      fluidPressureY->ExtractView( &fluidPressureYPtr );
      fluidPressureV->ExtractView( &fluidPressureVPtr );
      fluidPressureDeltaU->ExtractView( &fluidPressureDeltaUPtr );
      unknownsU->ExtractView( &unknownsUPtr );
      unknownsY->ExtractView( &unknownsYPtr );
      unknownsV->ExtractView( &unknownsVPtr );
      unknownsDeltaU->ExtractView( &unknownsDeltaUPtr );
  }

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
  double timeCurrent = timeSteps[0];
  double timePrevious = timeCurrent;

  // Apply BC at time zero
  PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
  boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");
  PeridigmNS::Timer::self().startTimer("Apply Body Forces");
  boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

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

    timePrevious = timeCurrent;
    timeCurrent = timeSteps[step];
    double timeIncrement = timeCurrent - timePrevious;
    workset->timeStep = timeIncrement;

    m_noxJacobianUpdateCounter = 0;

    soln->PutScalar(0.0);
    deltaU->PutScalar(0.0);

    // soln is already appropriately sized for multiphysics
    if(analysisHasMultiphysics){
    	unknownsDeltaU->PutScalar(0.0);
			fluidPressureDeltaU->PutScalar(0.0);
    }

    if(analysisHasMultiphysics){
			// Ensure that information from v and fluidPressureV is in unknownsV
			for(int i=0 ; i<initialGuess->MyLength() ; i+=(3+numMultiphysDoFs)){
				for(int j = 0; j<3; ++j){
					(*initialGuess)[i+j] = (*v)[i/(3+numMultiphysDoFs)*3+j]*timeIncrement;
				}
				(*initialGuess)[i+3] = (*fluidPressureV)[i/(3+numMultiphysDoFs)]*timeIncrement;
			}
    }
    else{
        // Use a predictor based on the velocity from the previous load step
				//
        for(int i=0 ; i<initialGuess->MyLength() ; ++i)
            (*initialGuess)[i] = (*v)[i]*timeIncrement;
    }

    v->PutScalar(0.0); 
    if(analysisHasMultiphysics){
    	unknownsV->PutScalar(0.0);
			fluidPressureV->PutScalar(0.0);
    }

    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(initialGuess);
    // Apply the displacement increment and update the temperature field
    // Note that the soln vector was created to be compatible with the tangent matrix, and hence needed to be
    // constructed with an Epetra_Map.  The mothership vectors were constructed with an Epetra_BlockMap, and it's
    // this map that the boundary and intial condition manager expects.  So, make sure that the boundary and initial
    // condition manager gets the right type of vector.

    PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent, timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");

    // For NOX, add the increment in displacement BC directly into the displacement vector
   //TEUCHOS_TEST_FOR_EXCEPT_MSG(initialGuess->MyLength() != unknownsV->MyLength(), "**** PeridigmNS::Peridigm::executeNOXQuasiStatic() initialGuess vector different length than unknownsV.\n");
    if(analysisHasMultiphysics){
			for(int i=0 ; i<u->MyLength() ; ++i)
				uPtr[i] += deltaUPtr[i];

			for(int i=0 ; i<fluidPressureU->MyLength() ; ++i){
				fluidPressureUPtr[i] += fluidPressureDeltaUPtr[i];
			}

      for(int i=0 ; i<unknownsU->MyLength() ; i+=(3+numMultiphysDoFs)){
				for(int j = 0; j<3; ++j){
			  	unknownsUPtr[i+j] = uPtr[i*3/(3+numMultiphysDoFs) + j];
				}
        unknownsUPtr[i+3] = fluidPressureUPtr[i/(3+numMultiphysDoFs)];
			}
    }
    else{
    	for(int i=0 ; i<u->MyLength() ; ++i)
        	uPtr[i] += deltaUPtr[i];
    }

    // Note:  applyBoundaryConditions() sets the velocity as well as the displacement.
    //        This needs to be stored because NOX returns a deltaU of zero for these dof
    *noxVelocityAtDOFWithKinematicBC = *v;

    if(analysisHasMultiphysics)
        *noxPressureVAtDOFWithKinematicBC = *fluidPressureV;

    // evaluate the external (body) forces:
    PeridigmNS::Timer::self().startTimer("Apply Body Forces");
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent, timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

    *soln = *initialGuess;

    double toleranceMultiplier = 1.0;
    if(!useAbsoluteTolerance){
      // compute the vector of reactions, i.e., the forces corresponding to degrees of freedom for which kinematic B.C. are applied
      if(analysisHasMultiphysics){
				for(int i=0 ; i<y->MyLength() ; ++i)
					yPtr[i] = xPtr[i] + uPtr[i];

				for(int i=0 ; i<fluidPressureY->MyLength() ; ++i){
					fluidPressureYPtr[i] = fluidPressureUPtr[i];
				}

        for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
					for(int j = 0; j<3; ++j){
						unknownsYPtr[i+j] = yPtr[i/(3+numMultiphysDoFs)*3 + j];
						unknownsUPtr[i+j] = uPtr[i/(3+numMultiphysDoFs)*3 + j];
					}
						unknownsYPtr[i+3] = fluidPressureYPtr[i/(3+numMultiphysDoFs)];
						unknownsUPtr[i+3] = fluidPressureUPtr[i/(3+numMultiphysDoFs)];
				}
      }
      else{
      	for(int i=0 ; i<y->MyLength() ; ++i)
         	yPtr[i] = xPtr[i] + uPtr[i];
      }

      computeQuasiStaticResidual(residual);

      if(analysisHasMultiphysics)
        boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(unknownsForce, reaction);
      else
        boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(force, reaction);

      // convert force density to force
      if(analysisHasMultiphysics){
        for(int i=0 ; i<reaction->MyLength() ; ++i)
            (*reaction)[i] *= (*volume)[i/(3 + numMultiphysDoFs)];
      }
      else{
        for(int i=0 ; i<reaction->MyLength() ; ++i)
            (*reaction)[i] *= (*volume)[i/3];
      }

      double reactionNorm2;
      reaction->Norm2(&reactionNorm2);
      toleranceMultiplier = reactionNorm2;
    }
    double residualTolerance = tolerance*toleranceMultiplier;

    // Get the linear solver parameters from the proper sublist
    // \todo Handle all allowable "Direction" settings
    Teuchos::RCP<Teuchos::ParameterList> linearSystemParams = Teuchos::rcp(new Teuchos::ParameterList);
    string directionMethod = noxQuasiStaticParams->sublist("Direction").get<string>("Method");
    if(directionMethod == "Newton")
      linearSystemParams = Teuchos::rcpFromRef( noxQuasiStaticParams->sublist("Direction").sublist("Newton").sublist("Linear Solver") );
    else if(directionMethod == "NonlinearCG")
      linearSystemParams = Teuchos::rcpFromRef( noxQuasiStaticParams->sublist("Direction").sublist("Nonlinear CG").sublist("Linear Solver") );
    else{
      TEUCHOS_TEST_FOR_EXCEPT_MSG(directionMethod != "Newton" && directionMethod != "NonlinearCG", "\n****Error:  User-supplied NOX Direction currently not supported by Peridigm.\n");
    }

    Material::JacobianType peridigmPreconditioner = Material::FULL_MATRIX;
    if(solverParams->isParameter("Peridigm Preconditioner")){
      std::string peridigmPreconditionerStr = solverParams->get<std::string>("Peridigm Preconditioner");
      if(peridigmPreconditionerStr == "Full Tangent")
        peridigmPreconditioner = Material::FULL_MATRIX;
      else if(peridigmPreconditionerStr == "Block 3x3")
        peridigmPreconditioner = Material::BLOCK_DIAGONAL;
      else if(peridigmPreconditionerStr == "None")
        peridigmPreconditioner = Material::BLOCK_DIAGONAL;
      else
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "\n****Error:  Unrecognized Peridigm Preconditioner, must be \"Full Tangent\", \"Block 3x3\", or \"None\".\n");
    }
    bool isMatrixFree(false);
    if(linearSystemParams->isParameter("Jacobian Operator")){
      std::string jacobianOperator = linearSystemParams->get<std::string>("Jacobian Operator");
      if(jacobianOperator == "Matrix-Free")
        isMatrixFree = true;
    }

    // Construct the NOX linear system
    // This object provides an interface that is used by the NOX nonlinear solvers to query the tangent matrix, or some approximation of that matrix
    Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO> linearSystemAztecOO;
    Teuchos::RCP<NOX::Epetra::Interface::Required> noxRequiredInterface = Teuchos::RCP<NOX::Epetra::Interface::Required>(this, false);;
    const NOX::Epetra::Vector& noxCloneVector = *soln;

    // If we are using the tangent matrix as the Jacobian and/or preconditioner, then we want to provide a residual
    // evaluation ("required" interface) and the tangent
    if(!isMatrixFree || peridigmPreconditioner != Material::BLOCK_DIAGONAL){
      if(step == 1 && peridigmComm->MyPID() == 0)
        cout << "NOX initialized with standard Jacobian operator\n" << endl;
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> noxJacobianInterface = Teuchos::RCP<NOX::Epetra::Interface::Jacobian>(this, false);
      Teuchos::RCP<Epetra_RowMatrix> noxJacobian = tangent;
      linearSystemAztecOO = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                                              *linearSystemParams,
                                                                              noxRequiredInterface,
                                                                              noxJacobianInterface,
                                                                              noxJacobian,
                                                                              noxCloneVector));
    }
    // If we are using the Jacobian Free Newton Krylov approach with the block 3x3 preconditioner, then we want to provide the
    // internal force operator (a.k.a. the "required" interface) and also the block 3x3 preconditioner (the use of the preconditioner
    // is a optional and is specified in the input deck via "Peridigm Preconditioner Type" = "Full Tangent", "Block 3x3", or "None")
    else{
      if(step == 1 && peridigmComm->MyPID() == 0)
        cout << "NOX initialized with matrix-free Jacobian operator\n" << endl;

      // We always want to reuse the preconditioner (it there is one) for matrix-free solves
      // Default to 200 times within a nonlinear solve
      maxAgePeridigmPreconditioner = noxQuasiStaticParams->get<int>("Max Age Of Prec", 200);
      if(noxQuasiStaticParams->isParameter("Preconditioner Reuse Policy")){
        std::string reusePolicy = noxQuasiStaticParams->get<std::string>("Preconditioner Reuse Policy");
        TEUCHOS_TEST_FOR_EXCEPT_MSG(reusePolicy != "Reuse", "\n****Error:  Peridigm only supports \"Preconditioner Reuse Policy = Reuse\" for NOX with Jacobian Free Newton Krylov.\n");
      }

      // For matrix-free solves with the block 3x3 preconditioner, the only valid NOX->Direction->Newton->Linear Solver->Preconditioner options
      // are "None" and "User Defined"
      if(linearSystemParams->isParameter("Preconditioner")){
        std::string linSysPreconditioner = linearSystemParams->get<std::string>("Preconditioner");
        TEUCHOS_TEST_FOR_EXCEPT_MSG(linSysPreconditioner != "User Defined" && linSysPreconditioner != "None",
                                  "\n****Error:  Peridigm only supports \"Preconditioner = User Defined\" and \"Preconditioner = None\" for NOX with Jacobian Free Newton Krylov and the Peridigm Block 3x3 preconditioner.\n");
        // THIS IS A HACK TO GET AROUND APPARENT NOX LINEAR SYSTEM PARAMETER GLITCH
        // We are providing a preconditioner via the LinearSystemAztecOO constructor
        // For whatever reason (user error, probably), this only works if Preconditioner is set to one of the standard
        // preconditioner types.  Logically, we want the preconditioner type to be "User Defined", but this causes a seg
        // fault, probably because NOX looks for certain Epetra_Operator methods that are not implemented for our preconditioner.
        // Setting the preconditioner to "AztecOO" seems to work; the Peridigm preconditioner is used directly (AztecOO does not
        // appear to actually be involved).
        if(linSysPreconditioner == "User Defined")
          linearSystemParams->set("Preconditioner", "AztecOO");
      }

      Teuchos::RCP<NOX::Epetra::MatrixFree> matrixFreeJacobianOperator
        = Teuchos::RCP<NOX::Epetra::MatrixFree>(new NOX::Epetra::MatrixFree(printParams,
                                                                            noxRequiredInterface,
                                                                            noxCloneVector));
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> noxJacobianInterface = matrixFreeJacobianOperator;
      Teuchos::RCP<Epetra_Operator> noxJacobian = matrixFreeJacobianOperator;
      Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> noxPreconditionerInterface = Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>(this, false);
      Teuchos::RCP<Epetra_Operator> noxPreconditioner = tangent;
      linearSystemAztecOO = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                                              *linearSystemParams,
                                                                              noxJacobianInterface,
                                                                              noxJacobian,
                                                                              noxPreconditionerInterface,
                                                                              noxPreconditioner,
                                                                              noxCloneVector));
    }

    // Create the Group
    NOX::Epetra::Vector noxInitialGuess(soln, NOX::Epetra::Vector::CreateView);
    Teuchos::RCP<NOX::Epetra::Group> noxGroup = Teuchos::rcp(new NOX::Epetra::Group(printParams, noxRequiredInterface, noxInitialGuess, linearSystemAztecOO));

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

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(noxGroup, combo, noxQuasiStaticParams);
    NOX::StatusTest::StatusType noxSolverStatus = NOX::StatusTest::Unevaluated;
    int solverIteration = 1;
    agePeridigmPreconditioner = 0;

    if(peridigmComm->MyPID() == 0)
      cout << "Load step " << step << ", initial time = " << timePrevious << ", final time = " << timeCurrent <<
        ", convergence criterion = " << residualTolerance << endl;

    // Solve!
    while(noxSolverStatus != NOX::StatusTest::Converged && noxSolverStatus != NOX::StatusTest::Failed){

      // Track the total number of iterations taken over the simulation
      *nonlinearSolverIterations += 1;

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
    if(analysisHasMultiphysics){
        for(int i=0 ; i<unknownsV->MyLength() ; ++i)
            (*unknownsV)[i] = finalSolution[i]/timeIncrement;

        // Add the converged displacement increment to the displacement
				//std::vector<double> whatWas(unknownsU->MyLength()); 
				//std::vector<double> theChange(unknownsU->MyLength()); 

        for(int i=0 ; i<unknownsU->MyLength() ; ++i){
					  //whatWas[i] = (*unknownsU)[i];
            (*unknownsU)[i] += finalSolution[i];
						//theChange[i] = whatWas[i] - (*unknownsU)[i];
				}

				/*  
				double changeSeen = 0.0;
        for(int i=0 ; i<unknownsU->MyLength() ; ++i){
					changeSeen += (theChange[i]*theChange[i]);
				}
				std::cout << "**** change seen is: " << changeSeen << std::endl;
				*/

				// Update the uncombined mothership vectors from the combined mothership vectors
				// store velocity, displacement and fluid pressure velocity and displacement 
				// analogues. 
  			for(int i=0 ; i<v->MyLength() ; i+=3){
					for(int j=0 ; j<3 ; ++j){
						(*v)[i+j] = (*unknownsV)[i/3*(3+numMultiphysDoFs) + j];
						(*u)[i+j] = (*unknownsU)[i/3*(3+numMultiphysDoFs) + j];
					}
					(*fluidPressureV)[i/3] = (*unknownsV)[i/3*(3+numMultiphysDoFs) + 3];
					(*fluidPressureU)[i/3] = (*unknownsU)[i/3*(3+numMultiphysDoFs) + 3];
				}

				//add nox velocity vectors to uncombined velocity vectors and store back in combined vector
				v->Update(1.0, *noxVelocityAtDOFWithKinematicBC, 1.0);
				fluidPressureV->Update(1.0, *noxPressureVAtDOFWithKinematicBC, 1.0);

				for(int i=0 ; i<v->MyLength() ; i+=3){
					for(int j=0 ; j<3 ; ++j){
						(*unknownsV)[i/3*(3+numMultiphysDoFs) + j] = (*v)[i+j];
					}
					(*unknownsV)[i/3*(3+numMultiphysDoFs) + 3] = (*fluidPressureV)[i/3];
				}

    }
    else{

			for(int i=0 ; i<v->MyLength() ; ++i)
				(*v)[i] = finalSolution[i]/timeIncrement;

			v->Update(1.0, *noxVelocityAtDOFWithKinematicBC, 1.0);

			// Add the converged displacement increment to the displacement
        for(int i=0 ; i<u->MyLength() ; ++i)
            (*u)[i] += finalSolution[i];
    }

    // Write output for completed load step
    PeridigmNS::Timer::self().startTimer("Output");
    synchDataManagers();
    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->updateState();
  }
  // if(peridigmComm->MyPID() == 0)
  //   cout << endl;
}

void PeridigmNS::Peridigm::executeQuasiStatic(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  // Create vectors that are specific to quasi-statics.
  // These must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction;

  if(analysisHasMultiphysics)
	 reaction = Teuchos::rcp(new Epetra_Vector(unknownsForce->Map()));
  else
	 reaction = Teuchos::rcp(new Epetra_Vector(force->Map()));

  const bool disableHeuristics = solverParams->get("Disable Heuristics", false);
  if(disableHeuristics && peridigmComm->MyPID() == 0)
    cout << "\nUser has disabled solver heuristics.\n" << endl;

  // Vector for predictor
  //Epetra_Vector predictor(v->Map());
  Teuchos::RCP<Epetra_Vector> predictor;
  if(analysisHasMultiphysics)
  	predictor = Teuchos::rcp(new Epetra_Vector(unknownsV->Map()));
  else
  	predictor = Teuchos::rcp(new Epetra_Vector(v->Map()));

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
  // Pointers into mothership vectors
  double *xPtr, *uPtr, *yPtr, *vPtr, *aPtr, *deltaUPtr;
  double *unknownsUPtr, *unknownsYPtr, *unknownsVPtr, *unknownsDeltaUPtr;
  double *fluidPressureUPtr, *fluidPressureDeltaUPtr, *fluidPressureYPtr, *fluidPressureVPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  deltaU->ExtractView( &deltaUPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  a->ExtractView( &aPtr );
  if(analysisHasMultiphysics){
      fluidPressureU->ExtractView( &fluidPressureUPtr );
      fluidPressureY->ExtractView( &fluidPressureYPtr );
      fluidPressureV->ExtractView( &fluidPressureVPtr );
      fluidPressureDeltaU->ExtractView( &fluidPressureDeltaUPtr );
      unknownsU->ExtractView( &unknownsUPtr );
      unknownsY->ExtractView( &unknownsYPtr );
      unknownsV->ExtractView( &unknownsVPtr );
      unknownsDeltaU->ExtractView( &unknownsDeltaUPtr );
  }

  // Initialize velocity to zero
  v->PutScalar(0.0);
  if(analysisHasMultiphysics){
	  unknownsV->PutScalar(0.0);
	  fluidPressureV->PutScalar(0.0);
  }

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

  // Adaptive load-stepping parameters
  double timeFinal = timeSteps[timeSteps.size()-1];
  bool solverFailedToConverge = false;
  bool adaptiveLoadStepping = false;
  bool switchToExplicit = false;
  int maxSolverFailureInOneStep;
  int maxTotalSolverFailure;
  bool reduceAllSteps = false;
  bool adaptiveOutputFrequency = false;
  Teuchos::RCP< Teuchos::ParameterList > adaptiveQSparams;
  Teuchos::RCP< Teuchos::ParameterList > verletSolverParams;
  if( quasiStaticParams->isSublist("Adaptive Load-Stepping") ){
    adaptiveLoadStepping = true;
    adaptiveQSparams = sublist(quasiStaticParams, "Adaptive Load-Stepping", true);
    maxSolverFailureInOneStep = adaptiveQSparams->get<int>("Maximum load step reductions in one step");
    maxTotalSolverFailure = adaptiveQSparams->get<int>("Maximum total load step reductions");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(maxTotalSolverFailure < maxSolverFailureInOneStep, "**** 'Maximum total load step reductions' cannot be smaller than 'Maximum load step reductions in one step'. ****");
    if( adaptiveQSparams->isParameter("Reduce all remaining load steps") ){
      //If one load step is reduced, the remainder of the steps is also reduced
      reduceAllSteps = adaptiveQSparams->get<bool>("Reduce all remaining load steps");
    }
    if( adaptiveQSparams->isParameter("Adaptive output frequency") ){
      //To multiply output frequency while reducing the load step size
      adaptiveOutputFrequency = adaptiveQSparams->get<bool>("Adaptive output frequency");
    }
    if( adaptiveQSparams->isSublist("Switch to Verlet") ){
      switchToExplicit = true;
      verletSolverParams = sublist(adaptiveQSparams, "Switch to Verlet", true);
    }
  }
  int solverFailureInThisStep = 0;
  int totalSolverFailure = 0;
  bool convergenceHasFailedInThisStep = false;
  bool failedQS = false;

  double timeCurrent = timeSteps[0];
  double timePrevious = timeCurrent;

  // Apply BC at time zero for the sake of I/O
  PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
  boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");
  PeridigmNS::Timer::self().startTimer("Apply Body Forces");
  boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  Epetra_Time loadStepCPUTime(*peridigmComm);
  double cumulativeLoadStepCPUTime = 0.0;

  for(int step=1 ; step<(int)timeSteps.size() ; step++){

    if(!adaptiveLoadStepping || !solverFailedToConverge){
      loadStepCPUTime.ResetStartTime();
      timePrevious = timeCurrent;
    }
    solverFailedToConverge = false;

    timeCurrent = timeSteps[step];
    double timeIncrement = timeCurrent - timePrevious;
    workset->timeStep = timeIncrement;

    deltaU->PutScalar(0.0);
		if(analysisHasMultiphysics){
			fluidPressureDeltaU->PutScalar(0.0);
			unknownsDeltaU->PutScalar(0.0);
		}

    // Update nodal positions for nodes with kinematic B.C.
    PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");

    // evaluate the external (body) forces:
    PeridigmNS::Timer::self().startTimer("Apply Body Forces");
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

    // Set the current position and velocity
		for(int i=0 ; i<y->MyLength() ; ++i){
      yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
      vPtr[i] = deltaUPtr[i]/timeIncrement;
    }

		//If true, then synch the content of the combined vectors with that of the uncombined.
    if(analysisHasMultiphysics){
      for(int i=0; i<unknownsDeltaU->MyLength(); i+=(3+numMultiphysDoFs)){
        for(int j=0; j<3; ++j){
          unknownsDeltaUPtr[i+j] = deltaUPtr[i*3/(3+numMultiphysDoFs)+j];
        }
        unknownsDeltaUPtr[i+3] = fluidPressureDeltaUPtr[i/(3+numMultiphysDoFs)];
      }

      for(int i=0 ; i<fluidPressureY->MyLength() ; ++i){
        fluidPressureYPtr[i] = fluidPressureUPtr[i] + fluidPressureDeltaUPtr[i];
        fluidPressureVPtr[i] = fluidPressureDeltaUPtr[i]/timeIncrement;
      }

      for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
        for(int j = 0; j<3; ++j){
          unknownsYPtr[i+j] = yPtr[i/(3+numMultiphysDoFs)*3 + j];
          unknownsVPtr[i+j] = vPtr[i/(3+numMultiphysDoFs)*3 + j];
        }
        unknownsYPtr[i+3] = fluidPressureYPtr[i/(3+numMultiphysDoFs)];
        unknownsVPtr[i+3] = fluidPressureVPtr[i/(3+numMultiphysDoFs)];
      }
    }

    // compute the residual
    double residualNorm = computeQuasiStaticResidual(residual);

    double toleranceMultiplier = 1.0;
    if(!useAbsoluteTolerance){
      // compute the vector of reactions, i.e., the forces corresponding to degrees of freedom for which kinematic B.C. are applied
      if(analysisHasMultiphysics){
        boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(unknownsForce, reaction);
      }
      else{
        boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(force, reaction);
      }
      // convert force density to force
      for(int i=0 ; i<reaction->MyLength() ; ++i){
        (*reaction)[i] *= (*volume)[i/(3+numMultiphysDoFs)];
      }
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
    int numPureNewtonSteps = 50;//8;
    int numPreconditionerSteps = 24;
    int dampedNewtonNumStepsBetweenTangentUpdates = 8;
    double alpha = 0.0;
    while(residualNorm > tolerance*toleranceMultiplier && solverIteration <= maxSolverIterations){

      // Track the total number of iterations taken over the simulation
      *nonlinearSolverIterations += 1;

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
      if(solverIteration == 1 && step > 1 && !disableHeuristics) {
        for(int i=0 ; i<lhs->MyLength() ; ++i)
          (*lhs)[i] = (*predictor)[i]*timeIncrement;
        boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(lhs);
        isConverged = Belos::Converged;
      }
      else{
        // If we reach the specified maximum number of iterations of the nonlinear solver, switch to a damped Newton approach.
        // This should hurt the convergence rate but improve robustness.
        if(solverIteration > numPureNewtonSteps && !dampedNewton && !disableHeuristics){
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

        // Disable the preconditioner if the user specifies disable heuristics
        if(disableHeuristics) usePreconditioner = false;

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

        if(isConverged == Belos::Unconverged && !disableHeuristics){
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
        alpha = disableHeuristics ? 1.0 : quasiStaticsLineSearch(residual, lhs, timeIncrement);
        PeridigmNS::Timer::self().stopTimer("Line Search");

        // Apply increment to nodal positions
				if(analysisHasMultiphysics){
					for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
							for(int j=0 ; j<3 ; ++j){
								unknownsDeltaUPtr[i+j] += alpha*(*lhs)[i+j];
								unknownsYPtr[i+j] = xPtr[i*3/(3+numMultiphysDoFs) +j] + unknownsUPtr[i+j] + unknownsDeltaUPtr[i+j];
								unknownsVPtr[i+j] = unknownsDeltaUPtr[i+j]/timeIncrement;
							}
							unknownsDeltaUPtr[i+3] += alpha*(*lhs)[i+3];
							unknownsYPtr[i+3] = unknownsUPtr[i+3] + unknownsDeltaUPtr[i+3];
							unknownsVPtr[i+3] = unknownsDeltaUPtr[i+3]/timeIncrement;
					}

					for(int i=0 ; i<y->MyLength() ; i+=3){
						for(int j=0 ; j<3 ; ++j){
							deltaUPtr[i+j] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + j];
							yPtr[i+j] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + j];
							vPtr[i+j] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + j];
						}
						fluidPressureDeltaUPtr[i/3] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs)+3];
						fluidPressureYPtr[i/3] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + 3];
						fluidPressureVPtr[i/3] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + 3];
					}
				}
				else{
					for(int i=0 ; i<y->MyLength() ; ++i){
						deltaUPtr[i] += alpha*(*lhs)[i];
						yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
						vPtr[i] = deltaUPtr[i]/timeIncrement;
					}
				}
        // Compute residual
        residualNorm = computeQuasiStaticResidual(residual);

        solverIteration++;
      }
      else{
        if(peridigmComm->MyPID() == 0){
          cout << "\nError:  Belos linear solver failed to converge." << endl;
        }
        residualNorm = 1.0e50;
        break;
      }
    } // end loop of nonlinear iterations

    // If the maximum allowable number of load step reductions has been reached and the residual
    // is within a reasonable tolerance, then just accept the solution and forge ahead.
    // If not, abort the analysis.
    if(residualNorm > tolerance*toleranceMultiplier){
      if(residualNorm < 100.0*tolerance*toleranceMultiplier){
        if(peridigmComm->MyPID() == 0)
          cout << "\nWarning:  Accepting current solution and progressing to next load step.\n" << endl;
      }
      else{

        if(adaptiveLoadStepping){
          solverFailureInThisStep++;
          totalSolverFailure++;

          if(totalSolverFailure > maxTotalSolverFailure || solverFailureInThisStep > maxSolverFailureInOneStep){
            failedQS = true;
            if(!switchToExplicit)
              if(peridigmComm->MyPID() == 0)
                cout << "\nError:  Aborting analysis.\n" << endl;

            break;
          }
          else{
            if(peridigmComm->MyPID() == 0)
              cout << "Warning: Cutting the load step in half and redo the step.\n" << endl;
          }
          solverFailedToConverge = true;
          convergenceHasFailedInThisStep = true;
        }
        else{
          if(peridigmComm->MyPID() == 0)
            cout << "\nError:  Aborting analysis.\n" << endl;
          break;
        }
      }
    }


    // If adaptive load step scheme is chosen by user and QS solver failed
    // to converge, load step is cut in half to retry with new increment.
    if(adaptiveLoadStepping && solverFailedToConverge){
      timeSteps.insert(timeSteps.begin()+step, timePrevious+timeIncrement/2.0);
      if(adaptiveOutputFrequency)
        outputManager->multiplyOutputFrequency(2.0);
      step--;
    }

    if(adaptiveLoadStepping && reduceAllSteps && !solverFailedToConverge && convergenceHasFailedInThisStep){
      // The smallest load step size so far is set to use for the rest of the simulation
      timeSteps.erase(timeSteps.begin()+step+1, timeSteps.end());
      int i=0;
      while(timeSteps[step+i]<timeFinal){
        timeSteps.push_back(timeCurrent + (i+1)*timeIncrement);
        i++;
      }
    }

    if(solverIteration >= maxSolverIterations && peridigmComm->MyPID() == 0)
      cout << "\nWarning:  Nonlinear solver failed to converge in maximum allowable iterations." << endl;

    if(!adaptiveLoadStepping || !solverFailedToConverge){
      // The load step is complete
      // Update internal data and move on to the next load step
      convergenceHasFailedInThisStep = false;
      solverFailureInThisStep = 0;

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
      // Make sure even the non participating vectors are updated
      if(analysisHasMultiphysics){
        for(int i=0 ; i<unknownsU->MyLength() ; i+=(3+numMultiphysDoFs)){
          for(int j=0 ; j<3 ; ++j){
            (*unknownsU)[i+j] += (*unknownsDeltaU)[i+j];
            (*u)[i/(3+numMultiphysDoFs)*3 + j] += (*unknownsDeltaU)[i+j];
            (*deltaU)[i/(3+numMultiphysDoFs)*3 + j] = (*unknownsDeltaU)[i+j];
          }
            (*unknownsU)[i+3] += (*unknownsDeltaU)[i+3];
            (*fluidPressureU)[i/(3+numMultiphysDoFs)] += (*unknownsDeltaU)[i+3];
            (*fluidPressureDeltaU)[i/(3+numMultiphysDoFs)] = (*unknownsDeltaU)[i+3];
        }
        // Store the velocity for use as a predictor in the next load step
        for(int i=0 ; i<unknownsV->MyLength() ; ++i){
          (*predictor)[i] = (*unknownsV)[i];
        }
      }
      else{
        for(int i=0 ; i<u->MyLength() ; ++i){
          (*u)[i] += (*deltaU)[i];
          (*predictor)[i] = (*v)[i];
        }
      }

      // Write output for completed load step
      PeridigmNS::Timer::self().startTimer("Output");
      synchDataManagers();
      outputManager->write(blocks, timeCurrent);
      PeridigmNS::Timer::self().stopTimer("Output");

      // swap state N and state NP1
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->updateState();
    }

  } // end loop over load steps

  if(!switchToExplicit && failedQS){
    if(peridigmComm->MyPID() == 0)
      cout << "\nError: Number of Quasi-Static solver convergence failures was more than the maximum allowable number of failures." << endl;
  }

  if(switchToExplicit && failedQS){
    if(peridigmComm->MyPID() == 0){
      cout << "\nWarning: Number of Quasi-Static solver convergence failures was more than the maximum allowable number of failures.";
      cout << "\nSwitching to explicit time-stepping for the rest of simulation." << endl;
    }
    Teuchos::RCP<Teuchos::ParameterList> explicitSolverParams = Teuchos::rcp(new Teuchos::ParameterList);
    Teuchos::ParameterList& verletParams = explicitSolverParams->sublist("Verlet");
    if(verletSolverParams->isParameter("Fixed dt")){
      double userDefinedTimeStep = verletSolverParams->get<double>("Fixed dt");
      verletParams.set("Fixed dt", userDefinedTimeStep);
    }
    else{
      double safetyFactor = 1.0;
      if(verletSolverParams->isParameter("Safety Factor"))
        safetyFactor = verletSolverParams->get<double>("Safety Factor");
      verletParams.set("Safety Factor", safetyFactor);
    }

    explicitSolverParams->set("Initial Time", timePrevious);
    explicitSolverParams->set("Final Time", timeFinal);

    if(verletSolverParams->isParameter("Output Frequency")){
      int output_frequency = verletSolverParams->get<int>("Output Frequency");
      outputManager->changeOutputFrequency(output_frequency);
    }

    // Restore the values to the converged ones from previous step
		for(int i=0 ; i<y->MyLength() ; ++i){
      yPtr[i] = xPtr[i] + uPtr[i];
      vPtr[i] = (*predictor)[i];
    }

    executeExplicit(explicitSolverParams);
  }

  if(peridigmComm->MyPID() == 0)
    cout << endl;
}

void PeridigmNS::Peridigm::executeImplicitDiffusion(Teuchos::RCP<Teuchos::ParameterList> solverParams) {

  // Create vectors that are specific to implicit diffusion
  // These must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> lhs = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> reaction = Teuchos::rcp(new Epetra_Vector(fluxDivergence->Map()));;

  bool solverVerbose = solverParams->get("Verbose", false);
  Teuchos::RCP<Teuchos::ParameterList> implicitSolverParams = sublist(solverParams, "ImplicitDiffusion", true);
  int maxSolverIterations = implicitSolverParams->get("Maximum Solver Iterations", 10);

  // Determine tolerance
  double tolerance = implicitSolverParams->get("Relative Tolerance", 1.0e-6);
  bool useAbsoluteTolerance = false;
  if(implicitSolverParams->isParameter("Absolute Tolerance")){
    useAbsoluteTolerance = true;
    tolerance = implicitSolverParams->get<double>("Absolute Tolerance");
  }

  // Data for Belos linear solver object
  Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator> linearProblem;
  string linearSolver = implicitSolverParams->get("Belos Linear Solver", "BlockCG");
  Belos::ReturnType isConverged;
  Teuchos::ParameterList belosList;
  belosList.set( "Block Size", 1 );  // Use single-vector iteration
  belosList.set( "Maximum Iterations", implicitSolverParams->get("Belos Maximum Iterations", tangent->NumGlobalRows()) ); // Maximum number of iterations allowed
  belosList.set( "Convergence Tolerance", implicitSolverParams->get("Belos Relative Tolerance", 1.0e-4) ); // Relative convergence tolerance requested
  belosList.set( "Output Frequency", -1 );
  int verbosity = Belos::Errors + Belos::Warnings;
  if( implicitSolverParams->get("Belos Print Status", false) == true ){
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
  if( solverParams->isParameter("Final Time") && implicitSolverParams->isParameter("Number of Load Steps") ){
    double timeInitial = solverParams->get("Initial Time", 0.0);
    double timeFinal = solverParams->get<double>("Final Time");
    int numLoadSteps = implicitSolverParams->get<int>("Number of Load Steps");
    timeSteps.push_back(timeInitial);
    for(int i=0 ; i<numLoadSteps ; ++i)
      timeSteps.push_back(timeInitial + (i+1)*(timeFinal-timeInitial)/numLoadSteps);
  }
  // Case 2:  User provided a list of time steps
  else if( implicitSolverParams->isParameter("Time Steps") ){
    string timeStepString = implicitSolverParams->get<string>("Time Steps");
    istringstream iss(timeStepString);
    copy(istream_iterator<double>(iss),
	 istream_iterator<double>(),
	 back_inserter<vector<double> >(timeSteps));
  }
  else{
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "\n****Error: No valid time step data provided.\n");
  }

  double timeCurrent = timeSteps[0];
  double timePrevious = timeCurrent;
  bool solverFailedToConverge = false;

  // Apply BC at time zero
  PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
  boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");
  PeridigmNS::Timer::self().startTimer("Apply Body Forces");
  boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Body Forces");
  PeridigmNS::Timer::self().startTimer("Apply Initial Conditions");
  boundaryAndInitialConditionManager->applyInitialConditions();
  PeridigmNS::Timer::self().stopTimer("Apply Initial Conditions");

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  Epetra_Time loadStepCPUTime(*peridigmComm);
  double cumulativeLoadStepCPUTime = 0.0;

  for(int step=1 ; step<(int)timeSteps.size() ; step++){

    if(!solverFailedToConverge){
      loadStepCPUTime.ResetStartTime();
      timePrevious = timeCurrent;
    }
    solverFailedToConverge = false;

    timeCurrent = timeSteps[step];
    double timeIncrement = timeCurrent - timePrevious;
    workset->timeStep = timeIncrement;

    //deltaU->PutScalar(0.0);

    boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);

    // copy data into the data manager
		for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
			blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
		}

    // compute the residual
    double residualNorm = computeQuasiStaticResidual(residual);

    double toleranceMultiplier = 1.0;
    if(!useAbsoluteTolerance){
      // compute the vector of reactions, i.e., the forces corresponding to degrees of freedom for which kinematic B.C. are applied
      boundaryAndInitialConditionManager->applyKinematicBC_ComputeReactions(fluxDivergence, reaction);
      for(int i=0 ; i<reaction->MyLength() ; ++i){
        (*reaction)[i] *= (*volume)[i];
      }
      double reactionNorm2;
      reaction->Norm2(&reactionNorm2);
      toleranceMultiplier = reactionNorm2;
    }

    if(peridigmComm->MyPID() == 0)
      cout << "Load step " << step << ", initial time = " << timePrevious << ", final time = " << timeCurrent << ", convergence criterion = " << tolerance*toleranceMultiplier << endl;

    int solverIteration = 1;
    double alpha = 0.0;

    while(residualNorm > tolerance*toleranceMultiplier && solverIteration <= maxSolverIterations){

      // Track the total number of iterations taken over the simulation
      *nonlinearSolverIterations += 1;

      if(peridigmComm->MyPID() == 0)
        cout << "  iteration " << solverIteration << ": residual = " << residualNorm << endl;

      // Compute the tangent
      tangent->PutScalar(0.0);
      PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
      modelEvaluator->evalJacobian(workset);
      int err = tangent->GlobalAssemble();
      TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::executeImplicitDiffusion(), GlobalAssemble() returned nonzero error code.\n");
      PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");
      boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);
      boundaryAndInitialConditionManager->applyKinematicBC_InsertZerosAndSetDiagonal(tangent);
      tangent->Scale(-1.0);

      // Solve linear system
      quasiStaticsSetPreconditioner(linearProblem);
      isConverged = quasiStaticsSolveSystem(residual, lhs, linearProblem, belosSolver);

      if(isConverged == Belos::Converged){

        // Zero out the entries corresponding to the kinematic boundary conditions
        // The solver should have returned zeros, but there may be small errors.
        boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(lhs);

        // placeholder for line search
        alpha = 1.0;

        // Apply increment to temperature
        for(int i=0 ; i<temperature->MyLength() ; ++i){
          (*temperature)[i] += alpha*(*lhs)[i];
        }

        // Compute residual
        residualNorm = computeQuasiStaticResidual(residual);

        solverIteration++;
      }
      else{
        if(peridigmComm->MyPID() == 0){
          cout << "\nError:  Belos linear solver failed to converge." << endl;
        }
        residualNorm = 1.0e50;
        break;
      }
    } // end loop of nonlinear iterations

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

    if(solverIteration >= maxSolverIterations && peridigmComm->MyPID() == 0)
      cout << "\nWarning:  Nonlinear solver failed to converge in maximum allowable iterations." << endl;

    if(!solverFailedToConverge){

      if(peridigmComm->MyPID() == 0)
        cout << "  iteration " << solverIteration << ": residual = " << residualNorm << endl;

      // Print load step timing information
      double CPUTime = loadStepCPUTime.ElapsedTime();
      cumulativeLoadStepCPUTime += CPUTime;
      if(peridigmComm->MyPID() == 0)
        cout << setprecision(2) << "  cpu time for load step = " << CPUTime << " sec., cumulative cpu time = " << cumulativeLoadStepCPUTime << " sec.\n" << endl;

      // Write output for completed load step
      PeridigmNS::Timer::self().startTimer("Output");
      synchDataManagers();
      outputManager->write(blocks, timeCurrent);
      PeridigmNS::Timer::self().stopTimer("Output");

      // swap state N and state NP1
      for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
        blockIt->updateState();
    }
  } // end loop over load steps

  if(peridigmComm->MyPID() == 0)
    cout << endl;
}

void PeridigmNS::Peridigm::quasiStaticsSetPreconditioner(Belos::LinearProblem<double,Epetra_MultiVector,Epetra_Operator>& linearProblem) {
  Ifpack IFPFactory;
  Teuchos::ParameterList ifpackList;

  std::string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1, param is ignored.
  Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( IFPFactory.Create(PrecType, &(*tangent), OverlapLevel) );
  // ifpackList.set("fact: drop tolerance", 1e-9);
  ifpackList.set("fact: ilut level-of-fill", 0);
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
  Teuchos::RCP<Epetra_Vector> tempVector;
  if(analysisHasMultiphysics)
  	tempVector = Teuchos::rcp(new Epetra_Vector(*unknownsDeltaU));
  else
  	tempVector = Teuchos::rcp(new Epetra_Vector(*deltaU));

  // Pointers into mothership vectors
  double *xPtr, *uPtr, *yPtr, *vPtr, *deltaUPtr, *lhsPtr, *residualPtr;
  double *unknownsUPtr, *unknownsYPtr, *unknownsVPtr, *unknownsDeltaUPtr;
  double *fluidPressureUPtr, *fluidPressureDeltaUPtr, *fluidPressureYPtr, *fluidPressureVPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  deltaU->ExtractView( &deltaUPtr );
  lhs->ExtractView(&lhsPtr);
  residual->ExtractView(&residualPtr);
  if(analysisHasMultiphysics){
      fluidPressureU->ExtractView( &fluidPressureUPtr );
      fluidPressureY->ExtractView( &fluidPressureYPtr );
      fluidPressureV->ExtractView( &fluidPressureVPtr );
      fluidPressureDeltaU->ExtractView( &fluidPressureDeltaUPtr );
      unknownsU->ExtractView( &unknownsUPtr );
      unknownsY->ExtractView( &unknownsYPtr );
      unknownsV->ExtractView( &unknownsVPtr );
      unknownsDeltaU->ExtractView( &unknownsDeltaUPtr );
  }

  // compute the current residual
  double unperturbedResidualNorm = computeQuasiStaticResidual(residual);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(unperturbedResidualNorm), "**** NaN detected in residual calculation in quasiStaticsLineSearch().\n");
  if(unperturbedResidualNorm == 0.0)
    return 0.0;

  double bestAlpha = 1.0;
  double bestResidual = 1.0e50;

  vector<double> candidateAlphas;

  // a systematic guess for alpha
  double epsilon = 1.0e-4;
  if(analysisHasMultiphysics){
		for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
			for(int j=0 ; j<3 ; ++j){
				unknownsDeltaUPtr[i+j] += epsilon*lhsPtr[i+j];
      	unknownsYPtr[i+j] = xPtr[i*3/(3+numMultiphysDoFs)+j] + unknownsUPtr[i+j] + unknownsDeltaUPtr[i+j];
				unknownsVPtr[i+j] = unknownsDeltaUPtr[i+j]/dt;
			}
			unknownsDeltaUPtr[i+3] += epsilon*lhsPtr[i+3];
      unknownsYPtr[i+3] = unknownsUPtr[i+3] + unknownsDeltaUPtr[i+3];
			unknownsVPtr[i+3] = unknownsDeltaUPtr[i+3]/dt;
		}

 		for(int i=0 ; i<y->MyLength() ; i+=3){
			for(int j=0 ; j<3 ; ++j){
				deltaUPtr[i+j] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + j];
				yPtr[i+j] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + j];
        vPtr[i+j] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + j];
			}
			fluidPressureDeltaUPtr[i/3] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + 3];
			fluidPressureYPtr[i/3] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + 3];
			fluidPressureVPtr[i/3] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + 3];
    }
  }
  else{

 		for(int i=0 ; i<y->MyLength() ; ++i){
    		deltaUPtr[i] += epsilon*lhsPtr[i];
    		yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
    		vPtr[i] = deltaUPtr[i]/dt;
  	}
  }

  Teuchos::RCP<Epetra_Vector> perturbedResidual = Teuchos::rcp(new Epetra_Vector(*residual));
  computeQuasiStaticResidual(perturbedResidual);
  double SR, SPerturbedR;
  lhs->Dot(*residual, &SR);
  lhs->Dot(*perturbedResidual, &SPerturbedR);
  double tempAlpha = -1.0*epsilon*SR/(SPerturbedR - SR);
  if(tempAlpha > -0.1 && tempAlpha < 10.0)
    candidateAlphas.push_back(tempAlpha);
  if(analysisHasMultiphysics)
	{
  	*unknownsDeltaU = *tempVector;
	}
  else{
  	*deltaU = *tempVector;
	}

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

    if(analysisHasMultiphysics){
			for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
				for(int j=0 ; j<3 ; ++j){
					unknownsDeltaUPtr[i+j] += epsilon*lhsPtr[i+j];
					unknownsYPtr[i+j] = xPtr[i*3/(3+numMultiphysDoFs)+j] + unknownsUPtr[i+j] + unknownsDeltaUPtr[i+j];
					unknownsVPtr[i+j] = unknownsDeltaUPtr[i+j]/dt;
				}
				unknownsDeltaUPtr[i+3] += epsilon*lhsPtr[i+3];
				unknownsYPtr[i+3] = unknownsUPtr[i+3] + unknownsDeltaUPtr[i+3];
				unknownsVPtr[i+3] = unknownsDeltaUPtr[i+3]/dt;
			}

			for(int i=0 ; i<y->MyLength() ; i+=3){
				for(int j=0 ; j<3 ; ++j){
					deltaUPtr[i+j] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + j];
					yPtr[i+j] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + j];
					vPtr[i+j] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + j];
				}
				fluidPressureDeltaUPtr[i/3] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + 3];
				fluidPressureYPtr[i/3] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + 3];
				fluidPressureVPtr[i/3] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + 3];
			}
		}
    else{
			for(int i=0 ; i<y->MyLength() ; ++i){
					deltaUPtr[i] += alpha*lhsPtr[i];
					yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
					vPtr[i] = deltaUPtr[i]/dt;
			}
    }

    double residualNorm = computeQuasiStaticResidual(residual);
		if(analysisHasMultiphysics){
			*unknownsDeltaU = *tempVector;
		}
		else{
    	*deltaU = *tempVector;
		}
    if(std::isfinite(residualNorm) && residualNorm < bestResidual){
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
		if(analysisHasMultiphysics){
			for(int i=0 ; i<unknownsY->MyLength() ; i+=(3+numMultiphysDoFs)){
				for(int j=0 ; j<3 ; ++j){
					unknownsDeltaUPtr[i+j] += alpha*lhsPtr[i+j];
					unknownsYPtr[i+j] = xPtr[i*3/(3+numMultiphysDoFs)+j] + unknownsUPtr[i+j] + unknownsDeltaUPtr[i+j];
					unknownsVPtr[i+j] = unknownsDeltaUPtr[i+j]/dt;
				}
				unknownsDeltaUPtr[i+3] += alpha*lhsPtr[i+3];
				unknownsYPtr[i+3] = unknownsUPtr[i+3] + unknownsDeltaUPtr[i+3];
				unknownsVPtr[i+3] = unknownsDeltaUPtr[i+3]/dt;
			}

			for(int i=0 ; i<y->MyLength() ; i+=3){
				for(int j=0 ; j<3 ; ++j){
					deltaUPtr[i+j] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + j];
					yPtr[i+j] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + j];
					vPtr[i+j] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + j];
				}
				fluidPressureDeltaUPtr[i/3] = unknownsDeltaUPtr[i/3*(3+numMultiphysDoFs) + 3];
				fluidPressureYPtr[i/3] = unknownsYPtr[i/3*(3+numMultiphysDoFs) + 3];
				fluidPressureVPtr[i/3] = unknownsVPtr[i/3*(3+numMultiphysDoFs) + 3];
			}
    }
    else{
			for(int i=0 ; i<y->MyLength() ; ++i){
					deltaUPtr[i] += alpha*lhsPtr[i];
					yPtr[i] = xPtr[i] + uPtr[i] + deltaUPtr[i];
					vPtr[i] = deltaUPtr[i]/dt;
			}
    }

    double residualNorm = computeQuasiStaticResidual(residual);
	if(analysisHasMultiphysics){
      *unknownsDeltaU = *tempVector;
	}
	else{
    	*deltaU = *tempVector;
	}
    if(std::isfinite(residualNorm)){
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
	//TODO: Eliminate fluid pressure V

  // Create vectors that are specific to implicit dynamics
  // The residual must use the same map as the tangent matrix, which is an Epetra_Map and is not consistent
  // with the Epetra_BlockMap used for the mothership multivector.
  Teuchos::RCP<Epetra_Vector> residual = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> displacementIncrement = Teuchos::rcp(new Epetra_Vector(tangent->Map()));
  Teuchos::RCP<Epetra_Vector> un = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> vn = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> an = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  Teuchos::RCP<Epetra_Vector> fluidPressureUn;
  Teuchos::RCP<Epetra_Vector> fluidPressureVn;

  //Don't bother allocating memory for the multiphysics variables unless we need to.
  if(analysisHasMultiphysics){
     fluidPressureUn = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
     fluidPressureVn = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  }

  Teuchos::RCP<Teuchos::ParameterList> implicitParams = sublist(solverParams, "Implicit", true);
  double timeInitial = solverParams->get("Initial Time", 0.0);
  double timeFinal = solverParams->get<double>("Final Time");
  double timeCurrent = timeInitial;
  double timePrevious = timeCurrent;
  double absoluteTolerance       = implicitParams->get("Absolute Tolerance", 1.0e-6);
  int maxSolverIterations        = implicitParams->get("Maximum Solver Iterations", 10);
  double dt                      = implicitParams->get<double>("Fixed dt");
  double beta                    = implicitParams->get("Beta", 0.25);
  double gamma                   = implicitParams->get("Gamma", 0.50);
  workset->timeStep = dt;
  double dt2 = dt*dt;
  int nsteps = (int)floor((timeFinal-timeInitial)/dt);

  // Pointer index into sub-vectors for use with BLAS
  double *xPtr, *uPtr, *yPtr, *vPtr, *aPtr;
  double *fluidPressureUPtr, *fluidPressureYPtr, *fluidPressureVPtr;
  x->ExtractView( &xPtr );
  u->ExtractView( &uPtr );
  y->ExtractView( &yPtr );
  v->ExtractView( &vPtr );
  a->ExtractView( &aPtr );
  if(analysisHasMultiphysics){
  	fluidPressureU->ExtractView( &fluidPressureUPtr );
	fluidPressureY->ExtractView( &fluidPressureYPtr );
	fluidPressureV->ExtractView( &fluidPressureVPtr );
  }

  // Data for linear solver object
  Belos::LinearProblem<double, Epetra_MultiVector, Epetra_Operator> linearProblem;
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
  Teuchos::RCP<Epetra_Vector> u2;
  Teuchos::RCP<Epetra_Vector> v2;
  Teuchos::RCP<Epetra_Vector> fluidPressureU2;

  u2 = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  v2 = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  //Don't bother allocating memory for multiphysics variables unless we need to.
  if(analysisHasMultiphysics){
	fluidPressureU2 = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  }

  // \todo Put in mothership.
  //Teuchos::RCP<Epetra_Vector> deltaU = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));

  // Apply BC at time zero
  PeridigmNS::Timer::self().startTimer("Apply Kinematic B.C.");
  boundaryAndInitialConditionManager->applyBoundaryConditions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Kinematic B.C.");
  PeridigmNS::Timer::self().startTimer("Apply Body Forces");
  boundaryAndInitialConditionManager->applyForceContributions(timeCurrent,timePrevious);
  PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

  // Write initial configuration to disk
  PeridigmNS::Timer::self().startTimer("Output");
  synchDataManagers();
  outputManager->write(blocks, timeCurrent);
  PeridigmNS::Timer::self().stopTimer("Output");

  for(int step=0; step<nsteps ; step++){

    if(peridigmComm->MyPID() == 0)
      cout << "Load step " << step << ", initial time = " << step*dt << ", final time = " << (step+1)*dt << endl;

    *un = *u;
    *vn = *v;
    *an = *a;

    if(analysisHasMultiphysics){
      *fluidPressureUn = *fluidPressureU;
      *fluidPressureVn = *fluidPressureV;
    }

    //u2 = un + dt*vn + 0.5*dt*dt*(1-2*beta)*an
    u2->Update(1.0, *un, 0.0);
    u2->Update(dt, *vn, 0.5*dt2*(1.0-2.0*beta), *an, 1.0);
    //v2 = vn + dt*(1-gamma)*an
    v2->Update(1.0, *vn, dt*(1.0-gamma), *an, 0.0);

    if(analysisHasMultiphysics){
      std::cout << "This happens and may not supposed to." << std::endl;
      // For the fluids part predictor, do backward Euler step.
      // rho*c*dp[x,t]/dt = fluidFlow[p[x,t],t]
      // sum both sides wrt t, from t_k to t_k+1, apply right side rectangle rule
      // p[x, t_k+1] - p[x, t_k] = delta_t*fluidFlow[p[x,t_k+1], t_k+1]/rho/c
      // Solve for updated pressure
      // p[x, t_k+1]_i+1 = delta_t*fluidFlow[p[x, t_k+1]_i, t_k+1]/rho/c + p[x, t_k]
      // where initially
      // p[x, t_k+1]_0 = p[x, t_k]

      // Use the rate of change of pressure now to predict fluid pressure
      fluidPressureU2->Update(1.0,*fluidPressureUn, 0.0);
      fluidPressureU2->Update(dt, *fluidPressureVn, 1.0);
    }

    // Fill the owned vectors with probe data
    // Assign predictor (use u2)
    u->Update(1.0, *u2, 0.0);
    // a = (1.0/(beta*dt*dt))*(u_np1 - u2);
    // a will be zero unless a different predictor is used, so do the following computation anyway
    a->Update(1.0, *u, -1.0, *u2, 0.0);
    a->Scale(1.0/(beta*dt2));
    // v = v2 + dt*gamma*an
    v->Update(1.0, *v2, dt*gamma, *a, 0.0);
    // Update y to be consistent with u
    y->Update(1.0, *x, 1.0, *u, 0.0);

    if(analysisHasMultiphysics){
      fluidPressureU->Update(1.0,*fluidPressureU2,0.0);
      fluidPressureY->Update(1.0, *fluidPressureU, 0.0);
    }

    // Copy data from mothership vectors to overlap vectors in data manager
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(a, accelerationFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
      if(analysisHasMultiphysics){
        blockIt->importData(fluidPressureU, fluidPressureUFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(fluidPressureY, fluidPressureYFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(fluidPressureV, fluidPressureVFieldId, PeridigmField::STEP_NP1, Insert);
      }
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

    // Update forces based on new positions
    PeridigmNS::Timer::self().startTimer("Internal Force");
    modelEvaluator->evalModel(workset);
    PeridigmNS::Timer::self().stopTimer("Internal Force");

    // Copy force from the data manager to the mothership vector
    force->PutScalar(0.0);
    if(analysisHasMultiphysics){
	fluidFlow->PutScalar(0.0);
	for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
 		scratch->PutScalar(0.0);
		scalarScratch->PutScalar(0.0);

		blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
		blockIt->exportData(scalarScratch, fluidFlowDensityFieldId, PeridigmField::STEP_NP1, Add);

		force->Update(1.0, *scratch, 1.0);
		fluidFlow->Update(1.0, *scalarScratch, 1.0);
	}
	scalarScratch->PutScalar(0.0);
    }
    else{
	for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
		scratch->PutScalar(0.0);

		blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);

		force->Update(1.0, *scratch, 1.0);
	}
    }
    scratch->PutScalar(0.0);

    // evaluate the external (body) forces:
    PeridigmNS::Timer::self().startTimer("Apply Body Forces");
    boundaryAndInitialConditionManager->applyForceContributions(timeCurrent, timePrevious);
    PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

    // Compute the residual
    // residual = beta*dt*dt*(M*a - force)
    // Note that due to restrictions to CrsMatrix, the residual has a different (but equivalent) map
    // than the force and acceleration
    if(analysisHasMultiphysics){
      for(int i=0 ; i<residual->MyLength() ; i+=(3+numMultiphysDoFs)){
        for(int j=0 ; j<3 ; ++j){
          (*residual)[i+j] = beta*dt2*( (*density)[i/(3+numMultiphysDoFs)] *
                                        (*a)[i/(3+numMultiphysDoFs)*3+j] -
                                        (*force)[i/(3+numMultiphysDoFs)*3+j] -
                                        (*externalForce)[i/(3+numMultiphysDoFs)*3+j]);
        }
        //      Expression based off of the backward Euler integration scheme as well
        //      as Equation (34) from "A state-based peridynamic formulation of diffusive mass
        //      transport" (Amit Katiyar, John T. Foster, and Mukul Sharma).
        //
        //	The residual measures how well the previous estimate of fluidPressureU_n+1
        //	predicted the current estimate of fluidPressureU_n+1. If the difference is
        //	minimal, little change change be expected with further iteration.
        (*residual)[i+3] = ((*fluidPressureU)[i/(3+numMultiphysDoFs)]-(*fluidPressureUn)[i/(3+numMultiphysDoFs)])/dt -
                            (*fluidFlow)[i/(3+numMultiphysDoFs)]/((*fluidDensity)[i/(3+numMultiphysDoFs)]*(*fluidCompressibility)[i/(3+numMultiphysDoFs)]);
      }
    }
    else{
	for(int i=0 ; i<residual->MyLength() ; ++i)
   	     (*residual)[i] = beta*dt2*( (*density)[i/3] * (*a)[i] - (*force)[i] - (*externalForce)[i]);
    }

    // Modify residual for kinematic BC
    boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);

    double residualNorm;
    residual->Norm2(&residualNorm);

    int NLSolverIteration = 0;
    while(residualNorm > absoluteTolerance && NLSolverIteration <= maxSolverIterations){

      // Track the total number of iterations taken over the simulation
      *nonlinearSolverIterations += 1;

      if(peridigmComm->MyPID() == 0)
        cout << "  iteration " << NLSolverIteration << ": residual = " << residualNorm << endl;

      // Fill the Jacobian
      computeImplicitJacobian(beta, dt);

      // Modify Jacobian for kinematic BC
      boundaryAndInitialConditionManager->applyKinematicBC_InsertZerosAndSetDiagonal(tangent);

      // Want to solve J*displacementIncrement = -residual
      residual->Scale(-1.0);

      // Solve linear system
      displacementIncrement->PutScalar(0.0);
      if(analysisHasMultiphysics){
	fluidPressureDeltaU->PutScalar(0.0);
      }
      linearProblem.setOperator(tangent);

      bool isSet = linearProblem.setProblem(displacementIncrement, residual);

      TEUCHOS_TEST_FOR_EXCEPT_MSG(!isSet, "**** Peridigm::executeImplicit(), failed to set linear problem.\n");
      PeridigmNS::Timer::self().startTimer("Solve Linear System");
      Belos::ReturnType isConverged = belosSolver->solve();
      if(isConverged != Belos::Converged && peridigmComm->MyPID() == 0)
        cout << "Warning:  Belos linear solver failed to converge!  Proceeding with nonconverged solution..." << endl;
      PeridigmNS::Timer::self().stopTimer("Solve Linear System");

     //TODO: Turn all mentions of unknownsDeltaU and deltaU to displacementIncrement where appropriate.
      //Update increments from combined vector
      if(analysisHasMultiphysics){
	for(int i=0 ; i<displacementIncrement->MyLength() ; i+=(3+numMultiphysDoFs)){
	  (*fluidPressureDeltaU)[i/(3+numMultiphysDoFs)] = (*displacementIncrement)[i+3];
	}

	//Apply increment to fluidPressure displacement and current value analogues.
	fluidPressureU->Update(1.0,*fluidPressureDeltaU,1.0);
	fluidPressureY->Update(1.0, *fluidPressureU, 0.0);
      }
      // Apply increment in u to u
      for(int i=0 ; i<u->MyLength() ; i+=3)
	for(int j=0 ; j<3 ; ++j)
          (*u)[i+j] += (*displacementIncrement)[(3+numMultiphysDoFs)*i/3+j];

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
        blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(a, accelerationFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
        blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);

        if(analysisHasMultiphysics){
          blockIt->importData(fluidPressureU, fluidPressureUFieldId, PeridigmField::STEP_NP1, Insert);
          blockIt->importData(fluidPressureY, fluidPressureYFieldId, PeridigmField::STEP_NP1, Insert);
        }
      }
      PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

      // Update forces based on new positions
      PeridigmNS::Timer::self().startTimer("Internal Force");
      modelEvaluator->evalModel(workset);
      PeridigmNS::Timer::self().stopTimer("Internal Force");

      // Copy force from the data manager to the mothership vector
      force->PutScalar(0.0);
      if(analysisHasMultiphysics){
        fluidFlow->PutScalar(0.0);
        for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
          scratch->PutScalar(0.0);
          scalarScratch->PutScalar(0.0);
          blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
          blockIt->exportData(scalarScratch, fluidFlowDensityFieldId, PeridigmField::STEP_NP1, Add);
          force->Update(1.0, *scratch, 1.0);
          fluidFlow->Update(1.0, *scalarScratch, 1.0);
        }
        scalarScratch->PutScalar(0.0);
      }
      else{
        for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
          scratch->PutScalar(0.0);
          blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
          force->Update(1.0, *scratch, 1.0);
        }
      }
      scratch->PutScalar(0.0);

      // evaluate the external (body) forces:
      PeridigmNS::Timer::self().startTimer("Apply Body Forces");
      boundaryAndInitialConditionManager->applyForceContributions(timeCurrent, timePrevious);
      PeridigmNS::Timer::self().stopTimer("Apply Body Forces");

      // Compute residual vector and its norm
      // residual = beta*dt*dt*(M*a - force)
      // Note that due to restrictions to CrsMatrix, the residual has a different (but equivalent) map
      // than the force and acceleration
      if(analysisHasMultiphysics){
        for(int i=0 ; i<residual->MyLength() ; i+=(3+numMultiphysDoFs)){
          for(int j=0 ; j<3 ; ++j){
            (*residual)[i+j] = beta*dt2*( (*density)[i/(3+numMultiphysDoFs)] *
                                          (*a)[i/(3+numMultiphysDoFs)*3+j] -
                                          (*force)[i/(3+numMultiphysDoFs)*3+j] -
                                          (*externalForce)[i/(3+numMultiphysDoFs)*3+j]);
          }
          //TODO: Values for fluid density and compressibility other than 1.0 and 1.0
          //      need to be added, otherwise this expression is the transient diffusion
          //      residual. It is based off of the backward Euler integration scheme as well
          //      as Equation (34) from "A state-based peridynamic formulation of diffusive mass
          //      transport" (Amit Katiyar, John T. Foster, and Mukul Sharma).
          (*residual)[i+3] = (((*fluidPressureU)[i/(3+numMultiphysDoFs)]-
                               (*fluidPressureUn)[i/(3+numMultiphysDoFs)])/dt -
                              (*fluidFlow)[i/(3+numMultiphysDoFs)]);
        }
      }
      else{
        for(int i=0 ; i<residual->MyLength() ; ++i)
          (*residual)[i] = beta*dt2*( (*density)[i/3] * (*a)[i] - (*force)[i] - (*externalForce)[i]);
      }

      // Modify residual for kinematic BC
      boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);

      residual->Norm2(&residualNorm);

      NLSolverIteration++;
    }

    if(peridigmComm->MyPID() == 0)
      cout << "  iteration " << NLSolverIteration << ": residual = " << residualNorm << "\n" << endl;

    timePrevious = timeCurrent;
    timeCurrent = timeInitial + (step+1)*dt;

    // Write output for completed time step
    PeridigmNS::Timer::self().startTimer("Output");
    synchDataManagers();
    outputManager->write(blocks, timeCurrent);
    PeridigmNS::Timer::self().stopTimer("Output");

    // swap state N and state NP1
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->updateState();
  }
}

void PeridigmNS::Peridigm::allocateJacobian() {

  // do not re-allocate if already allocated
  if (tangent != Teuchos::null) return;

  int numDofs = PeridigmNS::DegreesOfFreedomManager::self().totalNumberOfDegreesOfFreedom();

  // Construct map for global tangent matrix
  // Note that this must be an Epetra_Map, not an Epetra_BlockMap, so we can't use threeDimensionalMap directly
  int numGlobalElements = numDofs * oneDimensionalMap->NumGlobalElements();
  int numMyElements = numDofs * oneDimensionalMap->NumMyElements();
  vector<int> myGlobalElements(numMyElements);
  int* oneDimensionalMapGlobalElements = oneDimensionalMap->MyGlobalElements();
  for(int iElem=0 ; iElem<oneDimensionalMap->NumMyElements() ; ++iElem){
      for(int dof = 0; dof < numDofs; dof++){
        myGlobalElements[(numDofs * iElem + dof)] = numDofs * oneDimensionalMapGlobalElements[iElem] + dof;
      }
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
  map<int, std::unordered_set<int> > rowEntries;

  // Loop over the neighborhood for each locally-owned point and record non-zero entries in the matrix.
  // Entries will exist for any two points that are bonded, and any two points that are bonded to a common third point.
  int* neighborhoodList = globalNeighborhoodData->NeighborhoodList();
  int neighborhoodListIndex = 0;
  vector<int> globalIndices;
  int numOwnedPoints = globalNeighborhoodData->NumOwnedPoints();
  for(int LID=0 ; LID<numOwnedPoints ; ++LID){
    int GID = oneDimensionalOverlapMap->GID(LID);
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];

    unsigned int numEntries = numDofs*(numNeighbors+1);

    if(globalIndices.size() < numEntries)
      globalIndices.resize(numEntries);

    for(int dof = 0; dof < numDofs; dof++){
        globalIndices[dof] = numDofs * GID + dof;
    }

    for(int j=0 ; j<numNeighbors ; ++j){
      int neighborLocalID = neighborhoodList[neighborhoodListIndex++];
      int neighborGlobalID = oneDimensionalOverlapMap->GID(neighborLocalID);

      for(int dof = 0; dof < numDofs; dof++){
          globalIndices[numDofs * j + numDofs + dof] = numDofs * neighborGlobalID + dof;
      }
    }

    // The entries going into the tangent are a dense matrix of size numEntries by numEntries.
    // Each global ID in the list interacts with all other global IDs in the list.
    for(unsigned int i=0 ; i<numEntries ; ++i){
      for(unsigned int j=0 ; j<numEntries ; ++j){
        rowEntries[globalIndices[i]].insert(globalIndices[j]);
      }
    }
  }

  // Allocate space in the global Epetra_FECrsMatrix
  vector<int> indices;
  vector<double> zeros;
  for(map<int, std::unordered_set<int> >::iterator rowEntry=rowEntries.begin(); rowEntry!=rowEntries.end() ; ++rowEntry){
    unsigned int numRowNonzeros = rowEntry->second.size();
    if(zeros.size() < numRowNonzeros)
      zeros.resize(numRowNonzeros, 0.0);

    // Load indices into a sorted vector
    indices.resize(numRowNonzeros);
    int i=0;
    for(std::unordered_set<int>::const_iterator globalIndex=rowEntry->second.begin() ; globalIndex!=rowEntry->second.end() ; ++globalIndex)
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

  if(tangent != Teuchos::null) // the tangent matrix is already allocated, use it for the block diagonal as well
  {
    blockDiagonalTangent = tangent;
    blockDiagonalTangentMap = tangentMap;
    return;
  }

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
  int numEntriesPerRow = 3;
  bool ignoreNonLocalEntries = false;
  blockDiagonalTangent = Teuchos::rcp(new Epetra_FECrsMatrix(CV, *blockDiagonalTangentMap, numEntriesPerRow, ignoreNonLocalEntries));

  // Store nonzero columns for each row, with everything in global indices
  // The matrix is block 3x3, very straightforward
  int err, rowEntries[3];
  const double zeros[3] = {0.0, 0.0, 0.0};
  for(int row=0 ; row<blockDiagonalTangentMap->NumMyElements() ; row++){
    int globalId = blockDiagonalTangentMap->GID(row);
    rowEntries[0] = 3*(static_cast<int>(globalId)/3);
    rowEntries[1] = 3*(static_cast<int>(globalId)/3) + 1;
    rowEntries[2] = 3*(static_cast<int>(globalId)/3) + 2;
    err = blockDiagonalTangent->InsertGlobalValues(globalId, numEntriesPerRow, zeros, (const int*)rowEntries);
    TEUCHOS_TEST_FOR_EXCEPT_MSG(err < 0, "**** PeridigmNS::Peridigm::allocateblockDiagonalJacobian(), InsertGlobalValues() returned negative error code.\n");
  }
  err = blockDiagonalTangent->GlobalAssemble();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::allocateBlockDiagonalJacobian(), GlobalAssemble() returned nonzero error code.\n");

  // create the serial Jacobian
  overlapJacobian = Teuchos::rcp(new PeridigmNS::SerialMatrix(blockDiagonalTangent));
  workset->jacobian = overlapJacobian;

  PeridigmNS::Memstat * memstat = PeridigmNS::Memstat::Instance();
  const std::string statTag = "Block Diagonal Tangent";
  memstat->addStat(statTag);
}

double PeridigmNS::Peridigm::computeQuasiStaticResidual(Teuchos::RCP<Epetra_Vector> residual) {

  PeridigmNS::Timer::self().startTimer("Compute Residual");

  PeridigmNS::DegreesOfFreedomManager& dofManager = PeridigmNS::DegreesOfFreedomManager::self();
  int numDofs = dofManager.totalNumberOfDegreesOfFreedom();
  bool displacementTreatedAsUnknown = dofManager.displacementTreatedAsUnknown();
  int numDisplacementDofs = dofManager.numberOfDisplacementDegreesOfFreedom();
  int displacementDofOffset = dofManager.displacementDofOffset();
  bool temperatureTreatedAsUnknown = dofManager.temperatureTreatedAsUnknown();
  int temperatureDofOffset = dofManager.temperatureDofOffset();
  bool concentrationTreatedAsUnknown = dofManager.concentrationTreatedAsUnknown();
  int concentrationDofOffset = dofManager.concentrationDofOffset();
  bool pressureTreatedAsUnknown = dofManager.pressureTreatedAsUnknown();
  int pressureDofOffset = dofManager.pressureDofOffset();

  // Copy data from mothership vectors to overlap vectors in data manager
  PeridigmNS::Timer::self().startTimer("Gather/Scatter");
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(fluidPressureU, fluidPressureUFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(fluidPressureY, fluidPressureYFieldId, PeridigmField::STEP_NP1, Insert);
    blockIt->importData(fluidPressureV, fluidPressureVFieldId, PeridigmField::STEP_NP1, Insert);
  }
  PeridigmNS::Timer::self().stopTimer("Gather/Scatter");

  PeridigmNS::Timer::self().startTimer("Internal Force");
  modelEvaluator->evalModel(workset);
  PeridigmNS::Timer::self().stopTimer("Internal Force");

  // Copy force from the data manager to the mothership vector
  double* residual_ptr;
  residual->ExtractView(&residual_ptr);
  if(displacementTreatedAsUnknown) {
    force->PutScalar(0.0);
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scratch->PutScalar(0.0);
      blockIt->exportData(scratch, forceDensityFieldId, PeridigmField::STEP_NP1, Add);
      force->Update(1.0, *scratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
    double* force_ptr;
    force->ExtractView(&force_ptr);
    for(int i_node = 0; i_node < force->Map().NumMyElements(); ++i_node) {
      for (int dof = 0; dof < numDisplacementDofs; ++dof) {
        double value = force_ptr[i_node*numDisplacementDofs + dof];
        TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN returned by force evaluation.\n");
        residual_ptr[i_node*numDofs + displacementDofOffset + dof] = value;
      }
    }
  }
  if(temperatureTreatedAsUnknown) {
    fluxDivergence->PutScalar(0.0);
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scalarScratch->PutScalar(0.0);
      blockIt->exportData(scalarScratch, fluxDivergenceFieldId, PeridigmField::STEP_NP1, Add);
      fluxDivergence->Update(1.0, *scalarScratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
    double* flux_divergence_ptr;
    fluxDivergence->ExtractView(&flux_divergence_ptr);
    for(int i_node = 0; i_node < fluxDivergence->Map().NumMyElements(); ++i_node) {
      double value = flux_divergence_ptr[i_node];
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN returned by flux divergence evaluation.\n");
      residual_ptr[i_node*numDofs + temperatureDofOffset] = value;
    }
  }
  if(concentrationTreatedAsUnknown) {
    concentrationFluxDivergence->PutScalar(0.0);
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
      scalarScratch->PutScalar(0.0);
      blockIt->exportData(scalarScratch, concentrationFluxDivergenceFieldId, PeridigmField::STEP_NP1, Add);
      concentrationFluxDivergence->Update(1.0, *scalarScratch, 1.0);
    }
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
    double* concentration_flux_divergence_ptr;
    concentrationFluxDivergence->ExtractView(&concentration_flux_divergence_ptr);
    for(int i_node = 0; i_node < concentrationFluxDivergence->Map().NumMyElements(); ++i_node) {
      double value = concentration_flux_divergence_ptr[i_node];
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN returned by concentration flux divergence evaluation.\n");
      residual_ptr[i_node*numDofs + concentrationDofOffset] = value;
    }
  }
  if(pressureTreatedAsUnknown) {
    fluidFlow->PutScalar(0.0);
    PeridigmNS::Timer::self().startTimer("Gather/Scatter");
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
			scalarScratch->PutScalar(0.0);
			blockIt->exportData(scalarScratch, fluidFlowDensityFieldId, PeridigmField::STEP_NP1, Add);
			fluidFlow->Update(1.0, *scalarScratch, 1.0);
		}
    PeridigmNS::Timer::self().stopTimer("Gather/Scatter");
    double* fluid_flow_ptr;
    fluidFlow->ExtractView(&fluid_flow_ptr);
    for(int i_node = 0; i_node < fluidFlow->Map().NumMyElements(); ++i_node) {
      double value = fluid_flow_ptr[i_node];
      TEUCHOS_TEST_FOR_EXCEPT_MSG(!std::isfinite(value), "**** NaN returned by fluid flow evaluation.\n");
      residual_ptr[i_node*numDofs + pressureDofOffset] = value;
    }
  }

  // Add in the external force
  if(displacementTreatedAsUnknown) {
    double* external_force_ptr;
    externalForce->ExtractView(&external_force_ptr);
    for(int i_node = 0; i_node < externalForce->Map().NumMyElements(); ++i_node) {
      for (int dof = 0; dof < numDisplacementDofs; ++dof) {
        residual_ptr[i_node*numDofs + displacementDofOffset + dof] += external_force_ptr[i_node*numDisplacementDofs + dof];
      }
    }
  }

  // multiply the residual by volume (i.e., convert force density to force)
  double* volume_ptr;
  volume->ExtractView(&volume_ptr);
	for(int i_node=0 ; i_node<volume->Map().NumMyElements() ; ++i_node) {
    double v = volume_ptr[i_node];
    for (int dof = 0; dof < numDofs; ++dof) {
      residual_ptr[i_node*numDofs + dof] *= v;
    }
  }

  // zero out the rows corresponding to kinematic boundary conditions and compute the residual
  boundaryAndInitialConditionManager->applyKinematicBC_InsertZeros(residual);

  double residualNorm2;
  residual->Norm2(&residualNorm2);
  double residualNormInf;
  residual->NormInf(&residualNormInf);

  PeridigmNS::Timer::self().stopTimer("Compute Residual");

  return residualNorm2 + 20.0*residualNormInf;
}

void PeridigmNS::Peridigm::computeImplicitJacobian(double beta, double dt) {
//TODO make multiphysics
  // Compute the tangent
  tangent->PutScalar(0.0);
  PeridigmNS::Timer::self().startTimer("Evaluate Jacobian");
  modelEvaluator->evalJacobian(workset);
  int err = tangent->GlobalAssemble();
  TEUCHOS_TEST_FOR_EXCEPT_MSG(err != 0, "**** PeridigmNS::Peridigm::computeImplicitJacobian(), GlobalAssemble() returned nonzero error code.\n");
  PeridigmNS::Timer::self().stopTimer("Evaluate Jacobian");

  // tangent = M - beta*dt*dt*K
  tangent->Scale(-beta*dt*dt);

	//TODO: esp this part
  Epetra_Vector diagonal(tangent->RowMap());
  tangent->ExtractDiagonalCopy(diagonal);
  for(int i=0 ; i<diagonal.MyLength() ; ++i)
    diagonal[i] += (*density)[i/3];
  tangent->ReplaceDiagonalValues(diagonal);
}

void PeridigmNS::Peridigm::synchDataManagers() {

  // Copy data from mothership vectors to overlap vectors in blocks
  // Volume, Block_Id, and Model_Coordinates are synched at initialization and never change

  PeridigmNS::Timer::self().startTimer("Gather/Scatter");

	if(analysisHasMultiphysics){
		for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
			blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(force, forceDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(fluidFlow, fluidFlowDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(contactForce, contactForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(externalForce, externalForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(fluidPressureU, fluidPressureUFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(fluidPressureY, fluidPressureYFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(fluidPressureV, fluidPressureVFieldId, PeridigmField::STEP_NP1, Insert);
    }
	}
	else{
		for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
			blockIt->importData(u, displacementFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(y, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(v, velocityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(force, forceDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(contactForce, contactForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(externalForce, externalForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
			blockIt->importData(temperature, temperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(deltaTemperature, deltaTemperatureFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(concentration, concentrationFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(fluxDivergence, fluxDivergenceFieldId, PeridigmField::STEP_NP1, Insert);
      blockIt->importData(concentrationFluxDivergence, concentrationFluxDivergenceFieldId, PeridigmField::STEP_NP1, Insert);
		}
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
      blockIt->exportData(scratch, hourglassForceDensityFieldId, PeridigmField::STEP_NP1, Add);
      tempVector->Update(1.0, *scratch, 1.0);
    }
    for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
      blockIt->importData(tempVector, hourglassForceDensityFieldId, PeridigmField::STEP_NP1, Insert);
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

void PeridigmNS::Peridigm::writeRestart(Teuchos::RCP<Teuchos::ParameterList> solverParams){
//  system("date +"%m-%d-%Y-%H-%M-%S"");
  char createDirectory[100];
  char  path[100];
  int IterationNumber;

  if(peridigmComm->MyPID() == 0){
  IterationNumber = atoi(firstNumbersSring( restartFiles["path"]  ).c_str())+1;
  sprintf(path,"restart-%06d",IterationNumber);
  cout << "The restart folder is " << path  <<"." << endl;
  setRestartNames(path);
  sprintf(createDirectory,"mkdir %s",path);
  system(createDirectory);
  cout << "Writing restart files. \n" << endl;

  double timeInitial = solverParams->get("Initial Time", 0.0);
  if(peridigmComm->MyPID() == 0)
	  if (currentTime != timeInitial){
		  char timeError[251];
		  sprintf(timeError, "Error, Incompatible times:\nPrevious restart final time is %e, while initial time is %e.\n",currentTime,timeInitial);
		  TEUCHOS_TEST_FOR_EXCEPT_MSG(true,timeError);
		  MPI_Finalize();
		  exit(0);
	  }
  currentTime += solverParams->get("Final Time", 1.0)-timeInitial;
  ofstream outputFile;
  outputFile.open(restartFiles["currentTime"].c_str());
  outputFile << "Current time is " << "\n" << currentTime  << "\n";
  outputFile.close();
  }
  if(analysisHasMultiphysics){
	 cout << "Restart for Multiphysics is not implemented yet." << endl;
	 exit (0);
    }
  else {
	  //write block ID
	  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["blockIDs"].c_str(),
	  		*blockIDs,"blockIDs","",true);
      //write horizon for each point
	  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["horizon"].c_str(),
	  		*horizon,"horizon","",true);
	  //write cell volume
	  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["volume"].c_str(),
	  		*volume,"volume","",true);
	  //write density
	  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["density"].c_str(),
	  		*density,"density","",true);
	  //write change in temperature
	  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["deltaTemperature"].c_str(),
	  		*deltaTemperature,"deltaTemperature","",true);
  }
  //write initial positions
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["x"].c_str(),
  		*x,"x","",true);
  //write displacement
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["u"].c_str(),
  		*u,"u","",true);
  //write current positions
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["y"].c_str(),
  		*y,"y","",true);
  //write velocities
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["v"].c_str(),
  		*v,"v","",true);
  //write accelerations
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["a"].c_str(),
  		*a,"a","",true);
  //write force
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["force"].c_str(),
  		*force,"force","",true);
  //write contact force
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["contactForce"].c_str(),
  		*contactForce,"contactForce","",true);
  //write external force
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["externalForce"].c_str(),
  		*externalForce,"externalForce","",true);
  //write deltaU (increment in displacement)
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["deltaU"].c_str(),
  		*deltaU,"deltaU","",true);
  //write scratch
  EpetraExt::VectorToMatrixMarketFile 	(restartFiles["scratch"].c_str(),
  		*scratch,"scratch","",true);
  //write block data
	  std::vector<PeridigmNS::Block>::iterator blockIt;
	  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
     	  std::string blockName = blockIt->getName();
		  blockIt->writeBlocktoDisk(blockName,restartFiles["path"].c_str());
	  }
}

void PeridigmNS::Peridigm::readRestart(){
	  double* UpdatePtr;
	  double* oldPtr;
	  Epetra_Vector * vectorUpdate;
	  std::string trash, data;
	  if(peridigmComm->MyPID() == 0){
		  //read global current time
		  ifstream outputFile;
		  outputFile.open(restartFiles["currentTime"].c_str());
	      	  if (!outputFile.fail())
	      	  {
	      		  getline(outputFile,trash);
	      		  getline(outputFile,data);
	      		  outputFile.close();
	      	  }
	      currentTime = atof(data.c_str());
	  }
  if(peridigmComm->MyPID() == 0){
  	cout <<"Reading restart. \n"<< endl;
  	cout.flush();
  }
  if(analysisHasMultiphysics){
	  if(peridigmComm->MyPID() == 0){
		  TEUCHOS_TEST_FOR_EXCEPT_MSG(true,"Error: Restart for Multiphysics is not implemented yet.\n");
		  MPI_Finalize();
		  exit(0);
	  }
  }else{
	  //read block ID
	  EpetraExt::MatrixMarketFileToVector(restartFiles["blockIDs"].c_str(), *oneDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  blockIDs->ExtractView(&oldPtr);
	  blas.COPY(blockIDs->MyLength(), UpdatePtr, oldPtr);
      //read horizon for each point
	  EpetraExt::MatrixMarketFileToVector(restartFiles["horizon"].c_str(), *oneDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  horizon->ExtractView(&oldPtr);
	  blas.COPY(horizon->MyLength(), UpdatePtr, oldPtr);
	  //read cell volume
	  EpetraExt::MatrixMarketFileToVector(restartFiles["volume"].c_str(), *oneDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  volume->ExtractView(&oldPtr);
	  blas.COPY(volume->MyLength(), UpdatePtr, oldPtr);
	  //read density
	  EpetraExt::MatrixMarketFileToVector(restartFiles["density"].c_str(), *oneDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  density->ExtractView(&oldPtr);
	  blas.COPY(density->MyLength(), UpdatePtr, oldPtr);
	  //read change in temperature
	  EpetraExt::MatrixMarketFileToVector(restartFiles["deltaTemperature"].c_str(), *oneDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  deltaTemperature->ExtractView(&oldPtr);
	  blas.COPY(deltaTemperature->MyLength(), UpdatePtr, oldPtr);
  }
	  //read initial positions
	  EpetraExt::MatrixMarketFileToVector(restartFiles["x"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  x->ExtractView(&oldPtr);
	  blas.COPY(x->MyLength(), UpdatePtr, oldPtr);
	  //read displacement
	  EpetraExt::MatrixMarketFileToVector(restartFiles["u"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  u->ExtractView(&oldPtr);
	  blas.COPY(u->MyLength(), UpdatePtr, oldPtr);
	  //read current positions
	  EpetraExt::MatrixMarketFileToVector(restartFiles["y"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  y->ExtractView(&oldPtr);
	  blas.COPY(y->MyLength(), UpdatePtr, oldPtr);
	  //read velocities
	  EpetraExt::MatrixMarketFileToVector(restartFiles["v"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  v->ExtractView(&oldPtr);
	  blas.COPY(v->MyLength(), UpdatePtr, oldPtr);
	  //read accelerations
	  EpetraExt::MatrixMarketFileToVector(restartFiles["a"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  a->ExtractView(&oldPtr);
	  blas.COPY(a->MyLength(), UpdatePtr, oldPtr);
	  //read force
	  EpetraExt::MatrixMarketFileToVector(restartFiles["force"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  force->ExtractView(&oldPtr);
	  blas.COPY(force->MyLength(), UpdatePtr, oldPtr);
	  //read contact force
	  EpetraExt::MatrixMarketFileToVector(restartFiles["contactForce"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  contactForce->ExtractView(&oldPtr);
	  blas.COPY(contactForce->MyLength(), UpdatePtr, oldPtr);
	  vectorUpdate = 0;
	  UpdatePtr = 0;
	  oldPtr = 0;
	  //read external force
	  EpetraExt::MatrixMarketFileToVector(restartFiles["externalForce"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  externalForce->ExtractView(&oldPtr);
	  blas.COPY(externalForce->MyLength(), UpdatePtr, oldPtr);
	  //read deltaU (increment in displacement)
	  EpetraExt::MatrixMarketFileToVector(restartFiles["deltaU"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  deltaU->ExtractView(&oldPtr);
	  blas.COPY(deltaU->MyLength(), UpdatePtr, oldPtr);
	  //read scratch
	  EpetraExt::MatrixMarketFileToVector(restartFiles["scratch"].c_str(), *threeDimensionalMap, vectorUpdate);
	  vectorUpdate->ExtractView(&UpdatePtr);
	  scratch->ExtractView(&oldPtr);
	  blas.COPY(scratch->MyLength(), UpdatePtr, oldPtr);
	  //read block data
	  	  std::vector<PeridigmNS::Block>::iterator blockIt;
	  	  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
	     	  std::string blockName = blockIt->getName();
			  blockIt->readBlockfromDisk(blockName,restartFiles["path"].c_str());
	  	  }
}
