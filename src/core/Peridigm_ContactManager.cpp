/*! \file Peridigm_ContactManager.cpp */

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

#include "Peridigm_ContactManager.hpp"
#include "contact/Peridigm_ContactModelFactory.hpp"
#include "Peridigm_Timer.hpp"
#include <boost/algorithm/string/trim.hpp> // \todo Replace this include with correct include for istream_iterator.
#include "Peridigm_PdQuickGridDiscretization.hpp"

#include "pdneigh/PdZoltan.h"
#include "pdneigh/NeighborhoodList.h"

using namespace std;

//#define NEW_STUFF 1

PeridigmNS::ContactManager::ContactManager(const Teuchos::ParameterList& contactParams,
                                           Teuchos::RCP<Discretization> disc,
                                           Teuchos::RCP<Teuchos::ParameterList> peridigmParams)
  : verbose(false), myPID(-1), params(contactParams), contactRebalanceFrequency(0), contactSearchRadius(0.0),
    blockIdFieldId(-1), volumeFieldId(-1), coordinatesFieldId(-1), velocityFieldId(-1), contactForceDensityFieldId(-1)
{
  if(contactParams.isParameter("Verbose"))
    verbose = contactParams.get<bool>("Verbose");
  myPID = disc->getGlobalOwnedMap(1)->Comm().MyPID();

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  blockIdFieldId = fieldManager.getFieldId("Block_Id");
  volumeFieldId = fieldManager.getFieldId("Volume");
  coordinatesFieldId = fieldManager.getFieldId("Coordinates");
  velocityFieldId = fieldManager.getFieldId("Velocity");
  contactForceDensityFieldId = fieldManager.getFieldId("Contact_Force_Density");

  if(!contactParams.isParameter("Search Radius"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Radius\" not specified.");
  contactSearchRadius = contactParams.get<double>("Search Radius");
  if(!contactParams.isParameter("Search Frequency"))
    TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, "Contact parameter \"Search Frequency\" not specified.");
  contactRebalanceFrequency = contactParams.get<int>("Search Frequency");

  createContactInteractionsList(contactParams, disc);

  


  // Did user specify default blocks?
  bool defaultBlocks = false;
  // Parameterlist for default blocks
  Teuchos::ParameterList defaultBlockParams;

  // Create vector of blocks
  contactBlocks = Teuchos::rcp(new vector<PeridigmNS::ContactBlock>());

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

      PeridigmNS::ContactBlock contactBlock(*it, blockID, params);
      contactBlocks->push_back(contactBlock);
    }
  }

  // Add in all default blocks
  if (defaultBlocks) {
    vector<string> discretizationBlockNames = disc->getBlockNames();
    for(vector<string>::const_iterator it=discretizationBlockNames.begin() ; it!=discretizationBlockNames.end() ; ++it){
      bool blockMatch = false;
      for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
        // if name match, break
        if ((*it) == contactBlockIt->getName()) {
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

        PeridigmNS::ContactBlock contactBlock(*it, blockID, defaultBlockParams);
        contactBlocks->push_back(contactBlock);
      }
    }
  }
}

void PeridigmNS::ContactManager::createContactInteractionsList(const Teuchos::ParameterList& contactParams,
                                                               Teuchos::RCP<Discretization> disc)
{
  set< boost::tuple<int, int, string> > contactInteractionSet;
  Teuchos::ParameterList interactionParams = contactParams.sublist("Interactions");

  // First, check for general and self contact
  bool hasGeneralContact = interactionParams.isSublist("General Contact");
  string generalContactModelName, selfContactModelName;
  if(hasGeneralContact)
    generalContactModelName = interactionParams.sublist("General Contact").get<string>("Contact Model");
  bool hasSelfContact = interactionParams.isSublist("Self Contact");
  if(hasSelfContact)
    selfContactModelName = interactionParams.sublist("Self Contact").get<string>("Contact Model");
  if(hasGeneralContact || hasSelfContact){
    vector<string> blockNames = disc->getBlockNames() ;
    for(unsigned int i=0 ; i<blockNames.size() ; ++i){
      for(unsigned int j=0 ; j<blockNames.size() ; ++j){
        int blockId_1 = disc->blockNameToBlockId(blockNames[i]);
        int blockId_2 = disc->blockNameToBlockId(blockNames[j]);
        if(blockId_1 > blockId_2){
          int temp = blockId_1;
          blockId_1 = blockId_2;
          blockId_2 = temp;
        }
        if(hasSelfContact && i == j){
          boost::tuple<int, int, string> interactionTuple(blockId_1, blockId_2, selfContactModelName);
          contactInteractionSet.insert(interactionTuple);          
        }
        if(hasGeneralContact && i != j){
          boost::tuple<int, int, string> interactionTuple(blockId_1, blockId_2, generalContactModelName);
          contactInteractionSet.insert(interactionTuple);
        }
      }
    }
  }

  // Parse the individually-defined interactions
  for(Teuchos::ParameterList::ConstIterator it = interactionParams.begin() ; it != interactionParams.end() ; it++){
    const string& name = it->first;
    if(name != "General Contact" && name != "Self Contact" && interactionParams.isSublist(name)){
      Teuchos::ParameterList& interaction = interactionParams.sublist(name);
      int blockId_1 = disc->blockNameToBlockId( interaction.get<string>("First Block") );
      int blockId_2 = disc->blockNameToBlockId( interaction.get<string>("Second Block") );
      string contactModelName = interaction.get<string>("Contact Model");
      boost::tuple<int, int, string> interactionTuple(blockId_1, blockId_2, contactModelName);
      contactInteractionSet.insert(interactionTuple);
    }
  }

  // Elliminate duplicates and load the interactions into a vector
  set< boost::tuple<int, int, string> >::iterator it, nextIt;
  nextIt = contactInteractionSet.begin();
  nextIt++;
  for(it = contactInteractionSet.begin() ; it != contactInteractionSet.end() ; it++){
    // If a block pair has multiple contact definitions, the final definition is the one that will be used
    bool duplicate = false;
    if( nextIt != contactInteractionSet.end() && it->get<0>() == nextIt->get<0>() && it->get<1>() == nextIt->get<1>() )
      duplicate = true;
    if(!duplicate)
      contactInteractions.push_back(*it);
    if( nextIt != contactInteractionSet.end() )
      nextIt++;
  }

  // Print to stdout if verbose flag is set
  if(verbose && myPID == 0){
    cout << "--Contact Interactions--" << endl;
    cout << "  First block    Second block    Contact model" << endl;
    for(unsigned int i=0 ; i<contactInteractions.size() ; ++i)
      cout << "  " << contactInteractions[i].get<0>() << "              " <<  contactInteractions[i].get<1>() << "               " <<  contactInteractions[i].get<2>() << endl;
    cout << endl;
  }
}

void PeridigmNS::ContactManager::initialize(Teuchos::RCP<const Epetra_BlockMap> oneDimensionalMap_,
                                            Teuchos::RCP<const Epetra_BlockMap> threeDimensionalMap_,
                                            Teuchos::RCP<const Epetra_BlockMap> oneDimensionalOverlapMap_,
                                            Teuchos::RCP<const Epetra_BlockMap> bondMap_,
                                            Teuchos::RCP<const Epetra_Vector> blockIds_,
                                            map<string, double> blockHorizonValues)
{
  ContactModelFactory contactModelFactory;
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){

    // Obtain the horizon for this block
    string blockName = contactBlockIt->getName();
    double blockHorizon(0.0);
    if(blockHorizonValues.find(blockName) != blockHorizonValues.end())
      blockHorizon = blockHorizonValues[blockName];
    else if(blockHorizonValues.find("default") != blockHorizonValues.end())
      blockHorizon = blockHorizonValues["default"];
    else{
      string msg = "\n**** Error, no Horizon parameter supplied for block " + blockName + " and no default block parameter list provided.\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
    }

    // For the initial implementation, assume that there is only one contact model
    // \todo Refactor contact!
    Teuchos::ParameterList contactModelParams = params.sublist("Models", true);
    Teuchos::ParameterList contactParams = contactModelParams.sublist( contactModelParams.begin()->first );
    TEUCHOS_TEST_FOR_EXCEPT_MSG(contactParams.isParameter("Horizon") , "\n**** Error, Horizon is an invalid contact model parameter.\n");
    contactParams.set("Horizon", blockHorizon);
    if(!contactParams.isParameter("Friction Coefficient"))
      contactParams.set("Friction Coefficient", 0.0);
    Teuchos::RCP<const PeridigmNS::ContactModel> contactModel = contactModelFactory.create(contactParams);
    contactBlockIt->setContactModel(contactModel);
  }

  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(*oneDimensionalMap_));
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(*threeDimensionalMap_));
  oneDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(*oneDimensionalOverlapMap_));
  bondMap = Teuchos::rcp(new Epetra_BlockMap(*bondMap_));

  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(oneDimensionalOverlapMap->NumGlobalElements(),
                                                                oneDimensionalOverlapMap->NumMyElements(),
                                                                oneDimensionalOverlapMap->MyGlobalElements(),
                                                                3,
                                                                0,
                                                                oneDimensionalOverlapMap->Comm()));

  PeridigmNS::Timer::self().startTimer("Initialize Contact Maps");

#ifdef NEW_STUFF

  set<int> blocksWithContact;
  for(vector< boost::tuple<int, int, string> >::iterator it=contactInteractions.begin() ; it!=contactInteractions.end() ; it++){
    blocksWithContact.insert( it->get<0>() );
    blocksWithContact.insert( it->get<1>() );
  }

  vector<int> contactGlobalIds;
  contactGlobalIds.reserve(oneDimensionalMap->NumMyElements());
  for(int i=0 ; i<oneDimensionalMap->NumMyElements() ; ++i){
    int blockId = (*blockIds_)[i];
    if(blocksWithContact.find(blockId) != blocksWithContact.end())
      contactGlobalIds.push_back(oneDimensionalMap->GID(i));
  }

  // Instantiate the owned maps for the contact mothership vectors
  oneDimensionalContactMap = Teuchos::rcp(new Epetra_BlockMap(-1, contactGlobalIds.size(), &contactGlobalIds[0], 1, 0, oneDimensionalMap_->Comm()));
  threeDimensionalContactMap = Teuchos::rcp(new Epetra_BlockMap(-1, contactGlobalIds.size(), &contactGlobalIds[0], 3, 0, oneDimensionalMap_->Comm()));

  // Create an overlap vector with block ids
  Epetra_Import tempImporter(*oneDimensionalOverlapMap, *oneDimensionalMap);
  Epetra_Vector overlapBlockIds(*oneDimensionalOverlapMap);
  overlapBlockIds.Import(*blockIds_, tempImporter, Insert);

  vector<int> contactGlobalOverlapIds;
  contactGlobalOverlapIds.reserve(oneDimensionalOverlapMap->NumMyElements());
  for(int i=0 ; i<oneDimensionalOverlapMap->NumMyElements() ; ++i){
    int blockId = overlapBlockIds[i];
    if(blocksWithContact.find(blockId) != blocksWithContact.end())
      contactGlobalOverlapIds.push_back(oneDimensionalOverlapMap->GID(i));
  }

  // Instantiate the overlap maps for the contact mothership vectors
  oneDimensionalOverlapContactMap = Teuchos::rcp(new Epetra_BlockMap(-1, contactGlobalOverlapIds.size(), &contactGlobalOverlapIds[0], 1, 0, oneDimensionalMap_->Comm()));
  threeDimensionalOverlapContactMap = Teuchos::rcp(new Epetra_BlockMap(-1, contactGlobalOverlapIds.size(), &contactGlobalOverlapIds[0], 3, 0, oneDimensionalMap_->Comm()));

  // Create the contact bond map
  // Care must be taken because elements with zero neighbors are not included in the bond map (Epetra_BlockMap does not support elements with zero nodes)
  vector<int> contactBondMapGlobalIds;
  contactBondMapGlobalIds.reserve(bondMap_->NumMyElements());
  vector<int> contactBondMapElementSizeList;
  contactBondMapElementSizeList.reserve(bondMap_->NumMyElements());
  int* bondMapElementSizes = bondMap_->ElementSizeList();
  for(int i=0 ; i<oneDimensionalContactMap->NumMyElements() ; ++i){
    int globalId = oneDimensionalContactMap->GID(i);
    int bondMapLocalId = bondMap_->LID(globalId);
    if(bondMapLocalId != -1){
      contactBondMapGlobalIds.push_back(globalId);
      contactBondMapElementSizeList.push_back( bondMapElementSizes[bondMapLocalId] );
    }
  }
  bondContactMap = Teuchos::rcp(new Epetra_BlockMap(-1, contactBondMapGlobalIds.size(), &contactBondMapGlobalIds[0], &contactBondMapElementSizeList[0], 0, oneDimensionalMap_->Comm()));

#else

  oneDimensionalContactMap = Teuchos::rcp(new Epetra_BlockMap(*oneDimensionalMap));
  threeDimensionalContactMap = Teuchos::rcp(new Epetra_BlockMap(*threeDimensionalMap));
  oneDimensionalOverlapContactMap = Teuchos::rcp(new Epetra_BlockMap(*oneDimensionalOverlapMap));
  threeDimensionalOverlapContactMap = Teuchos::rcp(new Epetra_BlockMap(-1, oneDimensionalOverlapContactMap->NumMyElements(), oneDimensionalOverlapContactMap->MyGlobalElements(), 3, 0, oneDimensionalMap_->Comm()));
  bondContactMap = Teuchos::rcp(new Epetra_BlockMap(*bondMap));

#endif

  PeridigmNS::Timer::self().stopTimer("Initialize Contact Maps");

  // Instantiate the importers for passing data between the mothership and contact mothership vectors
  oneDimensionalMothershipToContactMothershipImporter = Teuchos::rcp(new Epetra_Import(*oneDimensionalContactMap, *oneDimensionalMap));
  threeDimensionalMothershipToContactMothershipImporter = Teuchos::rcp(new Epetra_Import(*threeDimensionalContactMap, *threeDimensionalMap));

  // Create the contact mothership multivectors
  oneDimensionalContactMothership = Teuchos::rcp(new Epetra_MultiVector(*oneDimensionalContactMap, 2));
  contactBlockIDs = Teuchos::rcp((*oneDimensionalContactMothership)(0), false);         // block ID
  contactVolume = Teuchos::rcp((*oneDimensionalContactMothership)(1), false);           // cell volume

  threeDimensionalContactMothership = Teuchos::rcp(new Epetra_MultiVector(*threeDimensionalContactMap, 4));
  contactY = Teuchos::rcp((*threeDimensionalContactMothership)(0), false);             // current positions
  contactV = Teuchos::rcp((*threeDimensionalContactMothership)(1), false);             // velocities
  contactContactForce = Teuchos::rcp((*threeDimensionalContactMothership)(2), false);  // contact force
  contactScratch = Teuchos::rcp((*threeDimensionalContactMothership)(3), false);       // scratch
}

void PeridigmNS::ContactManager::loadAllMothershipData(Teuchos::RCP<Epetra_Vector> blockIds,
                                                       Teuchos::RCP<Epetra_Vector> volume,
                                                       Teuchos::RCP<Epetra_Vector> y,
                                                       Teuchos::RCP<Epetra_Vector> v)
{
  contactBlockIDs->Import(*blockIds, *oneDimensionalMothershipToContactMothershipImporter, Insert);
  contactVolume->Import(*volume, *oneDimensionalMothershipToContactMothershipImporter, Insert);
  contactY->Import(*y, *threeDimensionalMothershipToContactMothershipImporter, Insert);
  contactV->Import(*v, *threeDimensionalMothershipToContactMothershipImporter, Insert);
  contactContactForce->PutScalar(0.0);
  contactScratch->PutScalar(0.0);
}

void PeridigmNS::ContactManager::loadNeighborhoodData(Teuchos::RCP<PeridigmNS::NeighborhoodData> globalNeighborhoodData,
                                                      Teuchos::RCP<const Epetra_BlockMap> globalNeighborhoodDataOneDimensionalMap,
                                                      Teuchos::RCP<const Epetra_BlockMap> globalNeighborhoodDataOneDimensionalOverlapMap)
{
#ifdef NEW_STUFF

  int peridigmNumOwned = globalNeighborhoodData->NumOwnedPoints();
  int* peridigmOwnedIds = globalNeighborhoodData->OwnedIDs();
  int peridigmNeighborhoodListSize = globalNeighborhoodData->NeighborhoodListSize();
  int* peridigmNeighborhoodList = globalNeighborhoodData->NeighborhoodList();

  vector<int> ownedIds;
  vector<int> neighborhoodList;
  vector<int> neighborhoodPtr;
  ownedIds.reserve(peridigmNumOwned);
  neighborhoodList.reserve(peridigmNeighborhoodListSize);
  neighborhoodPtr.reserve(peridigmNumOwned);

  int peridigmNeighborhoodIndex = 0;
  for(int i=0 ; i<peridigmNumOwned; ++i){

    int peridigmNumNeighbors = peridigmNeighborhoodList[peridigmNeighborhoodIndex++];

    int globalId = globalNeighborhoodDataOneDimensionalMap->GID(peridigmOwnedIds[i]);
    int localId = oneDimensionalContactMap->LID(globalId);
    if(localId != -1){
      ownedIds.push_back(localId);
      vector<int> tempNeighborList;
      for(int j=0 ; j<peridigmNumNeighbors ; ++j){
        globalId = globalNeighborhoodDataOneDimensionalOverlapMap->GID(peridigmNeighborhoodList[peridigmNeighborhoodIndex + j]);
        localId = oneDimensionalOverlapContactMap->LID(globalId);
        if(localId != -1)
          tempNeighborList.push_back(localId);
      }
      neighborhoodPtr.push_back(neighborhoodList.size());
      neighborhoodList.push_back(tempNeighborList.size());
      for(unsigned int k=0 ; k<tempNeighborList.size() ; ++k)
        neighborhoodList.push_back(tempNeighborList[k]);
    }

    peridigmNeighborhoodIndex += peridigmNumNeighbors;
  }

  // Create the neighborhood list
  neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
  neighborhoodData->SetNumOwned(ownedIds.size());
  memcpy(neighborhoodData->OwnedIDs(), 
 		 &ownedIds[0],
 		 ownedIds.size()*sizeof(int));
  memcpy(neighborhoodData->NeighborhoodPtr(), 
 		 &neighborhoodPtr[0],
 		 neighborhoodPtr.size()*sizeof(int));
  neighborhoodData->SetNeighborhoodListSize(neighborhoodList.size());
  memcpy(neighborhoodData->NeighborhoodList(),
 		 &neighborhoodList[0],
 		 neighborhoodList.size()*sizeof(int));

#else
  neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData(*globalNeighborhoodData));
#endif
  contactNeighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData(*neighborhoodData));
}

void PeridigmNS::ContactManager::initializeContactBlocks()
{
  // Initialize the contact blocks (creates maps, neighborhoods, DataManager)
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++)
    contactBlockIt->initialize(oneDimensionalContactMap,
                               oneDimensionalOverlapContactMap,
                               threeDimensionalContactMap,
                               threeDimensionalOverlapContactMap,
                               bondContactMap,
                               contactBlockIDs,
                               contactNeighborhoodData);

  // Load data from the contact manager's mothership vectors into the contact blocks
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
    contactBlockIt->importData(*contactBlockIDs, blockIdFieldId, PeridigmField::STEP_NONE, Insert);
    contactBlockIt->importData(*contactVolume, volumeFieldId, PeridigmField::STEP_NONE, Insert);
    contactBlockIt->importData(*contactY, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    contactBlockIt->importData(*contactV, velocityFieldId, PeridigmField::STEP_NP1, Insert);
  }
}

void PeridigmNS::ContactManager::importData(Teuchos::RCP<Epetra_Vector> volume,
                                            Teuchos::RCP<Epetra_Vector> coordinates,
                                            Teuchos::RCP<Epetra_Vector> velocity)
{
  // Import data to the contact manager's mothership vectors
  contactY->Import(*coordinates, *threeDimensionalMothershipToContactMothershipImporter, Insert);
  contactV->Import(*velocity, *threeDimensionalMothershipToContactMothershipImporter, Insert);

  // Distribute data to the contact blocks
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
    contactBlockIt->importData(*contactY, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    contactBlockIt->importData(*contactV, velocityFieldId, PeridigmField::STEP_NP1, Insert);
  }
}

void PeridigmNS::ContactManager::exportData(Teuchos::RCP<Epetra_Vector> contactForce)
{
  contactContactForce->PutScalar(0.0);
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
    contactScratch->PutScalar(0.0);
    contactBlockIt->exportData(*contactScratch, contactForceDensityFieldId, PeridigmField::STEP_NP1, Add);
    contactContactForce->Update(1.0, *contactScratch, 1.0);
  }
  // Copy data from the contact mothership vector to the mothership vector
  contactForce->Export(*contactContactForce, *threeDimensionalMothershipToContactMothershipImporter, Insert);
}

void PeridigmNS::ContactManager::rebalance(int step)
{
  if( step%contactRebalanceFrequency != 0)
    return;

  const Epetra_Comm& comm = oneDimensionalMap->Comm();

  // \todo Handle serial case.  We don't need to rebalance, but we still want to update the contact search.
  QUICKGRID::Data rebalancedDecomp = currentConfigurationDecomp();

  Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGridDiscretization::getOwnedMap(comm, rebalancedDecomp, 1)));
  Teuchos::RCP<const Epetra_Import> oneDimensionalMapImporter = Teuchos::rcp(new Epetra_Import(*rebalancedOneDimensionalMap, *oneDimensionalContactMap));

  Teuchos::RCP<Epetra_BlockMap> rebalancedThreeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGridDiscretization::getOwnedMap(comm, rebalancedDecomp, 3)));
  Teuchos::RCP<const Epetra_Import> threeDimensionalMapImporter = Teuchos::rcp(new Epetra_Import(*rebalancedThreeDimensionalMap, *threeDimensionalContactMap));

  Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap = createRebalancedBondMap(rebalancedOneDimensionalMap, oneDimensionalMapImporter);
  Teuchos::RCP<const Epetra_Import> bondMapImporter = Teuchos::rcp(new Epetra_Import(*rebalancedBondMap, *bondContactMap));

  // create a list of neighbors in the rebalanced configuration
  // this list has the global ID for each neighbor of each on-processor point (that is, on processor in the rebalanced configuration)
  Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs = createRebalancedNeighborGlobalIDList(rebalancedBondMap, bondMapImporter);

  // create a list of all the off-processor IDs that will need to be ghosted
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
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 1, indexBase, comm));
  Teuchos::RCP<Epetra_BlockMap> rebalancedThreeDimensionalOverlapMap =
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, 3, indexBase, comm));
  delete[] myGlobalElements;

  // update the current-configuration neighborhood data
  neighborhoodData = createRebalancedNeighborhoodData(rebalancedOneDimensionalMap,
                                                      rebalancedOneDimensionalOverlapMap,
                                                      rebalancedBondMap,
                                                      rebalancedNeighborGlobalIDs);

  // create a new NeighborhoodData object for contact
  contactNeighborhoodData = createRebalancedContactNeighborhoodData(contactNeighborGlobalIDs,
                                                                    rebalancedOneDimensionalMap,
                                                                    rebalancedOneDimensionalOverlapMap);
  
  // rebalance the mothership (global) contact vectors
  Teuchos::RCP<Epetra_MultiVector> rebalancedOneDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*rebalancedOneDimensionalMap, oneDimensionalContactMothership->NumVectors()));
  rebalancedOneDimensionalMothership->Import(*oneDimensionalContactMothership, *oneDimensionalMapImporter, Insert);
  oneDimensionalContactMothership = rebalancedOneDimensionalMothership;
  contactBlockIDs = Teuchos::rcp((*oneDimensionalContactMothership)(0), false);         // block ID
  contactVolume = Teuchos::rcp((*oneDimensionalContactMothership)(1), false);           // cell volume

  Teuchos::RCP<Epetra_MultiVector> rebalancedThreeDimensionalMothership = Teuchos::rcp(new Epetra_MultiVector(*rebalancedThreeDimensionalMap, threeDimensionalContactMothership->NumVectors()));
  rebalancedThreeDimensionalMothership->Import(*threeDimensionalContactMothership, *threeDimensionalMapImporter, Insert);
  threeDimensionalContactMothership = rebalancedThreeDimensionalMothership;
  contactY = Teuchos::rcp((*threeDimensionalContactMothership)(0), false);             // current positions
  contactV = Teuchos::rcp((*threeDimensionalContactMothership)(1), false);             // velocities
  contactContactForce = Teuchos::rcp((*threeDimensionalContactMothership)(2), false);  // contact force
  contactScratch = Teuchos::rcp((*threeDimensionalContactMothership)(3), false);       // scratch

  // rebalance the contact blocks
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++)
    contactBlockIt->rebalance(rebalancedOneDimensionalMap,
                              rebalancedOneDimensionalOverlapMap,
                              rebalancedThreeDimensionalMap,
                              rebalancedThreeDimensionalOverlapMap,
                              rebalancedBondMap,
                              contactBlockIDs,
                              contactNeighborhoodData);

  // Reload data from the contact manager's mothership vectors into the contact blocks
  // \todo Cut back on loading data, it's probably correct already depending on how rebalance is handled above.
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){
    contactBlockIt->importData(*contactBlockIDs, blockIdFieldId, PeridigmField::STEP_NONE, Insert);
    contactBlockIt->importData(*contactVolume, volumeFieldId, PeridigmField::STEP_NONE, Insert);
    contactBlockIt->importData(*contactY, coordinatesFieldId, PeridigmField::STEP_NP1, Insert);
    contactBlockIt->importData(*contactV, velocityFieldId, PeridigmField::STEP_NP1, Insert);
  }

  // set all the pointers to the new maps
  oneDimensionalContactMap = rebalancedOneDimensionalMap;
  oneDimensionalOverlapContactMap = rebalancedOneDimensionalOverlapMap;
  threeDimensionalOverlapContactMap = rebalancedThreeDimensionalOverlapMap;
  threeDimensionalContactMap = rebalancedThreeDimensionalMap;
  bondContactMap = rebalancedBondMap;

  // Reset the importers for passing data between the mothership and contact mothership vectors
  oneDimensionalMothershipToContactMothershipImporter = Teuchos::rcp(new Epetra_Import(*oneDimensionalContactMap, *oneDimensionalMap));
  threeDimensionalMothershipToContactMothershipImporter = Teuchos::rcp(new Epetra_Import(*threeDimensionalContactMap, *threeDimensionalMap));
}

QUICKGRID::Data PeridigmNS::ContactManager::currentConfigurationDecomp() {

  // Create a decomp object and fill necessary data for rebalance
  int myNumElements = oneDimensionalContactMap->NumMyElements();
  int dimension = 3;
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(myNumElements, dimension);

  decomp.globalNumPoints = oneDimensionalContactMap->NumGlobalElements();

  // fill myGlobalIDs
  UTILITIES::Array<int> myGlobalIDs(myNumElements);
  int* myGlobalIDsPtr = myGlobalIDs.get();
  int* gIDs = oneDimensionalContactMap->MyGlobalElements();
  memcpy(myGlobalIDsPtr, gIDs, myNumElements*sizeof(int));
  decomp.myGlobalIDs = myGlobalIDs.get_shared_ptr();

  // fill myX
  // use current positions for x
  UTILITIES::Array<double> myX(myNumElements*dimension);
  double* myXPtr = myX.get();
  double* yPtr;
  contactY->ExtractView(&yPtr);
  memcpy(myXPtr, yPtr, myNumElements*dimension*sizeof(double));
  decomp.myX = myX.get_shared_ptr();

  // fill cellVolume
  UTILITIES::Array<double> cellVolume(myNumElements);
  double* cellVolumePtr = cellVolume.get();
  double* volumePtr;
  contactVolume->ExtractView(&volumePtr);
  memcpy(cellVolumePtr, volumePtr, myNumElements*sizeof(double));
  decomp.cellVolume = cellVolume.get_shared_ptr();

  // call the rebalance function on the current-configuration decomp
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

  return decomp;
}

Teuchos::RCP<Epetra_BlockMap> PeridigmNS::ContactManager::createRebalancedBondMap(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                  Teuchos::RCP<const Epetra_Import> oneDimensionalMapToRebalancedOneDimensionalMapImporter) {

  const Epetra_Comm& comm = oneDimensionalContactMap->Comm();

  // communicate the number of bonds for each point so that space for bond data can be allocated
  Teuchos::RCP<Epetra_Vector> numberOfBonds = Teuchos::rcp(new Epetra_Vector(*oneDimensionalContactMap));
  for(int i=0 ; i<oneDimensionalContactMap->NumMyElements() ; ++i){
    int globalID = oneDimensionalContactMap->GID(i);
    int bondMapLocalID = bondContactMap->LID(globalID);
    if(bondMapLocalID != -1)
      (*numberOfBonds)[i] = (double)( bondContactMap->ElementSize(i) );
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
    Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, comm));
  delete[] myGlobalElements;
  delete[] elementSizeList;

  return rebalancedBondMap;
}

template<class T>
struct NonDeleter{
	void operator()(T* d) {}
};

void PeridigmNS::ContactManager::contactSearch(Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap, 
                                               Teuchos::RCP<const Epetra_BlockMap> rebalancedBondMap,
                                               Teuchos::RCP<const Epetra_Vector> rebalancedNeighborGlobalIDs,
                                               QUICKGRID::Data& rebalancedDecomp,
                                               Teuchos::RCP< map<int, vector<int> > > contactNeighborGlobalIDs,
                                               Teuchos::RCP< set<int> > offProcessorContactIDs)
{
  const Epetra_Comm& comm = oneDimensionalMap->Comm();

  std::tr1::shared_ptr<const Epetra_Comm> comm_shared_ptr(&comm,NonDeleter<const Epetra_Comm>());
  QUICKGRID::Data d = rebalancedDecomp;

  // TEMPORARY PLACEHOLDER FOR PER-NODE SEARCH RADII
  Teuchos::RCP<Epetra_Vector> contactSearchRadii = Teuchos::rcp(new Epetra_Vector(*rebalancedOneDimensionalMap));
  contactSearchRadii->PutScalar(contactSearchRadius);

  PDNEIGH::NeighborhoodList neighList(comm_shared_ptr,d.zoltanPtr.get(),d.numPoints,d.myGlobalIDs,d.myX,contactSearchRadii);

  int* searchNeighborhood = neighList.get_neighborhood().get();

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
      list<int>::iterator it = find(bondedNeighbors.begin(), bondedNeighbors.end(), globalNeighborID);  // \todo Don't consider broken bonds here
      if(it == bondedNeighbors.end()){
        contactNeighborGlobalIDList.push_back(globalNeighborID);
        if(rebalancedOneDimensionalMap->LID(globalNeighborID) == -1)
          offProcessorContactIDs->insert(globalNeighborID);
      }
    }
  }
}

Teuchos::RCP<Epetra_Vector> PeridigmNS::ContactManager::createRebalancedNeighborGlobalIDList(Teuchos::RCP<Epetra_BlockMap> rebalancedBondMap,
                                                                                             Teuchos::RCP<const Epetra_Import> bondMapToRebalancedBondMapImporter) {
  // construct a globalID neighbor list in the static global decomposition
  Teuchos::RCP<Epetra_Vector> neighborGlobalIDs = Teuchos::rcp(new Epetra_Vector(*bondContactMap));
  int* neighborhoodList = neighborhoodData->NeighborhoodList();
  int neighborhoodListIndex = 0;
  int neighborGlobalIDIndex = 0;
  for(int i=0 ; i<neighborhoodData->NumOwnedPoints() ; ++i){
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(int j=0 ; j<numNeighbors ; ++j){
      int neighborLocalID = neighborhoodList[neighborhoodListIndex++];
      TEUCHOS_TEST_FOR_EXCEPTION(neighborGlobalIDIndex >= neighborGlobalIDs->MyLength(), Teuchos::RangeError, "ContactManager::createRebalancedNeighborGlobalIDList(), Invalid index into neighborGlobalIDs\n");
      (*neighborGlobalIDs)[neighborGlobalIDIndex++] = oneDimensionalOverlapContactMap->GID(neighborLocalID);
    }
  }

  // redistribute the globalID neighbor list to the rebalanced configuration
  Teuchos::RCP<Epetra_Vector> rebalancedNeighborGlobalIDs = Teuchos::rcp(new Epetra_Vector(*rebalancedBondMap));
  rebalancedNeighborGlobalIDs->Import(*neighborGlobalIDs, *bondMapToRebalancedBondMapImporter, Insert);

  return rebalancedNeighborGlobalIDs;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> PeridigmNS::ContactManager::createRebalancedNeighborhoodData(Teuchos::RCP<Epetra_BlockMap> rebalancedOneDimensionalMap,
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


Teuchos::RCP<PeridigmNS::NeighborhoodData> PeridigmNS::ContactManager::createRebalancedContactNeighborhoodData(Teuchos::RCP<map<int, vector<int> > > contactNeighborGlobalIDs,
                                                                                                               Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalMap,
                                                                                                               Teuchos::RCP<const Epetra_BlockMap> rebalancedOneDimensionalOverlapMap)
{
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

void PeridigmNS::ContactManager::evaluateContactForce(double dt)
{
  for(contactBlockIt = contactBlocks->begin() ; contactBlockIt != contactBlocks->end() ; contactBlockIt++){

    Teuchos::RCP<PeridigmNS::NeighborhoodData> nData = contactBlockIt->getNeighborhoodData();
    const int numOwnedPoints = nData->NumOwnedPoints();
    const int* ownedIDs = nData->OwnedIDs();
    const int* neighborhoodList = nData->NeighborhoodList();
    Teuchos::RCP<PeridigmNS::DataManager> dataManager = contactBlockIt->getDataManager();
    Teuchos::RCP<const PeridigmNS::ContactModel> contactModel = contactBlockIt->getContactModel();

    if(!contactModel.is_null())
      contactModel->computeForce(dt, 
                                 numOwnedPoints,
                                 ownedIDs,
                                 neighborhoodList,
                                 *dataManager);
  }
}
