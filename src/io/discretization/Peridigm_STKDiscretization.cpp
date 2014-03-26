/*! \file Peridigm_STKDiscretization.cpp */

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

#include "Peridigm_STKDiscretization.hpp"
#include "Peridigm_ProximitySearch.hpp"
#include "Peridigm_HorizonManager.hpp"
#include "Peridigm_GeometryUtils.hpp"
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>
#include <Ionit_Initializer.h>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <sstream>
#include <boost/math/constants/constants.hpp>

using namespace std;

PeridigmNS::STKDiscretization::STKDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) :
  minElementRadius(1.0e50),
  maxElementRadius(0.0),
  storeExodusMesh(false),
  computeIntersections(false), /* todo, set this only if needed based on params */
  maxElementDimension(0.0),
  numBonds(0),
  maxNumBondsPerElem(0),
  myPID(epetra_comm->MyPID()),
  numPID(epetra_comm->NumProc()),
  bondFilterCommand("None"),
  comm(epetra_comm)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Exodus", "Invalid Type in STKDiscretization");

  if(params->isParameter("Omit Bonds Between Blocks"))
    bondFilterCommand = params->get<string>("Omit Bonds Between Blocks");
  string meshFileName = params->get<string>("Input Mesh File");

  // Store exodus mesh for intersection calculations, or if it was specifically requested (e.g., unit tests)
  if(params->isParameter("Store Exodus Mesh"))
    storeExodusMesh = params->get<bool>("Store Exodus Mesh");
  if(computeIntersections)
    storeExodusMesh = true;

  // Set up bond filters
  createBondFilters(params);

  // Load data from mesh file
  loadData(meshFileName);
  
  if(computeIntersections)
    maxElementDimension = computeMaxElementDimension();

  // Assign the correct horizon to each node
  PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();
  horizonForEachPoint = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  for(map<string, vector<int> >::const_iterator it = elementBlocks->begin() ; it != elementBlocks->end() ; it++){
    const string& blockName = it->first;
    const vector<int>& globalIds = it->second;

    bool hasConstantHorizon = horizonManager.blockHasConstantHorizon(blockName);
    double constantHorizonValue(0.0);
    if(hasConstantHorizon)
      constantHorizonValue = horizonManager.getBlockConstantHorizonValue(blockName);

    for(unsigned int i=0 ; i<globalIds.size() ; ++i){
      int localId = oneDimensionalMap->LID(globalIds[i]);
      if(hasConstantHorizon){
        (*horizonForEachPoint)[localId] = constantHorizonValue;
      }
      else{
        double x = (*initialX)[localId*3];
        double y = (*initialX)[localId*3 + 1];
        double z = (*initialX)[localId*3 + 2];
        double horizon = horizonManager.evaluateHorizon(blockName, x, y, z);
        (*horizonForEachPoint)[localId] = horizon;
      }
    }
  }

  // Execute the neighbor search
  int neighborListSize;
  int* neighborList;
  if(!computeIntersections){
    ProximitySearch::GlobalProximitySearch(initialX, horizonForEachPoint, oneDimensionalOverlapMap, neighborListSize, neighborList, bondFilters);
  }
  else{
    ProximitySearch::GlobalProximitySearch(initialX, horizonForEachPoint, oneDimensionalOverlapMap, neighborListSize, neighborList, bondFilters, maxElementDimension);
    removeNonintersectingNeighborsFromNeighborList(initialX, horizonForEachPoint, oneDimensionalOverlapMap, neighborListSize, neighborList);
  }

  createNeighborhoodData(neighborListSize, neighborList);

  // Create the three-dimensional overlap map based on the one-dimensional overlap map
  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(-1, 
                                                                oneDimensionalOverlapMap->NumMyElements(),
                                                                oneDimensionalOverlapMap->MyGlobalElements(),
                                                                3,
                                                                0,
                                                                oneDimensionalOverlapMap->Comm()));

  // \todo Move this functionality to base class, it's currently duplicated in PdQuickGridDiscretization.
  // Create the bondMap, a local map used for constitutive data stored on bonds.
  // Due to Epetra_BlockMap restrictions, there can not be any entries with length zero.
  // This means that points with no neighbors can not appear in the bondMap.
  int numMyElementsUpperBound = oneDimensionalMap->NumMyElements();
  int numGlobalElements = -1; 
  int numMyElements = 0;
  int maxNumBonds = 0;
  int* oneDimensionalMapGlobalElements = oneDimensionalMap->MyGlobalElements();
  int* myGlobalElements = new int[numMyElementsUpperBound];
  int* elementSizeList = new int[numMyElementsUpperBound];
  int* const neighborhood = neighborhoodData->NeighborhoodList();
  int neighborhoodIndex = 0;
  int numPointsWithZeroNeighbors = 0;
  for(int i=0 ; i<neighborhoodData->NumOwnedPoints() ; ++i){
    int numNeighbors = neighborhood[neighborhoodIndex];
    if(numNeighbors > 0){
      numMyElements++;
      myGlobalElements[i-numPointsWithZeroNeighbors] = oneDimensionalMapGlobalElements[i];
      elementSizeList[i-numPointsWithZeroNeighbors] = numNeighbors;
    }
    else{
      numPointsWithZeroNeighbors++;
    }
    numBonds += numNeighbors;
    if(numNeighbors>maxNumBonds) maxNumBonds = numNeighbors;
    neighborhoodIndex += 1 + numNeighbors;
  }
  maxNumBondsPerElem = maxNumBonds;
  int indexBase = 0;
  bondMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, *comm));
  delete[] myGlobalElements;
  delete[] elementSizeList;

  // find the minimum element radius
  for(int i=0 ; i<cellVolume->MyLength() ; ++i){
    double radius = pow(0.238732414637843*(*cellVolume)[i], 0.33333333333333333);
    if(radius < minElementRadius)
      minElementRadius = radius;
    if(radius > maxElementRadius)
      maxElementRadius = radius;
  }
  vector<double> localMin(1);
  vector<double> globalMin(1);
  localMin[0] = minElementRadius;
  epetra_comm->MinAll(&localMin[0], &globalMin[0], 1);
  minElementRadius = globalMin[0];
  localMin[0] = maxElementRadius;
  epetra_comm->MaxAll(&localMin[0], &globalMin[0], 1);
  maxElementRadius = globalMin[0];
}

PeridigmNS::STKDiscretization::~STKDiscretization() {}

void PeridigmNS::STKDiscretization::loadData(const string& meshFileName)
{

  string meshType = "exodusii";
  string workingDirectory = "";
  Teuchos::RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(comm, true);
  metaData = Teuchos::rcp(new stk::mesh::fem::FEMMetaData);

  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  meshData = Teuchos::rcp(new stk::io::MeshData);
  #else
  meshData = Teuchos::rcp(new stk::io::util::MeshData);
  #endif
  Ioss::Init::Initializer io;

  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  stk::io::create_input_mesh(meshType,
                             meshFileName,
                             mpiComm->Comm(),
                             *metaData,
                             *meshData);
  #else
  stk::io::util::create_input_mesh(meshType,
                                   meshFileName,
                                   workingDirectory,
                                   mpiComm->Comm(),
                                   *metaData,
                                   *meshData);
  #endif
  
  int numberOfDimensions = metaData->spatial_dimension();
  TEUCHOS_TEST_FOR_EXCEPTION(numberOfDimensions != 3, std::invalid_argument, "Peridigm operates only on three-dimensional meshes.");

  // This assigns a null Ioss::GroupingEntity attribute to the universal part
  stk::io::put_io_part_attribute(metaData->universal_part());

  // Loop over the parts and store the element parts and node parts.
  const stk::mesh::PartVector& stkParts = metaData->get_parts();
  stk::mesh::PartVector stkElementBlocks;
//   stk::mesh::PartVector stkSideSets;
  stk::mesh::PartVector stkNodeSets;
  for(stk::mesh::PartVector::const_iterator it = stkParts.begin(); it != stkParts.end(); ++it){
    stk::mesh::Part* const part = *it;
    if(part->name()[0] == '{')
      continue;
    if(part->primary_entity_rank() == metaData->element_rank())
      stkElementBlocks.push_back(part);
//     else if(part->primary_entity_rank() == metaData->side_rank())
//       stkSideSets.push_back(part);
    else if(part->primary_entity_rank() == metaData->node_rank())
      stkNodeSets.push_back(part);
    else
      if(myPID == 0)
        cout << "Warning, unknown part type for part " << part->name() << endl;
  }
  if(myPID == 0){
    stringstream ss;
    ss << "Converting input file " << meshFileName << " to sphere mesh:" << endl;
    ss << "  Element blocks:";
    for(stk::mesh::PartVector::const_iterator it = stkElementBlocks.begin(); it != stkElementBlocks.end(); ++it)
      ss << " " << (*it)->name();
    ss << endl;
//     ss << "  Side sets:";
//     for(stk::mesh::PartVector::const_iterator it = stkSideSets.begin(); it != stkSideSets.end(); ++it)
//       ss << " " << (*it)->name();
//     ss << endl;
    ss << "  Node sets:";
    for(stk::mesh::PartVector::const_iterator it = stkNodeSets.begin(); it != stkNodeSets.end(); ++it)
      ss << " " << (*it)->name();
    ss << endl;
    cout << ss.str() << endl;
  }

  if (!metaData->is_FEM_initialized())
    metaData->FEM_initialize(numberOfDimensions);

  stk::mesh::BulkData bulkData(stk::mesh::fem::FEMMetaData::get_meta_data(*metaData), mpiComm->Comm());

  metaData->commit();

  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  stk::io::populate_bulk_data(bulkData, *meshData);
  #else
  stk::io::util::populate_bulk_data(bulkData, *meshData, "exodusii");
  #endif
  bulkData.modification_end();

  stk::mesh::Field<double, stk::mesh::Cartesian>* coordinatesField = 
    metaData->get_field< stk::mesh::Field<double, stk::mesh::Cartesian> >("coordinates");

  // The volume field is present only for sphere meshes
  // volumeField will be a null pointer for tet or hex meshes
  stk::mesh::Field<double, stk::mesh::Cartesian>* volumeField = 
    metaData->get_field< stk::mesh::Field<double, stk::mesh::Cartesian> >("volume");

  // Create a selector to select everything in the universal part that is either locally owned or globally shared
  stk::mesh::Selector selector = 
    stk::mesh::Selector( metaData->universal_part() ) & ( stk::mesh::Selector( metaData->locally_owned_part() ) | stk::mesh::Selector( metaData->globally_shared_part() ) );

  // Select element mesh entities that match the selector
  std::vector<stk::mesh::Entity*> elements;
  stk::mesh::get_selected_entities(selector, bulkData.buckets(metaData->element_rank()), elements);

  // Select node mesh entities that match the selector
  std::vector<stk::mesh::Entity*> nodes;
  stk::mesh::get_selected_entities(selector, bulkData.buckets(metaData->node_rank()), nodes);

  // Determine the total number of elements in the model
  // \todo There must be a cleaner way to determine the number of elements in a model.
  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  int localElemCount = elements.size();
  int globalElemCount(0);
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, int(1), &localElemCount, &globalElemCount);

  // Store the positions of the nodes in the initial mesh
  // This is used for the calculation of partial neighbor volumes
  if(storeExodusMesh){
    vector<int> myGlobalNodeIds(nodes.size());
    for(unsigned int i=0 ; i<nodes.size() ; ++i)
      myGlobalNodeIds[i] = nodes[i]->identifier() - 1;
    Epetra_BlockMap exodusMeshNodePositionsMap(-1, nodes.size(), &myGlobalNodeIds[0], 3, 0, *comm);
    exodusMeshNodePositions = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(exodusMeshNodePositionsMap));
    for(unsigned int iNode=0 ; iNode<nodes.size() ; ++iNode){
      int globalID = nodes[iNode]->identifier() - 1;
      double* coordinates = stk::mesh::field_data(*coordinatesField, *nodes[iNode]);
      int localID = exodusMeshNodePositionsMap.LID(globalID);
      (*exodusMeshNodePositions)[3*localID]   = coordinates[0];
      (*exodusMeshNodePositions)[3*localID+1] = coordinates[1];
      (*exodusMeshNodePositions)[3*localID+2] = coordinates[2];
    }
  }

  // Create a list of owned global element ids
  vector<int> globalElementIds(elements.size());
  for(unsigned int iElem=0 ; iElem<elements.size() ; ++iElem){
    stk::mesh::PairIterRelation nodeRelations = elements[iElem]->node_relations();
    int elementId = elements[iElem]->identifier() - 1;
    globalElementIds[iElem] = elementId;
  }

  // Create the owned maps
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(-1, elements.size(), &globalElementIds[0], 1, 0, *comm));
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(-1, elements.size(), &globalElementIds[0], 3, 0, *comm));

  // Create Epetra_Vectors for the initial positions, volumes, and block_ids
  initialX = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  cellVolume = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  blockID = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));

  double* initialXPtr;
  initialX->ExtractView(&initialXPtr);
  double* cellVolumePtr;
  cellVolume->ExtractView(&cellVolumePtr);

  // warning flags
  bool tenNodedTetWarningGiven = false;
  bool twentyNodedHexWarningGiven = false;

  // loop over the elements and load the data into the discretization's data structures
  for(unsigned int iElem=0 ; iElem<elements.size() ; ++iElem){
    stk::mesh::PairIterRelation nodeRelations = elements[iElem]->node_relations();
    int elementID = elements[iElem]->identifier() - 1;
    if(storeExodusMesh)
      exodusMeshElementConnectivity[elementID] = vector<int>();
    initialXPtr[iElem*3] = 0.0;
    initialXPtr[iElem*3+1] = 0.0;
    initialXPtr[iElem*3+2] = 0.0;
    std::vector<double*> nodeCoordinates(nodeRelations.size());
    int iNode = 0;
    for(stk::mesh::PairIterRelation::iterator it=nodeRelations.begin() ; it!=nodeRelations.end() ; ++it){
      stk::mesh::Entity* node = it->entity();
      if(storeExodusMesh)
        exodusMeshElementConnectivity[elementID].push_back(node->identifier() - 1);
      double* coordinates = stk::mesh::field_data(*coordinatesField, *node);
      initialXPtr[iElem*3] += coordinates[0];
      initialXPtr[iElem*3+1] += coordinates[1];
      initialXPtr[iElem*3+2] += coordinates[2];
      nodeCoordinates[iNode++] = coordinates;
    }
    initialXPtr[iElem*3] /= nodeRelations.size();
    initialXPtr[iElem*3+1] /= nodeRelations.size();
    initialXPtr[iElem*3+2] /= nodeRelations.size();
    if(nodeRelations.size() == 1){
      double* exodusVolume = stk::mesh::field_data(*volumeField, *elements[iElem]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(exodusVolume == NULL, "**** Volume attribute not found for sphere element.\n");
      cellVolumePtr[iElem] = exodusVolume[0];
    }
    else if(nodeRelations.size() == 4){
      // 4-noded tet
      cellVolumePtr[iElem] = tetVolume(nodeCoordinates);
    }
    else if(nodeRelations.size() == 8){
      // 8-noded hex
      cellVolumePtr[iElem] = hexVolume(nodeCoordinates);
    }
    else if(nodeRelations.size() == 10){
      // 10-noded tet, treat as 4-noded tet
      if(!tenNodedTetWarningGiven){
        cout << "**** Warning on processor " << myPID << ", side nodes being discarded for 10-node tetrahedron element, will be treated as 4-node tetrahedron element." << endl;
        tenNodedTetWarningGiven = true;
      }
      cellVolumePtr[iElem] = hexVolume(nodeCoordinates);
    }
    else if(nodeRelations.size() == 20){
      // 20-noded hex, treat as 8-noded tet
      if(!twentyNodedHexWarningGiven){
        cout << "**** Warning on processor " << myPID << ", side nodes being discarded for 20-node hexahedron element, will be treated as 8-node hexahedron element." << endl;
        twentyNodedHexWarningGiven = true;
      }
      cellVolumePtr[iElem] = hexVolume(nodeCoordinates);
    }
    else{
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Invalid element topology.\n");
    }
  }

  // loop over the element blocks
  for(unsigned int iBlock=0 ; iBlock<stkElementBlocks.size() ; iBlock++){

    const std::string blockName = stkElementBlocks[iBlock]->name();
    TEUCHOS_TEST_FOR_EXCEPT_MSG(elementBlocks->find(blockName) != elementBlocks->end(), "**** Duplicate block found: " + blockName + "\n");
    (*elementBlocks)[blockName] = std::vector<int>();
    std::vector<int>& elementBlock = (*elementBlocks)[blockName];

    // Create a selector for all locally-owned elements in the block
    stk::mesh::Selector selector = 
      stk::mesh::Selector( *stkElementBlocks[iBlock] ) & stk::mesh::Selector( metaData->locally_owned_part() );

    // Select the mesh entities that match the selector
    std::vector<stk::mesh::Entity*> elementsInElementBlock;
    stk::mesh::get_selected_entities(selector, bulkData.buckets(metaData->element_rank()), elementsInElementBlock);

    // Loop over the elements in this block
    std::vector<int> elementGlobalIds(elementsInElementBlock.size());

    for(unsigned int iElement=0 ; iElement<elementsInElementBlock.size() ; iElement++)
      elementGlobalIds[iElement] = elementsInElementBlock[iElement]->identifier() - 1;

    // Load the element block into the elementBlock container
    for(unsigned int i=0 ; i<elementGlobalIds.size() ; ++i)
      elementBlock.push_back( elementGlobalIds[i] );
  }

  // loop over the node sets
  for(unsigned int iNodeSet=0 ; iNodeSet<stkNodeSets.size() ; iNodeSet++){

    const std::string nodeSetName = stkNodeSets[iNodeSet]->name();
    TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSetName) != nodeSets->end(), "**** Duplicate node set found: " + nodeSetName + "\n");
    (*nodeSets)[nodeSetName] = std::vector<int>();
    std::vector<int>& nodeSet = (*nodeSets)[nodeSetName];

    // Create a selector for all nodes in the node set
    stk::mesh::Selector selector = stk::mesh::Selector( *stkNodeSets[iNodeSet] );

    // Select the mesh entities that match the selector
    std::vector<stk::mesh::Entity*> nodesInNodeSet;
    stk::mesh::get_selected_entities(selector, bulkData.buckets(metaData->node_rank()), nodesInNodeSet);
    
    // Loop over the nodes in this node set
    std::set<int> elementGlobalIds;
    for(unsigned int iNode=0 ; iNode<nodesInNodeSet.size() ; iNode++){

      // Get all the elements relations for this node and record the element ID in the node set
      // This works because the original element ID is the same as the node ID in the sphere mesh
      stk::mesh::PairIterRelation elementRelations = nodesInNodeSet[iNode]->relations(metaData->element_rank()); 
      for(stk::mesh::PairIterRelation::iterator it=elementRelations.begin() ; it!=elementRelations.end() ; ++it){
        stk::mesh::Entity* element = it->entity();
        int globalId = element->identifier() - 1;
        // Determine if the element is on-processor
        bool onProcessor = false;
        for(std::map< std::string, std::vector<int> >::iterator b_it = elementBlocks->begin() ; b_it != elementBlocks->end() && !onProcessor ; b_it++){
          if( find(b_it->second.begin(), b_it->second.end(), globalId) != b_it->second.end() )
            onProcessor = true;
        }
        // If the element is on processor, add the corresponding node to the node list
        if(onProcessor)
          elementGlobalIds.insert(globalId);
      }
    }

    // Load the node set into the nodeSets container
    for(std::set<int>::iterator it=elementGlobalIds.begin() ; it!=elementGlobalIds.end() ; it++)
      nodeSet.push_back( *it );
  }

  // Record the element block ids
  double* blockIDPtr;
  blockID->ExtractView(&blockIDPtr);
  std::map< std::string, std::vector<int> >::const_iterator it;
  for(it = elementBlocks->begin() ; it != elementBlocks->end() ; it++){
    const std::string& blockName = it->first;
    int bID = blockNameToBlockId(blockName);
    const std::vector<int>& elementIDs = it->second;
    for(unsigned int i=0 ; i<elementIDs.size() ; ++i){
      int globalID = elementIDs[i];
      int localID = oneDimensionalMap->LID(globalID);
      blockIDPtr[localID] = bID;
    }
  }

  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  // free the meshData
  meshData = Teuchos::RCP<stk::io::MeshData>();
  #else
  // free the meshData
  meshData = Teuchos::RCP<stk::io::util::MeshData>();
  #endif

  return;
}

void
PeridigmNS::STKDiscretization::createNeighborhoodData(int neighborListSize, int* neighborList)
{

   int numOwnedIds = oneDimensionalMap->NumMyElements();
   
   int* ownedGlobalIds = oneDimensionalMap->MyGlobalElements();

   vector<int> ownedLocalIds(numOwnedIds);
   vector<int> neighborhoodPtr(numOwnedIds);

   int numNeighbors(0), neighborListIndex(0);
   for(int i=0 ; i<numOwnedIds ; ++i){
     ownedLocalIds[i] = oneDimensionalMap->LID(ownedGlobalIds[i]); // \todo This seems unnecessary, is it just i?  What if MyGlobalElements is not sorted, then is ownedLocalIds also not sorted?
     neighborhoodPtr[i] = neighborListIndex;     
     numNeighbors = neighborList[neighborListIndex++];
     neighborListIndex += numNeighbors;
   }

   neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
   neighborhoodData->SetNumOwned(numOwnedIds);
   memcpy(neighborhoodData->OwnedIDs(), &ownedLocalIds[0], numOwnedIds*sizeof(int));
   memcpy(neighborhoodData->NeighborhoodPtr(), &neighborhoodPtr[0], numOwnedIds*sizeof(int));
   neighborhoodData->SetNeighborhoodListSize(neighborListSize);
   memcpy(neighborhoodData->NeighborhoodList(), neighborList, neighborListSize*sizeof(int));
   neighborhoodData = filterBonds(neighborhoodData);
}

Teuchos::RCP<PeridigmNS::NeighborhoodData>
PeridigmNS::STKDiscretization::filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData)
{
  // Set up a block bonding matrix, which defines whether or not bonds should be formed across blocks
  int numBlocks = getNumBlocks();
  std::vector< std::vector<bool> > blockBondingMatrix(numBlocks);
  for(int i=0 ; i<numBlocks ; ++i){
    blockBondingMatrix[i].resize(numBlocks, true);
  }

  if(bondFilterCommand == "None"){
    // All blocks are bonded, the blockBondingMatrix is unchanged
    return unfilteredNeighborhoodData;
  }
  else if(bondFilterCommand == "All"){
    // No blocks are bonded, the blockBondingMatrix is the identity matrix
    for(int i=0 ; i<numBlocks ; ++i){
      for(int j=0 ; j<numBlocks ; ++j){
        if(i != j)
          blockBondingMatrix[i][j] = false;
      }
    }
  }
  else{
    string msg = "**** Error, unrecognized value for \"Omit Bonds Between Blocks\":  ";
    msg += bondFilterCommand + "\n";
    msg += "**** Valid options are:  All, None\n";
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
  }

  // Create an overlap vector containing the block IDs of each cell
  Teuchos::RCP<const Epetra_BlockMap> ownedMap = getGlobalOwnedMap(1);
  Teuchos::RCP<const Epetra_BlockMap> overlapMap = getGlobalOverlapMap(1);
  Epetra_Vector blockIDs(*overlapMap);
  Epetra_Import importer(*overlapMap, *ownedMap);
  Teuchos::RCP<Epetra_Vector> ownedBlockIDs = getBlockID();
  blockIDs.Import(*ownedBlockIDs, importer, Insert);

  // Apply the block bonding matrix and create a new NeighborhoodData
  Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
  neighborhoodData->SetNumOwned(unfilteredNeighborhoodData->NumOwnedPoints());
  memcpy(neighborhoodData->OwnedIDs(), unfilteredNeighborhoodData->OwnedIDs(), neighborhoodData->NumOwnedPoints()*sizeof(int));
  vector<int> neighborhoodListVec;
  neighborhoodListVec.reserve(unfilteredNeighborhoodData->NeighborhoodListSize());
  int* const neighborhoodPtr = neighborhoodData->NeighborhoodPtr();

  int numOwnedPoints = neighborhoodData->NumOwnedPoints();
  int* const unfilteredNeighborhoodList = unfilteredNeighborhoodData->NeighborhoodList();
  int unfilteredNeighborhoodListIndex(0);
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    int blockID = static_cast<int>(blockIDs[iID]);
	int numUnfilteredNeighbors = unfilteredNeighborhoodList[unfilteredNeighborhoodListIndex++];
    unsigned int numNeighborsIndex = neighborhoodListVec.size();
    neighborhoodListVec.push_back(-1); // placeholder for number of neighbors
    int numNeighbors = 0;
	for(int iNID=0 ; iNID<numUnfilteredNeighbors ; ++iNID){
      int unfilteredNeighborID = unfilteredNeighborhoodList[unfilteredNeighborhoodListIndex++];
      int unfilteredNeighborBlockID = static_cast<int>(blockIDs[unfilteredNeighborID]);
      if(blockBondingMatrix[blockID-1][unfilteredNeighborBlockID-1] == true){
        neighborhoodListVec.push_back(unfilteredNeighborID);
        numNeighbors += 1;
      }
    }
    neighborhoodListVec[numNeighborsIndex] = numNeighbors;
    neighborhoodPtr[iID] = numNeighborsIndex;
  }

  neighborhoodData->SetNeighborhoodListSize(neighborhoodListVec.size());
  memcpy(neighborhoodData->NeighborhoodList(), &neighborhoodListVec[0], neighborhoodListVec.size()*sizeof(int));

  return neighborhoodData;
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::STKDiscretization::getGlobalOwnedMap(int d) const
{
  switch (d) {
    case 1:
      return oneDimensionalMap;
      break;
    case 3:
      return threeDimensionalMap;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         std::endl << "STKDiscretization::getGlobalOwnedMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::STKDiscretization::getGlobalOverlapMap(int d) const
{
  switch (d) {
    case 1:
      return oneDimensionalOverlapMap;
      break;
    case 3:
      return threeDimensionalOverlapMap;
      break;
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         std::endl << "STKDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::STKDiscretization::getGlobalBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::STKDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::STKDiscretization::getHorizon() const
{
  return horizonForEachPoint;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::STKDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::STKDiscretization::getBlockID() const
{
  return blockID;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::STKDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::STKDiscretization::getNumBonds() const
{
  return numBonds;
}

unsigned int
PeridigmNS::STKDiscretization::getMaxNumBondsPerElem() const
{
  return maxNumBondsPerElem;
}

void PeridigmNS::STKDiscretization::getExodusMeshNodePositions(int globalNodeID,
                                                               vector<double>& nodePositions)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!storeExodusMesh, "**** Error:  getExodusMeshNodePositions() called, but exodus information not stored.\n");
  vector<int>& elementConnectivity = exodusMeshElementConnectivity[globalNodeID];
  unsigned int numNodes = elementConnectivity.size();
  if(nodePositions.size() != 3*numNodes){
    nodePositions.resize(3*numNodes);
  }
  Epetra_Vector& exodusNodePositions = *exodusMeshNodePositions;
  for(unsigned int i=0 ; i<numNodes ; ++i){
    int globalID = elementConnectivity[i];
    int localID = exodusMeshNodePositions->Map().LID(globalID);
    nodePositions[3*i]   = exodusNodePositions[3*localID];
    nodePositions[3*i+1] = exodusNodePositions[3*localID+1];
    nodePositions[3*i+2] = exodusNodePositions[3*localID+2];
  }
  return;
}

double PeridigmNS::STKDiscretization::computeMaxElementDimension()
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(!storeExodusMesh, "**** Error:  computeMaxElementDimension() called, but exodus information not stored.\n");
  double length, volume, localMaxElementDimension(0.0);
  double x, y, z;
  int globalId, numNodes, numNodesInElement;
  vector<double> nodeCoordinates;
  const double pi = boost::math::constants::pi<double>();

  for(int iElem=0 ; iElem<oneDimensionalMap->NumMyElements() ; ++iElem){
    globalId = oneDimensionalMap->GID(iElem);
    getExodusMeshNodePositions(globalId, nodeCoordinates);
    numNodes = nodeCoordinates.size()/3;
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numNodes != 8 &&
                                numNodes != 4 &&
                                numNodes != 1 &&
                                numNodes != 10 &&
                                numNodes != 20,
                                "**** Invalid element topology.\n");
    if(numNodes == 1){
      // Sphere element
      volume = (*cellVolume)[iElem];
      length = std::pow( (3.0*volume)/(4.0*pi), 1.0/3.0 );
    }
    else{

      // Consider 4- and 10-noded tets (treat either as 4-noded)
      // Consider 10- and 20-noded hexes (treat either as 10-noded)
      numNodesInElement = 4;
      if(numNodes == 8 || numNodes == 20)
        numNodesInElement = 8;

      for(int i=0 ; i<numNodesInElement ; ++i){
        for(int j=i ; j<numNodesInElement ; ++j){
          if(i != j){
            x = nodeCoordinates[i*3]   - nodeCoordinates[j*3];
            y = nodeCoordinates[i*3+1] - nodeCoordinates[j*3+1];
            z = nodeCoordinates[i*3+2] - nodeCoordinates[j*3+2];
            length = x*x + y*y + z*z;
            if(length > localMaxElementDimension)
              localMaxElementDimension = length;
          }
        }
      }
    }
  }

  localMaxElementDimension = std::sqrt(localMaxElementDimension);

  // Parallel communication to determine overall maximum element dimension
  double globalMaxElementDimension;
  comm->MaxAll(&localMaxElementDimension, &globalMaxElementDimension, 1);

  return globalMaxElementDimension;
}

void PeridigmNS::STKDiscretization::removeNonintersectingNeighborsFromNeighborList(Teuchos::RCP<Epetra_Vector> x,
                                                                                   Teuchos::RCP<Epetra_Vector> searchRadii,
                                                                                   Teuchos::RCP<Epetra_BlockMap>& overlapMap,
                                                                                   int& neighborListSize,
                                                                                   int*& neighborList)
{
  int refinedNumNeighbors, numNeighbors, neighborLocalId, neighborGlobalId;
  unsigned int refinedNumNeighborsIndex;
  double horizon;
  SphereIntersection sphereIntersection;
  vector<double> exodusNodePositions;
  vector<double> sphereCenter(3);
  vector<int> refinedNeighborGlobalIdList;
  refinedNeighborGlobalIdList.reserve(neighborListSize);

  int index = 0;
  int elemLocalId = 0;
  while(index < neighborListSize){
    numNeighbors = neighborList[index++];
    refinedNumNeighborsIndex = refinedNeighborGlobalIdList.size();
    refinedNumNeighbors = 0;
    refinedNeighborGlobalIdList.push_back(refinedNumNeighbors);
    for(int iNeighbor=0 ; iNeighbor<numNeighbors ; ++iNeighbor){
      neighborLocalId = neighborList[index++];
      neighborGlobalId = overlapMap->GID(neighborLocalId);

      // Determine if the element intersects the sphere
      sphereCenter[0] = (*x)[elemLocalId];
      sphereCenter[1] = (*x)[elemLocalId+1];
      sphereCenter[2] = (*x)[elemLocalId+2];
      getExodusMeshNodePositions(neighborGlobalId, exodusNodePositions);
      horizon = (*searchRadii)[elemLocalId];
   
      TEUCHOS_TEST_FOR_EXCEPT_MSG(exodusNodePositions.size()/3 != 8,
                                  "\n**** Error:  Element-horizon intersection calculations currently enabled only for hexahedron elements.\n");

      sphereIntersection = hexahedronSphereIntersection(&exodusNodePositions[0], sphereCenter, horizon);

      if(sphereIntersection != OUTSIDE_SPHERE){
        refinedNeighborGlobalIdList.push_back(neighborGlobalId);
        refinedNumNeighbors += 1;
      }
    }
    refinedNeighborGlobalIdList[refinedNumNeighborsIndex] = refinedNumNeighbors;
    elemLocalId += 1;
    cout << "DEBUGGING numNeighbors = " << numNeighbors << ", refinedNumNeighbors = " << refinedNumNeighbors << endl;
  }

  // TODO
  // Debug above code (way too many neighbors are being excluded, s.f. WaveInBar test)
  // Create new overlap map and neighborlist based on refinedNeighborGlobalIdList
}

double PeridigmNS::STKDiscretization::scalarTripleProduct(std::vector<double>& a,
                                                          std::vector<double>& b,
                                                          std::vector<double>& c) const
{
  double tripleProduct = 
    a[0]*(b[1]*c[2] - b[2]*c[1]) + a[1]*(b[2]*c[0] - b[0]*c[2]) + a[2]*(b[0]*c[1] - b[1]*c[0]);

  return tripleProduct;
}

double PeridigmNS::STKDiscretization::hexVolume(std::vector<double*>& nodeCoordinates) const
{
  int map[] = {0,1,3,2,4,5,7,6};
  std::vector<double> x17(3), x27(3), x47(3), x06(3), x05(3), x03(3), A(3), B(3), C(3);  
  for(int dof=0 ; dof<3 ; ++dof){
    x17[dof] = nodeCoordinates[map[7]][dof] - nodeCoordinates[map[1]][dof];
    x27[dof] = nodeCoordinates[map[7]][dof] - nodeCoordinates[map[2]][dof];
    x47[dof] = nodeCoordinates[map[7]][dof] - nodeCoordinates[map[4]][dof];
    x06[dof] = nodeCoordinates[map[6]][dof] - nodeCoordinates[map[0]][dof];
    x05[dof] = nodeCoordinates[map[5]][dof] - nodeCoordinates[map[0]][dof];
    x03[dof] = nodeCoordinates[map[3]][dof] - nodeCoordinates[map[0]][dof];
  }

  for(int dof=0 ; dof<3 ; ++dof){
    A[dof] = x17[dof] + x06[dof];
    B[dof] = x27[dof] + x05[dof];
    C[dof] = x47[dof] + x03[dof];
  }

  double v1 = fabs( scalarTripleProduct(A,x27,x03) );
  double v2 = fabs( scalarTripleProduct(x06,B,x47) );
  double v3 = fabs( scalarTripleProduct(x17,x05,C) );

  return (v1+v2+v3)/12.0;
}

double PeridigmNS::STKDiscretization::tetVolume(std::vector<double*>& nodeCoordinates) const
{
  // Change the coordinate system such that the first point is at the origin.
  // The other three points are labeled a, b, and c.
  std::vector<double> a(3), b(3), c(3);
  for(int dof=0 ; dof<3 ; ++dof){
    a[dof] = nodeCoordinates[1][dof] - nodeCoordinates[0][dof];
    b[dof] = nodeCoordinates[2][dof] - nodeCoordinates[0][dof];
    c[dof] = nodeCoordinates[3][dof] - nodeCoordinates[0][dof];
  }

  // The volume is then | a . (b x c) | / 6
  double volume = scalarTripleProduct(a, b, c) / 6.0;

  return volume;
}

// double PeridigmNS::STKDiscretization::hexMaxElementDimension(std::vector<double*>& nodeCoordinates) const
// {
//   double maxDimension = 0.0;
//   double dx, dy, dz, diagonal;

//   // Check edges

//   // Check face diagonals

//   // Check element diagonals
  
//   // Exodus nodes 1 7
//   dx = nodeCoordinates[0][0] - nodeCoordinates[6][0];
//   dy = nodeCoordinates[0][1] - nodeCoordinates[6][1];
//   dz = nodeCoordinates[0][2] - nodeCoordinates[6][2];
//   diagonal = sqrt(dx*dx + dy*dy + dz*dz);
//   if(diagonal > maxDimension)
//     maxDimension = diagonal;

//   // Exodus nodes 2 8
//   dx = nodeCoordinates[1][0] - nodeCoordinates[7][0];
//   dy = nodeCoordinates[1][1] - nodeCoordinates[7][1];
//   dz = nodeCoordinates[1][2] - nodeCoordinates[7][2];
//   diagonal = sqrt(dx*dx + dy*dy + dz*dz);
//   if(diagonal > maxDimension)
//     maxDimension = diagonal;

//   // Exodus nodes 3 5
//   dx = nodeCoordinates[2][0] - nodeCoordinates[4][0];
//   dy = nodeCoordinates[2][1] - nodeCoordinates[4][1];
//   dz = nodeCoordinates[2][2] - nodeCoordinates[4][2];
//   diagonal = sqrt(dx*dx + dy*dy + dz*dz);
//   if(diagonal > maxDimension)
//     maxDimension = diagonal;

//   // Exodus nodes 4 6
//   dx = nodeCoordinates[3][0] - nodeCoordinates[5][0];
//   dy = nodeCoordinates[3][1] - nodeCoordinates[5][1];
//   dz = nodeCoordinates[3][2] - nodeCoordinates[5][2];
//   diagonal = sqrt(dx*dx + dy*dy + dz*dz);
//   if(diagonal > maxDimension)
//     maxDimension = diagonal;

//   return maxDimension;
// }

// double PeridigmNS::STKDiscretization::tetMaxElementDimension(std::vector<double*>& nodeCoordinates) const
// {
//   double maxDimension = 0.0;
//   double dx, dy, dz, edgeLength;

//   // Exodus nodes 1 2
//   dx = nodeCoordinates[0][0] - nodeCoordinates[1][0];
//   dy = nodeCoordinates[0][1] - nodeCoordinates[1][1];
//   dz = nodeCoordinates[0][2] - nodeCoordinates[1][2];
//   edgeLength = sqrt(dx*dx + dy*dy + dz*dz);
//   if(edgeLength > maxDimension)
//     maxDimension = edgeLength;

//   // Exodus nodes 1 3
//   dx = nodeCoordinates[0][0] - nodeCoordinates[2][0];
//   dy = nodeCoordinates[0][1] - nodeCoordinates[2][1];
//   dz = nodeCoordinates[0][2] - nodeCoordinates[2][2];
//   edgeLength = sqrt(dx*dx + dy*dy + dz*dz);
//   if(edgeLength > maxDimension)
//     maxDimension = edgeLength;

//   // Exodus nodes 1 4
//   dx = nodeCoordinates[0][0] - nodeCoordinates[3][0];
//   dy = nodeCoordinates[0][1] - nodeCoordinates[3][1];
//   dz = nodeCoordinates[0][2] - nodeCoordinates[3][2];
//   edgeLength = sqrt(dx*dx + dy*dy + dz*dz);
//   if(edgeLength > maxDimension)
//     maxDimension = edgeLength;

//   // Exodus nodes 2 3
//   dx = nodeCoordinates[1][0] - nodeCoordinates[2][0];
//   dy = nodeCoordinates[1][1] - nodeCoordinates[2][1];
//   dz = nodeCoordinates[1][2] - nodeCoordinates[2][2];
//   edgeLength = sqrt(dx*dx + dy*dy + dz*dz);
//   if(edgeLength > maxDimension)
//     maxDimension = edgeLength;

//   // Exodus nodes 2 4
//   dx = nodeCoordinates[1][0] - nodeCoordinates[3][0];
//   dy = nodeCoordinates[1][1] - nodeCoordinates[3][1];
//   dz = nodeCoordinates[1][2] - nodeCoordinates[3][2];
//   edgeLength = sqrt(dx*dx + dy*dy + dz*dz);
//   if(edgeLength > maxDimension)
//     maxDimension = edgeLength;

//   // Exodus nodes 3 4
//   dx = nodeCoordinates[2][0] - nodeCoordinates[3][0];
//   dy = nodeCoordinates[2][1] - nodeCoordinates[3][1];
//   dz = nodeCoordinates[2][2] - nodeCoordinates[3][2];
//   edgeLength = sqrt(dx*dx + dy*dy + dz*dz);
//   if(edgeLength > maxDimension)
//     maxDimension = edgeLength;

//   return maxDimension;
// }
