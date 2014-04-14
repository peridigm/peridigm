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

//#define DEBUGGING_BACKWARDS_COMPATIBILITY_NEIGHBORHOOD_LIST

PeridigmNS::STKDiscretization::STKDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) :
  minElementRadius(1.0e50),
  maxElementRadius(0.0),
  storeExodusMesh(false),
  computeIntersections(false),
  constructInterfaces(false),
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
  if(params->isParameter("Store Exodus Mesh")){
    storeExodusMesh = params->get<bool>("Store Exodus Mesh");
  }
  if(params->isParameter("Compute Element-Horizon Intersections")){
    computeIntersections = params->get<bool>("Compute Element-Horizon Intersections");
    storeExodusMesh = true;
  }
  // also store the mesh if interfaces will be used (i.e. for the aperture calcs)
  if(params->isParameter("Construct Interfaces")){
    constructInterfaces = params->get<bool>("Construct Interfaces");
    storeExodusMesh = constructInterfaces;
  }

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

  int neighborListSize;
  int* neighborList;

  // Execute the neighbor search
  // When computing element-horizon intersections, the search is expanded by the maximum element dimension
  if(computeIntersections)
    ProximitySearch::GlobalProximitySearch(initialX, horizonForEachPoint, oneDimensionalOverlapMap, neighborListSize, neighborList, bondFilters, maxElementDimension);
  else
    ProximitySearch::GlobalProximitySearch(initialX, horizonForEachPoint, oneDimensionalOverlapMap, neighborListSize, neighborList, bondFilters);

  // Ghost exodus data so that element-horizon intersections can be calculated for ghosted neighbors
  if(storeExodusMesh)
    ghostExodusMeshData();

  // Remove elements from neighbor lists that are outside the horizon
  // Some will have been picked up in the initial neighbor search when computing element-horizon intersections
  if(computeIntersections)
    removeNonintersectingNeighborsFromNeighborList(initialX, horizonForEachPoint, oneDimensionalMap, oneDimensionalOverlapMap, neighborListSize, neighborList);

  createNeighborhoodData(neighborListSize, neighborList);

  // if interfaces are requested construct the interfaces after the neighborhood data is known:
  if(constructInterfaces)
    constructInterfaceData();



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

PeridigmNS::STKDiscretization::~STKDiscretization() {
}

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
  // Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  // int localElemCount = elements.size();
  // int globalElemCount(0);
  // reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, int(1), &localElemCount, &globalElemCount);

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
    int elementId = elements[iElem]->identifier() - 1;
    globalElementIds[iElem] = elementId;
  }

  // Create a vector for storing exodus element connectivity
  if(storeExodusMesh){
    int numGlobalElements = -1;
    int numMyElements = static_cast<int>(globalElementIds.size());
    vector<int> myGlobalElements(numMyElements);
    vector<int> elementSizeList(numMyElements);
    int indexBase = 0;

    for(unsigned int iElem=0 ; iElem<elements.size() ; ++iElem){
      stk::mesh::PairIterRelation nodeRelations = elements[iElem]->node_relations();
      int elementId = elements[iElem]->identifier() - 1;
      int numNodesInElement = nodeRelations.size();
      myGlobalElements[iElem] = elementId;
      elementSizeList[iElem] = numNodesInElement;
    }

    // Epetra_BlockMap exodusMeshNodePositionsMap(-1, nodes.size(), &myGlobalNodeIds[0], 3, 0, *comm);
    //  bondMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, *comm));
    Epetra_BlockMap exodusMeshElementConnectivityMap(numGlobalElements, numMyElements, &myGlobalElements[0], &elementSizeList[0], indexBase, *comm);
    exodusMeshElementConnectivity = Teuchos::rcp(new Epetra_Vector(exodusMeshElementConnectivityMap));
    exodusMeshElementConnectivity->PutScalar(-1.0);
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
    int elementId = elements[iElem]->identifier() - 1;
    initialXPtr[iElem*3] = 0.0;
    initialXPtr[iElem*3+1] = 0.0;
    initialXPtr[iElem*3+2] = 0.0;
    std::vector<double*> nodeCoordinates(nodeRelations.size());
    int iNode = 0;

    int exodusMeshElementIndex = 0;
    if(storeExodusMesh){
      int localExodusElementId = exodusMeshElementConnectivity->Map().LID(elementId);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(localExodusElementId == -1, "**** Invalid local Exodus element Id.\n");
      exodusMeshElementIndex = exodusMeshElementConnectivity->Map().FirstPointInElement(localExodusElementId);
    }

    for(stk::mesh::PairIterRelation::iterator it=nodeRelations.begin() ; it!=nodeRelations.end() ; ++it){
      stk::mesh::Entity* node = it->entity();
      if(storeExodusMesh){
        (*exodusMeshElementConnectivity)[exodusMeshElementIndex++] = node->identifier() - 1;
      }
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
      cellVolumePtr[iElem] = tetVolume(nodeCoordinates);
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
    const std::vector<int>& elementIds = it->second;
    for(unsigned int i=0 ; i<elementIds.size() ; ++i){
      int globalID = elementIds[i];
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
PeridigmNS::STKDiscretization::constructInterfaceData()
{
  int* const neighPtr = neighborhoodData->NeighborhoodList();
  const int listSize = neighborhoodData->NeighborhoodListSize();

  // faces of a cube as stored by exodus ordering:
  int faces[6][4] = {{4,5,6,7},{7,6,2,3},{0,4,7,3},{0,1,2,3},{0,4,5,1},{1,5,6,2}};

  interfaceData = Teuchos::rcp(new PeridigmNS::InterfaceData);

  TEUCHOS_TEST_FOR_EXCEPTION(storeExodusMesh!=true,std::logic_error," Exodus mesh should have been stored if this is called.");
  std::vector<int> leftElements;
  std::vector<int> rightElements;
  std::vector<int> numNodes;
  std::vector<std::vector<int> > interfaceNodesVec;

  int elemIndex = 0;
  for(int i=0 ; i<listSize; i++){
    const int selfGID = oneDimensionalMap->GID(elemIndex);
    int numNeighbors = neighPtr[i];
    const int numNodesPerElem = exodusMeshElementConnectivity->Map().ElementSize(elemIndex);
    const int myIndex = exodusMeshElementConnectivity->Map().FirstPointInElement(elemIndex);
    // get the connectivity for this element:
    std::vector<int> selfNodeIds(numNodesPerElem);
    for(unsigned n=0;n<numNodesPerElem;++n){
      selfNodeIds[n] = (*exodusMeshElementConnectivity)[myIndex+n];
    }

    for(unsigned j=0;j<numNeighbors;j++){
      i++;
      int numNodesFound = 0;
      const int GID = oneDimensionalOverlapMap->GID(neighPtr[i]);
      const int neighIndex = exodusMeshElementConnectivity->Map().FirstPointInElement(neighPtr[i]);

      bool shareFace = false;
      // the two elements must have a face in common to be an interface pair:
      std::vector<int> neighNodeIds(numNodesPerElem);
      std::vector<int> foundNodeIds(4); // 4 works for tris or quads
      for(unsigned n=0;n<numNodesPerElem;++n)
        neighNodeIds[n] = (*exodusMeshElementConnectivity)[neighIndex + n];

      if(numNodesPerElem==8){ // cube (order needed to prevent folding the interfaces by mixing up the node numbering)
        for(unsigned ni=0;ni<6;++ni){
          numNodesFound = 0;
          for(unsigned nj=0;nj<4;++nj){
            for(unsigned nk=0;nk<numNodesPerElem;++nk){
              if(selfNodeIds[faces[ni][nj]]==neighNodeIds[nk]){
                numNodesFound++;
              }
            }
          }
          if(numNodesFound>2){ // this is a shared face
            for(unsigned kn=0;kn<4;++kn)
              foundNodeIds[kn] = selfNodeIds[faces[ni][kn]];
            if((GID!=-1&&GID>selfGID)||(selfGID<GID)){
              leftElements.push_back(selfGID);
              rightElements.push_back(GID);
              numNodes.push_back(numNodesFound);
              interfaceNodesVec.push_back(foundNodeIds);
              break;
            }
          }
        }
      }
      else { // tets (order doesn't matter)
        for(unsigned ni=0;ni<numNodesPerElem;++ni){
          for(unsigned nj=0;nj<numNodesPerElem;++nj){
            if(selfNodeIds[ni]==neighNodeIds[nj]){
              foundNodeIds[numNodesFound] = selfNodeIds[ni];
              numNodesFound++;
            }
          }
        }
        shareFace = numNodesFound>2; // works for tri and quad
        if(((GID!=-1&&GID>selfGID)||(selfGID<GID)) && shareFace){
          leftElements.push_back(selfGID);
          rightElements.push_back(GID);
          numNodes.push_back(numNodesFound);
          interfaceNodesVec.push_back(foundNodeIds);
        }
      }
    }
    elemIndex++;
  }
  // now create the interface data and populate it
  interfaceData->Initialize(leftElements, rightElements, numNodes, interfaceNodesVec, comm);
  // generate an exodus file for output:
  interfaceData->InitializeExodusOutput(exodusMeshElementConnectivity,exodusMeshNodePositions);
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

  int localId = exodusMeshElementConnectivity->Map().LID(globalNodeID);
  TEUCHOS_TEST_FOR_EXCEPT_MSG(localId == -1, "**** Invalid local Exodus element Id in getExodusMeshNodePositions().\n");
  unsigned int numNodes = exodusMeshElementConnectivity->Map().ElementSize(localId);
  int exodusMeshElementIndex = exodusMeshElementConnectivity->Map().FirstPointInElement(localId);
  vector<int> elementConnectivity(numNodes);
  for(unsigned int i=0 ; i<numNodes ; ++i)
    elementConnectivity[i] = static_cast<int>( (*exodusMeshElementConnectivity)[exodusMeshElementIndex++] );

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
                                                                                   Teuchos::RCP<Epetra_BlockMap> ownedMap,
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
  set<int> refinedGlobalIds;
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
      sphereCenter[0] = (*x)[3*elemLocalId];
      sphereCenter[1] = (*x)[3*elemLocalId+1];
      sphereCenter[2] = (*x)[3*elemLocalId+2];
      getExodusMeshNodePositions(neighborGlobalId, exodusNodePositions);
      horizon = (*searchRadii)[elemLocalId];
   
      TEUCHOS_TEST_FOR_EXCEPT_MSG(exodusNodePositions.size()/3 != 8,
                                  "\n**** Error:  Element-horizon intersection calculations currently enabled only for hexahedron elements.\n");

#ifdef DEBUGGING_BACKWARDS_COMPATIBILITY_NEIGHBORHOOD_LIST
      double centroid[3];
      centroid[0] = 0.0;
      centroid[1] = 0.0;
      centroid[2] = 0.0;
      for(int i=0 ; i<8 ; ++i){
        centroid[0] += exodusNodePositions[3*i];
        centroid[1] += exodusNodePositions[3*i+1];
        centroid[2] += exodusNodePositions[3*i+2];
      }
      centroid[0] /= 8.0;
      centroid[1] /= 8.0;
      centroid[2] /= 8.0;

      double distanceSquared = (sphereCenter[0] - centroid[0])*(sphereCenter[0] - centroid[0])
        + (sphereCenter[1] - centroid[1])*(sphereCenter[1] - centroid[1])
        + (sphereCenter[2] - centroid[2])*(sphereCenter[2] - centroid[2]);

      if(distanceSquared > horizon*horizon)
        sphereIntersection = OUTSIDE_SPHERE;
      else
        sphereIntersection = INSIDE_SPHERE;
#else
      sphereIntersection = hexahedronSphereIntersection(&exodusNodePositions[0], sphereCenter, horizon);
#endif

      if(sphereIntersection != OUTSIDE_SPHERE){
        refinedNeighborGlobalIdList.push_back(neighborGlobalId);
        refinedGlobalIds.insert(neighborGlobalId);
        refinedNumNeighbors += 1;
      }
    }
    refinedNeighborGlobalIdList[refinedNumNeighborsIndex] = refinedNumNeighbors;
    elemLocalId += 1;
  }

  // Create new overlap map and neighborlist based on refinedNeighborGlobalIdList
  vector<int> refinedGlobalIdVector;
  refinedGlobalIdVector.reserve(refinedGlobalIds.size());
  // The non-ghost portion of the overlapMap must match the ownedMap
  for(int i=0 ; i<ownedMap->NumMyElements() ; ++i)
    refinedGlobalIdVector.push_back(ownedMap->GID(i));
  // Add the ghosts to the overlapMap
  for(set<int>::const_iterator it=refinedGlobalIds.begin() ; it!=refinedGlobalIds.end() ; it++){
    if(!ownedMap->MyGID(*it))
      refinedGlobalIdVector.push_back(*it);
  }
  overlapMap = Teuchos::rcp(new Epetra_BlockMap(-1,
                                                static_cast<int>( refinedGlobalIdVector.size() ),
                                                &refinedGlobalIdVector[0],
                                                1,
                                                0,
                                                *comm));
  neighborListSize = static_cast<int>(refinedNeighborGlobalIdList.size());
  delete[] neighborList;
  neighborList = new int[neighborListSize];
  index = 0;
  while(index < neighborListSize){
    numNeighbors = refinedNeighborGlobalIdList[index];
    neighborList[index] = numNeighbors;
    index += 1;
    for(int i=0 ; i<numNeighbors ; ++i){
      neighborLocalId = overlapMap->LID(refinedNeighborGlobalIdList[index]);
      TEUCHOS_TEST_FOR_EXCEPT_MSG(neighborLocalId == -1, "\n**** Error:  Invalid local ID in removeNonintersectingNeighborsFromNeighborList().\n");
      neighborList[index] = neighborLocalId;
      index += 1;
    }
  }
}

void PeridigmNS::STKDiscretization::ghostExodusMeshData()
{
  int numGlobalElements = -1;
  int numMyElements = oneDimensionalOverlapMap->NumMyElements();
  int* myGlobalElements = oneDimensionalOverlapMap->MyGlobalElements();
  vector<int> elementSizeList(numMyElements);
  int indexBase = 0;

  // Obtain the element sizes for off-processor elements that will be ghosted
  Epetra_Vector elementSize(*oneDimensionalMap);
  for(int i=0 ; i<elementSize.MyLength() ; ++i)
    elementSize[i] = exodusMeshElementConnectivity->Map().ElementSize(i);
  Epetra_Vector elementSizeOverlap(*oneDimensionalOverlapMap);
  Epetra_Import elementSizeImporter(*oneDimensionalOverlapMap, *oneDimensionalMap);
  elementSizeOverlap.Import(elementSize, elementSizeImporter, Insert);
  for(int i=0 ; i<elementSizeOverlap.MyLength() ; ++i)
    elementSizeList[i] = static_cast<int>(elementSizeOverlap[i]);

  Epetra_BlockMap overlapExodusMeshElementConnectivityMap(numGlobalElements, numMyElements, &myGlobalElements[0], &elementSizeList[0], indexBase, *comm);

  // Import the element connectivities into the overlap vector
  Teuchos::RCP<Epetra_Vector> overlapExodusMeshElementConnectivity = Teuchos::rcp(new Epetra_Vector(overlapExodusMeshElementConnectivityMap));
  overlapExodusMeshElementConnectivity->PutScalar(-1.0);
  Epetra_Import elementConnectivityImporter(overlapExodusMeshElementConnectivityMap, exodusMeshElementConnectivity->Map());
  overlapExodusMeshElementConnectivity->Import(*exodusMeshElementConnectivity, elementConnectivityImporter, Insert);

  // Create a list of global Exodus node ids (including ghosts)
  set<int> globalNodeIdsSet;
  for(int i=0 ; i<overlapExodusMeshElementConnectivity->MyLength() ; ++i)
    globalNodeIdsSet.insert( static_cast<int>( (*overlapExodusMeshElementConnectivity)[i] ) );
  vector<int> globalNodeIds(globalNodeIdsSet.size());
  int index = 0;
  for(set<int>::const_iterator it=globalNodeIdsSet.begin() ; it!=globalNodeIdsSet.end() ; it++)
    globalNodeIds[index++] = *it;
  std::sort(globalNodeIds.begin(), globalNodeIds.end());

  Epetra_BlockMap overlapExodusMeshNodePositionsMap(numGlobalElements, globalNodeIds.size(), &globalNodeIds[0], 3, indexBase, *comm);
  Teuchos::RCP<Epetra_Vector> overlapExodusMeshNodePositions = Teuchos::rcp(new Epetra_Vector(overlapExodusMeshNodePositionsMap));
  Epetra_Import exodusMeshNodePositionsImporter(overlapExodusMeshNodePositionsMap, exodusMeshNodePositions->Map());
  overlapExodusMeshNodePositions->Import(*exodusMeshNodePositions, exodusMeshNodePositionsImporter, Insert);

  // Set the exodus node positions vector and the connectivity vector to the new overlap vectors
  exodusMeshNodePositions = overlapExodusMeshNodePositions;
  exodusMeshElementConnectivity = overlapExodusMeshElementConnectivity;
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
