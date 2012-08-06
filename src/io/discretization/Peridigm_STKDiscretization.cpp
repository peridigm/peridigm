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
#include "pdneigh/NeighborhoodList.h"
#include "pdneigh/BondFilter.h"
#include "pdneigh/PdZoltan.h"

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

using namespace std;

PeridigmNS::STKDiscretization::STKDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) :
  minElementRadius(1.0e50),
  maxElementDimension(0.0),
  numBonds(0),
  myPID(epetra_comm->MyPID()),
  numPID(epetra_comm->NumProc()),
  bondFilterCommand("None"),
  comm(epetra_comm)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Exodus", "Invalid Type in STKDiscretization");

  horizon = params->get<double>("Horizon");
  searchHorizon = params->get<double>("Search Horizon");
  if(params->isParameter("Omit Bonds Between Blocks"))
    bondFilterCommand = params->get<string>("Omit Bonds Between Blocks");
  string meshFileName = params->get<string>("Input Mesh File");

  QUICKGRID::Data decomp = getDecomp(meshFileName, searchHorizon);

  // \todo Refactor; the createMaps() call is currently inside getDecomp() due to order-of-operations issues with tracking element blocks.
  //createMaps(decomp);
  createNeighborhoodData(decomp);

  // \todo Move this functionality to base class, it's currently duplicated in PdQuickGridDiscretization.
  // Create the bondMap, a local map used for constitutive data stored on bonds.
  // Due to Epetra_BlockMap restrictions, there can not be any entries with length zero.
  // This means that points with no neighbors can not appear in the bondMap.
  int numMyElementsUpperBound = oneDimensionalMap->NumMyElements();
  int numGlobalElements = -1; 
  int numMyElements = 0;
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
    neighborhoodIndex += 1 + numNeighbors;
  }
  int indexBase = 0;
  bondMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, *comm));
  delete[] myGlobalElements;
  delete[] elementSizeList;

  // 3D only
  TEUCHOS_TEST_FOR_EXCEPT_MSG(decomp.dimension != 3, "Invalid dimension in decomposition (only 3D is supported)");

  // fill the x vector with the current positions (owned positions only)
  initialX = Teuchos::rcp(new Epetra_Vector(Copy,*threeDimensionalMap,decomp.myX.get()) );

  // fill cell volumes
  cellVolume = Teuchos::rcp(new Epetra_Vector(Copy,*oneDimensionalMap,decomp.cellVolume.get()) );

  // find the minimum element radius
  for(int i=0 ; i<cellVolume->MyLength() ; ++i){
    double radius = pow(0.238732414637843*(*cellVolume)[i], 0.33333333333333333);
    if(radius < minElementRadius)
      minElementRadius = radius;
  }
}

PeridigmNS::STKDiscretization::~STKDiscretization() {}


QUICKGRID::Data PeridigmNS::STKDiscretization::getDecomp(const string& meshFileName,
                                                         double horizon) {

  string meshType = "exodusii";
  string workingDirectory = "";
  Teuchos::RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(comm, true);
  metaData = Teuchos::rcp(new stk::mesh::fem::FEMMetaData);
  // \todo Remove backwards compatibility after next Trilinos release
  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  meshData = Teuchos::rcp(new stk::io::MeshData);
  #else
  meshData = Teuchos::rcp(new stk::io::util::MeshData);
  #endif
  Ioss::Init::Initializer io;
  // \todo Remove backwards compatibility after next Trilinos release
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

  // Loop over the parts and store the element parts, side parts, and node parts.
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
  // \todo Remove backwards compatibility after next Trilinos release
  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  stk::io::populate_bulk_data(bulkData, *meshData);
  #else
  stk::io::util::populate_bulk_data(bulkData, *meshData, "exodusii");
  #endif
  bulkData.modification_end();

  stk::mesh::Field<double, stk::mesh::Cartesian>* coordinatesField = 
    metaData->get_field< stk::mesh::Field<double, stk::mesh::Cartesian> >("coordinates");

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

  // Copy data from stk into a decomp object
  int myNumElements = elements.size();
  int dimension = 3;
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(myNumElements, dimension);
  decomp.globalNumPoints = globalElemCount;
  int* globalIds = decomp.myGlobalIDs.get();
  double* volumes = decomp.cellVolume.get();
  double* centroids = decomp.myX.get();

  // Store the positions of the nodes in the initial mesh
  // This is used for the calculation of partial neighbor volumes
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

  // loop over the elements and fill the decomp data structure
  for(unsigned int iElem=0 ; iElem<elements.size() ; ++iElem){
    stk::mesh::PairIterRelation nodeRelations = elements[iElem]->node_relations();
    int elementID = elements[iElem]->identifier() - 1;
    exodusMeshElementConnectivity[elementID] = vector<int>();
    centroids[iElem*3] = 0.0;
    centroids[iElem*3+1] = 0.0;
    centroids[iElem*3+2] = 0.0;
    std::vector<double*> nodeCoordinates(nodeRelations.size());
    int iNode = 0;
    for(stk::mesh::PairIterRelation::iterator it=nodeRelations.begin() ; it!=nodeRelations.end() ; ++it){
      stk::mesh::Entity* node = it->entity();
      exodusMeshElementConnectivity[elementID].push_back(node->identifier() - 1);
      double* coordinates = stk::mesh::field_data(*coordinatesField, *node);
      centroids[iElem*3] += coordinates[0];
      centroids[iElem*3+1] += coordinates[1];
      centroids[iElem*3+2] += coordinates[2];
      nodeCoordinates[iNode++] = coordinates;
    }
    centroids[iElem*3] /= nodeRelations.size();
    centroids[iElem*3+1] /= nodeRelations.size();
    centroids[iElem*3+2] /= nodeRelations.size();
    volumes[iElem] = hexVolume(nodeCoordinates);
    globalIds[iElem] = elements[iElem]->identifier() - 1; 
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
    std::vector<int> elementGlobalIds;
    for(unsigned int iElement=0 ; iElement<elementsInElementBlock.size() ; iElement++)
      elementGlobalIds.push_back(elementsInElementBlock[iElement]->identifier() - 1);

    // \todo Just load elementBlock directly.
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

    // Create a selector for all locally-owned nodes in the node set
    stk::mesh::Selector selector = 
      stk::mesh::Selector( *stkNodeSets[iNodeSet] ) & stk::mesh::Selector( metaData->locally_owned_part() );

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
        elementGlobalIds.insert(globalId);
      }
    }

    //
    //  \todo TEMPORARY SOLUTION THAT WILL NOT SCALE.
    //
    // Get the full node list on every processor.
    // This is important because we're about to rebalance, and after we do so nodes
    // that are locally owned may not be locally owned.  If a node that
    // changes processor is in a node set, this will cause problems downstream
    // when assigning initial/boundary conditions.  A simple solution is to
    // have all processors know about all the nodes in each node set.  But this
    // won't scale.
    // 
    // \todo Figure out why it doesn't work to simpily use selector = stk::mesh::Selector( *stkNodeSets[iNodeSet] ) to get the complete node set.
    //
    vector<int> localNodeSetLength(numPID, 0);
    vector<int> tempLocalNodeSetLength(numPID, 0);
    tempLocalNodeSetLength[myPID] = (int)elementGlobalIds.size();
    reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, (int)numPID, &tempLocalNodeSetLength[0], &localNodeSetLength[0]);

    int totalNodeSetLength = 0;
    int offset = 0;
    for(unsigned int i=0 ; i<localNodeSetLength.size() ; ++i){
      if(i < myPID)
        offset += localNodeSetLength[i];
      totalNodeSetLength += localNodeSetLength[i];
    }
    vector<int> completeElementGlobalIds(totalNodeSetLength, 0);
    vector<int> tempCompleteElementGlobalIds(totalNodeSetLength, 0);
    int index = 0;
    for(std::set<int>::iterator it=elementGlobalIds.begin() ; it!=elementGlobalIds.end() ; it++){
      tempCompleteElementGlobalIds[index+offset] = *it;
      index++;
    }
    reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, totalNodeSetLength, &tempCompleteElementGlobalIds[0], &completeElementGlobalIds[0]);

    // Load the node set into the nodeSets container
    for(unsigned int i=0 ; i<completeElementGlobalIds.size() ; ++i)
      nodeSet.push_back( completeElementGlobalIds[i] );
  }

  // Create a blockID vector in the current configuration
  // That is, the configuration prior to load balancing
  Epetra_BlockMap tempOneDimensionalMap(decomp.globalNumPoints,
                                        decomp.numPoints,
                                        decomp.myGlobalIDs.get(),
                                        1,
                                        0,
                                        *comm);
  Epetra_Vector tempBlockID(tempOneDimensionalMap);
  double* tempBlockIDPtr;
  tempBlockID.ExtractView(&tempBlockIDPtr);
  std::map< std::string, std::vector<int> >::const_iterator it;
  for(it = elementBlocks->begin() ; it != elementBlocks->end() ; it++){
    const std::string& blockName = it->first;

    size_t loc = blockName.find_last_of('_');
    TEUCHOS_TEST_FOR_EXCEPT_MSG(loc == string::npos, "\n**** Parse error, invalid block name.\n");
    stringstream blockIDSS(blockName.substr(loc+1, blockName.size()));
    int blockID;
    blockIDSS >> blockID;

    const std::vector<int>& elementIDs = it->second;
    for(unsigned int i=0 ; i<elementIDs.size() ; ++i){
      int globalID = elementIDs[i];
      int localID = tempOneDimensionalMap.LID(globalID);
      tempBlockIDPtr[localID] = blockID;
    }
  }

  // \todo Remove backwards compatibility after next Trilinos release
  #if TRILINOS_MAJOR_MINOR_VERSION > 101002
  // free the meshData
  meshData = Teuchos::RCP<stk::io::MeshData>();
  #else
  // free the meshData
  meshData = Teuchos::RCP<stk::io::util::MeshData>();
  #endif

  // call the rebalance function on the current-configuration decomp
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);
  
  // execute neighbor search and update the decomp to include resulting ghosts
  shared_ptr<PdBondFilter::BondFilter> bondFilterPtr(new PdBondFilter::BondFilterDefault(false));
  shared_ptr<const Epetra_Comm> commSp(comm.getRawPtr(),NonDeleter<const Epetra_Comm>());
  PDNEIGH::NeighborhoodList list(commSp,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,horizon,bondFilterPtr);
  decomp.neighborhood=list.get_neighborhood();
  decomp.sizeNeighborhoodList=list.get_size_neighborhood_list();
  decomp.neighborhoodPtr=list.get_neighborhood_ptr();

  // Create all the maps.
  // \todo This call really should be outside this function, but we need the maps here; should somehow remove the element id stuff from this function.
  createMaps(decomp);

  // Create a blockID vector corresponding to the load balanced decomposition
  // \todo Refactor to avoid duplication creation of 1D map, order of operations is currently jumbled
  blockID = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  Epetra_Import tempImporter(blockID->Map(), tempBlockID.Map());
  blockID->Import(tempBlockID, tempImporter, Insert);

  return decomp;
}

void
PeridigmNS::STKDiscretization::createMaps(const QUICKGRID::Data& decomp)
{
  int dimension;

  // oneDimensionalMap
  // used for global IDs and scalar data
  dimension = 1;
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(AbstractDiscretization::getOwnedMap(*comm, decomp, dimension)));

  // oneDimensionalOverlapMap
  // used for global IDs and scalar data, includes ghosts
  dimension = 1;
  oneDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(AbstractDiscretization::getOverlapMap(*comm, decomp, dimension)));

  // threeDimensionalMap
  // used for R3 vector data, e.g., u, v, etc.
  dimension = 3;
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(AbstractDiscretization::getOwnedMap(*comm, decomp, dimension)));

  // threeDimensionalOverlapMap
  // used for R3 vector data, e.g., u, v, etc.,  includes ghosts
  dimension = 3;
  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(AbstractDiscretization::getOverlapMap(*comm, decomp, dimension)));
}

void
PeridigmNS::STKDiscretization::createNeighborhoodData(const QUICKGRID::Data& decomp)
{
   neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
   neighborhoodData->SetNumOwned(decomp.numPoints);
   memcpy(neighborhoodData->OwnedIDs(), 
 		 AbstractDiscretization::getLocalOwnedIds(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.numPoints*sizeof(int));
   memcpy(neighborhoodData->NeighborhoodPtr(), 
 		 decomp.neighborhoodPtr.get(),
 		 decomp.numPoints*sizeof(int));
   neighborhoodData->SetNeighborhoodListSize(decomp.sizeNeighborhoodList);
   memcpy(neighborhoodData->NeighborhoodList(),
 		 AbstractDiscretization::getLocalNeighborList(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.sizeNeighborhoodList*sizeof(int));
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
    int blockID = blockIDs[iID];
	int numUnfilteredNeighbors = unfilteredNeighborhoodList[unfilteredNeighborhoodListIndex++];
    unsigned int numNeighborsIndex = neighborhoodListVec.size();
    neighborhoodListVec.push_back(-1); // placeholder for number of neighbors
    int numNeighbors = 0;
	for(int iNID=0 ; iNID<numUnfilteredNeighbors ; ++iNID){
      int unfilteredNeighborID = unfilteredNeighborhoodList[unfilteredNeighborhoodListIndex++];
      int unfilteredNeighborBlockID = blockIDs[unfilteredNeighborID];
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

Teuchos::RCP< std::vector<double> > PeridigmNS::STKDiscretization::getExodusMeshNodePositions(int globalNodeID)
{
  vector<int>& elementConnectivity = exodusMeshElementConnectivity[globalNodeID];
  unsigned int numNodes = elementConnectivity.size();
  Teuchos::RCP< vector<double> > nodePositionsPtr = Teuchos::RCP< vector<double> >(new vector<double>(3*numNodes));
  vector<double>& nodePositions = *nodePositionsPtr;
  Epetra_Vector& exodusNodePositions = *exodusMeshNodePositions;
  for(unsigned int i=0 ; i<numNodes ; ++i){
    int globalID = elementConnectivity[i];
    int localID = exodusMeshNodePositions->Map().LID(globalID);
    nodePositions[3*i]   = exodusNodePositions[3*localID];
    nodePositions[3*i+1] = exodusNodePositions[3*localID+1];
    nodePositions[3*i+2] = exodusNodePositions[3*localID+2];
  }
  return nodePositionsPtr;
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

double PeridigmNS::STKDiscretization::hexMaxElementDimension(std::vector<double*>& nodeCoordinates) const
{
  double maxDimension = 0.0;
  double dx, dy, dz, diagonal;

  // Check edges

  // Check face diagonals

  // Check element diagonals
  
  // Exodus nodes 1 7
  dx = nodeCoordinates[0][0] - nodeCoordinates[6][0];
  dy = nodeCoordinates[0][1] - nodeCoordinates[6][1];
  dz = nodeCoordinates[0][2] - nodeCoordinates[6][2];
  diagonal = sqrt(dx*dx + dy*dy + dz*dz);
  if(diagonal > maxDimension)
    maxDimension = diagonal;

  // Exodus nodes 2 8
  dx = nodeCoordinates[1][0] - nodeCoordinates[7][0];
  dy = nodeCoordinates[1][1] - nodeCoordinates[7][1];
  dz = nodeCoordinates[1][2] - nodeCoordinates[7][2];
  diagonal = sqrt(dx*dx + dy*dy + dz*dz);
  if(diagonal > maxDimension)
    maxDimension = diagonal;

  // Exodus nodes 3 5
  dx = nodeCoordinates[2][0] - nodeCoordinates[4][0];
  dy = nodeCoordinates[2][1] - nodeCoordinates[4][1];
  dz = nodeCoordinates[2][2] - nodeCoordinates[4][2];
  diagonal = sqrt(dx*dx + dy*dy + dz*dz);
  if(diagonal > maxDimension)
    maxDimension = diagonal;

  // Exodus nodes 4 6
  dx = nodeCoordinates[3][0] - nodeCoordinates[5][0];
  dy = nodeCoordinates[3][1] - nodeCoordinates[5][1];
  dz = nodeCoordinates[3][2] - nodeCoordinates[5][2];
  diagonal = sqrt(dx*dx + dy*dy + dz*dz);
  if(diagonal > maxDimension)
    maxDimension = diagonal;

  return maxDimension;
}
