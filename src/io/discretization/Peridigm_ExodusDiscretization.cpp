/*! \file Peridigm_ExodusDiscretization.cpp */

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

#include "Peridigm_ExodusDiscretization.hpp"
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
#include <sstream>
#include <boost/math/constants/constants.hpp>
#include <boost/algorithm/string.hpp>
#include <exodusII.h>

using namespace std;

//#define DEBUGGING_BACKWARDS_COMPATIBILITY_NEIGHBORHOOD_LIST

PeridigmNS::ExodusDiscretization::ExodusDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params) :
  verbose(false),
  minElementRadius(1.0e50),
  maxElementRadius(0.0),
  storeExodusMesh(false),
  constructInterfaces(false),
  computeIntersections(false),
  maxElementDimension(0.0),
  numBonds(0),
  maxNumBondsPerElem(0),
  myPID(epetra_comm->MyPID()),
  numPID(epetra_comm->NumProc()),
  bondFilterCommand("None"),
  comm(epetra_comm)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Exodus", "Invalid Type in ExodusDiscretization");

  if(params->isParameter("Omit Bonds Between Blocks"))
    bondFilterCommand = params->get<string>("Omit Bonds Between Blocks");
  string meshFileName = params->get<string>("Input Mesh File");

  if(params->isParameter("Verbose")){
    verbose = params->get<bool>("Verbose");
  }

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

PeridigmNS::ExodusDiscretization::~ExodusDiscretization() {
}

void PeridigmNS::ExodusDiscretization::loadData(const string& meshFileName)
{
  // Append processor id information to the file name, if necessary
  string fileName = meshFileName;
  if(numPID != 1){
    stringstream ss;
    ss << "." << numPID << "." << myPID;
    fileName += ss.str();
  }

  // Open the genesis file
  int compWordSize = sizeof(double);
  int ioWordSize = 0;
  float exodusVersion;
  int exodusFileId = ex_open(fileName.c_str(), EX_READ, &compWordSize, &ioWordSize, &exodusVersion);
  if (exodusFileId < 0) reportExodusError(exodusFileId, "ExodusDiscretization::loadData()", "ex_open");

  // Read the initialization parameters
  int numDim, numNodes, numElem, numElemBlocks, numNodeSets, numSideSets;
  char title[MAX_LINE_LENGTH];
  int retval = ex_get_init(exodusFileId, title, &numDim, &numNodes, &numElem, &numElemBlocks, &numNodeSets, &numSideSets);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_init");

  // Node coordinates
  vector<double> exodusNodeCoordX(numNodes), exodusNodeCoordY(numNodes), exodusNodeCoordZ(numNodes);
  retval = ex_get_coord(exodusFileId, &exodusNodeCoordX[0], &exodusNodeCoordY[0], &exodusNodeCoordZ[0]);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_coord");

  // Global node numbering
  vector<int> nodeIdMap(numNodes);
  retval = ex_get_id_map(exodusFileId, EX_NODE_MAP, &nodeIdMap[0]);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_id_map");    
  for(int i=0 ; i<numNodes ; ++i)
    nodeIdMap[i] -= 1; // Note the switch from 1-based indexing to 0-based indexing

  // Global element numbering
  vector<int> elemIdMap(numElem);
  retval = ex_get_id_map(exodusFileId, EX_ELEM_MAP, &elemIdMap[0]);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_id_map");
  for(int i=0 ; i<numElem ; ++i)
    elemIdMap[i] -= 1; // Note the switch from 1-based indexing to 0-based indexing

  // Check for auxiliary node maps and element maps
  int numNodeMaps, numElemMaps;
  retval = ex_get_map_param(exodusFileId, &numNodeMaps, &numElemMaps);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_map_param");

  // DJL
  // This block of code handles the case where an extra elem or node map
  // called "original_global_id_map" is supplied.  This can be the case
  // for parallel decompositions created with decomp or loadbal.
  // If there is an auxiliary map provided that has a different name, throw an
  // error because I don't know what to do with it.
  if(numElemMaps > 0){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numElemMaps > 1,
                                "**** Error in ExodusDiscretization::loadData(), genesis file contains invalid number of auxiliary element maps (>1).\n");
    char mapName[MAX_STR_LENGTH];
    retval = ex_get_name(exodusFileId, EX_ELEM_MAP, 1, mapName);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_name");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(string(mapName) != string("original_global_id_map"),
                                "**** Error in ExodusDiscretization::loadData(), unknown exodus EX_ELEM_MAP: " + string(mapName) + ".\n");
    vector<int> auxMap(numElem);
    retval = ex_get_num_map(exodusFileId, EX_ELEM_MAP, 1, &auxMap[0]);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_num_map");
    for(int i=0 ; i<numElem ; ++i)
      auxMap[i] -= 1; // Note the switch from 1-based indexing to 0-based indexing
    // Use original_global_id_map instead of the map returned by ex_get_id_map()
    elemIdMap = auxMap;
  }
  if(numNodeMaps > 0){
    TEUCHOS_TEST_FOR_EXCEPT_MSG(numNodeMaps > 1,
                                "**** Error in ExodusDiscretization::loadData(), genesis file contains invalid number of auxiliary node maps (>1).\n");
    char mapName[MAX_STR_LENGTH];
    retval = ex_get_name(exodusFileId, EX_NODE_MAP, 1, mapName);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_name");
    TEUCHOS_TEST_FOR_EXCEPT_MSG(string(mapName) != string("original_global_id_map"),
                                "**** Error in ExodusDiscretization::loadData(), unknown exodus EX_NODE_MAP: " + string(mapName) + ".\n");
    vector<int> auxMap(numNodes);
    retval = ex_get_num_map(exodusFileId, EX_NODE_MAP, 1, &auxMap[0]);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_num_map");
    for(int i=0 ; i<numNodes ; ++i)
      auxMap[i] -= 1; // Note the switch from 1-based indexing to 0-based indexing
    // Use original_global_id_map instead of the map returned by ex_get_id_map()
    nodeIdMap = auxMap;
  }

  // Store original exodus-mesh node positions
  // This is for calculation of element-horizon intersetions
  if(storeExodusMesh){
    Epetra_BlockMap exodusMeshNodePositionsMap(-1, numNodes, &nodeIdMap[0], 3, 0, *comm);
    exodusMeshNodePositions = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(exodusMeshNodePositionsMap));
    for(int iNode=0 ; iNode<numNodes ; ++iNode){
      (*exodusMeshNodePositions)[3*iNode]   = exodusNodeCoordX[iNode];
      (*exodusMeshNodePositions)[3*iNode+1] = exodusNodeCoordY[iNode];
      (*exodusMeshNodePositions)[3*iNode+2] = exodusNodeCoordZ[iNode];
    }
  }

  // Create the owned maps
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(-1, numElem, &elemIdMap[0], 1, 0, *comm));
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(-1, numElem, &elemIdMap[0], 3, 0, *comm));

  // Create Epetra_Vectors for the initial positions, volumes, and block_ids
  initialX = Teuchos::rcp(new Epetra_Vector(*threeDimensionalMap));
  cellVolume = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  blockID = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));

  // Process the element blocks
  vector<int> elemBlockIds(numElemBlocks);
  retval = ex_get_elem_blk_ids(exodusFileId, &elemBlockIds[0]);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_blk_ids");

  // Print a warning if the input mesh has side nodes (they will be ignored)
  bool tenNodedTetWarningGiven(false), twentyNodedHexWarningGiven(false);

  int localElemId(0);
  for(int iElemBlock=0 ; iElemBlock<numElemBlocks ; iElemBlock++){

    int elemBlockId = elemBlockIds[iElemBlock];

    // Get the block name, if there is one
    char exodusElemBlockName[MAX_STR_LENGTH];
    retval = ex_get_name(exodusFileId, EX_ELEM_BLOCK, elemBlockId, exodusElemBlockName);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_name");
    // If the block name came back blank, create one that looks like "block_1", "block_2", etc.
    string elemBlockName(exodusElemBlockName);
    if(elemBlockName.size() == 0){
      stringstream ss;
      ss << "block_" << elemBlockId;
      elemBlockName = ss.str();
    }
    TEUCHOS_TEST_FOR_EXCEPT_MSG(elementBlocks->find(elemBlockName) != elementBlocks->end(), "**** Duplicate block found: " + elemBlockName + "\n");
    // Create a list for storing the element ids in this block
    (*elementBlocks)[elemBlockName] = vector<int>();
    vector<int>& elementBlock = (*elementBlocks)[elemBlockName];

    // Get the block parameters and the element connectivity
    char elemType[MAX_STR_LENGTH];
    int numElemThisBlock, numNodesPerElem, numAttributes;
    retval = ex_get_elem_block(exodusFileId, elemBlockId, elemType, &numElemThisBlock, &numNodesPerElem, &numAttributes);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_block");
    ExodusElementType exodusElementType(UNKNOWN_ELEMENT);
    vector<int> conn;
    vector<double> attributes;
    if(numElemThisBlock > 0){
      string elemTypeString(elemType);
      boost::to_upper(elemTypeString);
      if(elemTypeString == string("SPHERE"))
        exodusElementType = SPHERE_ELEMENT;
      else if(elemTypeString == string("TET") || elemTypeString == string("TETRA") || elemTypeString == string("TET4") || elemTypeString == string("TET10"))
        exodusElementType = TET_ELEMENT;
      else if(elemTypeString == string("HEX") || elemTypeString == string("HEX8") || elemTypeString == string("HEX20"))
        exodusElementType = HEX_ELEMENT;
      else{
        string msg = "\n**** Error in loadData(), unknown element type " + elemTypeString + ".\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
      }
      conn.resize(numElemThisBlock*numNodesPerElem);
      retval = ex_get_elem_conn(exodusFileId, elemBlockId, &conn[0]);
      if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_conn");
      if(exodusElementType == SPHERE_ELEMENT){
        attributes.resize(numElemThisBlock*numAttributes);
        retval = ex_get_elem_attr(exodusFileId, elemBlockId, &attributes[0]);
        if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_attr");
      }
    }

    // Loop over the elements in the block
    vector<double> nodeCoordinates(3*numNodesPerElem);
    for(int iElem=0 ; iElem<numElemThisBlock ; iElem++, localElemId++){

      // Node coordinates for all the nodes in this element
      for(int i=0 ; i<numNodesPerElem ; ++i){
        int nodeId = conn[iElem*numNodesPerElem + i] - 1;
        nodeCoordinates[3*i] = exodusNodeCoordX[nodeId];
        nodeCoordinates[3*i+1] = exodusNodeCoordY[nodeId];
        nodeCoordinates[3*i+2] = exodusNodeCoordZ[nodeId];
      }

      int globalElemId = elemIdMap[localElemId];
      elementBlock.push_back(globalElemId);

      // convert elements to spheres and store the initial position and volume
      double volume(0.0);
      vector<double> coord(3);
      if(exodusElementType == SPHERE_ELEMENT){
        // The first attribute is the sph radius (for sph and contact in Sierra/SolidMechanics)
        // The second attribute is the sphere volume
        coord[0] = nodeCoordinates[0];
        coord[1] = nodeCoordinates[1];
        coord[2] = nodeCoordinates[2];
        volume = attributes[iElem*numAttributes + 1];
      }
      else if(exodusElementType == TET_ELEMENT){

        // Warn the user if the tet has ten nodes
        if(numNodesPerElem == 10){
          if(!tenNodedTetWarningGiven){
            cout << "**** Warning on processor " << myPID
                 << ", side nodes being discarded for 10-node tetrahedron element, will be treated as 4-node tetrahedron element." << endl;
            tenNodedTetWarningGiven = true;
          }
        }

        tetCentroidAndVolume(&nodeCoordinates[0], &coord[0], &volume);
      }
      else if(exodusElementType == HEX_ELEMENT){

        // Warn the user if the hex has twenty nodes
        if(numNodesPerElem == 20){
          // 20-noded hex, treat as 8-noded tet
          if(!twentyNodedHexWarningGiven){
            cout << "**** Warning on processor " << myPID
                 << ", side nodes being discarded for 20-node hexahedron element, will be treated as 8-node hexahedron element." << endl;
            twentyNodedHexWarningGiven = true;
          }
        }

        hexCentroidAndVolume(&nodeCoordinates[0], &coord[0], &volume);
      }

      // Store the data in mothership-style vectors
      int epetraLocalId = oneDimensionalMap->LID(globalElemId);
      (*blockID)[epetraLocalId] = elemBlockId;
      (*cellVolume)[epetraLocalId] = volume;
      (*initialX)[3*epetraLocalId] = coord[0];
      (*initialX)[3*epetraLocalId+1] = coord[1];
      (*initialX)[3*epetraLocalId+2] = coord[2];
    }
  }

  // Store element connectivity for original exodus mesh, if needed
  if(storeExodusMesh){
    int* myGlobalElements = oneDimensionalMap->MyGlobalElements();
    vector<int> elementSizeList(numElem);
    char elemType[MAX_STR_LENGTH];
    int numElemThisBlock, numNodesPerElem, numAttributes;
    int index(0);
    for(int iElemBlock=0 ; iElemBlock<numElemBlocks ; ++iElemBlock){
      int elemBlockId = elemBlockIds[iElemBlock];
      retval = ex_get_elem_block(exodusFileId, elemBlockId, elemType, &numElemThisBlock, &numNodesPerElem, &numAttributes);
      if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_block");
      for(int iElem=0 ; iElem<numElemThisBlock ; iElem++, index++){
        elementSizeList[index] = numNodesPerElem;
      }
    }
    int indexBase(0);
    Epetra_BlockMap exodusMeshElementConnectivityMap(-1, numElem, myGlobalElements, &elementSizeList[0], indexBase, *comm);
    exodusMeshElementConnectivity = Teuchos::rcp(new Epetra_Vector(exodusMeshElementConnectivityMap));
    exodusMeshElementConnectivity->PutScalar(-1.0);
    const Epetra_BlockMap& exodusMeshNodePositionsMap = exodusMeshNodePositions->Map();
    index = 0;
    for(int iElemBlock=0 ; iElemBlock<numElemBlocks ; ++iElemBlock){
      int elemBlockId = elemBlockIds[iElemBlock];
      retval = ex_get_elem_block(exodusFileId, elemBlockId, elemType, &numElemThisBlock, &numNodesPerElem, &numAttributes);
      if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_block");
      if(numElemThisBlock > 0){
        vector<int> conn(numElemThisBlock*numNodesPerElem);
        retval = ex_get_elem_conn(exodusFileId, elemBlockId, &conn[0]);
        if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_conn");
        for(int iElem=0 ; iElem<numElemThisBlock ; ++iElem){
          for(int i=0 ; i<numNodesPerElem ; i++, index++){
            int localNodeId = conn[iElem*numNodesPerElem + i] - 1;
            int globalNodeId = exodusMeshNodePositionsMap.GID(localNodeId);
            (*exodusMeshElementConnectivity)[index] = globalNodeId;
          }
        }
      }
    }
  }

  // For each node, record the elements that it belongs to
  vector< vector<int> > elementsThatNodeBelongsTo(numNodes);
  localElemId = 0;
  for(int iElemBlock=0 ; iElemBlock<numElemBlocks ; ++iElemBlock){
    int elemBlockId = elemBlockIds[iElemBlock];
    char elemType[MAX_STR_LENGTH];
    int numElemThisBlock, numNodesPerElem, numAttributes;
    retval = ex_get_elem_block(exodusFileId, elemBlockId, elemType, &numElemThisBlock, &numNodesPerElem, &numAttributes);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_block");
    vector<int> conn(numElemThisBlock*numNodesPerElem);
    if(numElemThisBlock > 0){
      retval = ex_get_elem_conn(exodusFileId, elemBlockId, &conn[0]);
      if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_elem_conn");
      for(int iElem=0 ; iElem<numElemThisBlock ; iElem++, localElemId++){
        for(int i=0 ; i<numNodesPerElem ; i++){
          int localNodeId = conn[iElem*numNodesPerElem + i] - 1;
          elementsThatNodeBelongsTo[localNodeId].push_back(localElemId);
        }
      }
    }
  }

  // Node sets must be converted to new sphere mesh and stored
  nodeSets = Teuchos::rcp< map<string, vector<int> > >(new map<string, vector<int> >() );
  nodeSetIds = Teuchos::rcp< map<string, int> >(new map<string, int>() );
  if(numNodeSets > 0){
    vector<int> exodusNodeSetIds(numNodeSets);
    retval = ex_get_node_set_ids(exodusFileId, &exodusNodeSetIds[0]);
    if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_node_set_ids");
    for(int i=0 ; i<numNodeSets ; ++i){
      int nodeSetId = exodusNodeSetIds[i];
      char exodusNodeSetName[MAX_STR_LENGTH];
      retval = ex_get_name(exodusFileId, EX_NODE_SET, nodeSetId, exodusNodeSetName);
      if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_name");
      // If the node set name came back blank, create one that looks like "nodelist_1", "nodelist_2", etc.
      string nodeSetName(exodusNodeSetName);
      if(nodeSetName.size() == 0){
        stringstream ss;
        ss << "nodelist_" << nodeSetId;
        nodeSetName = ss.str();
      }
      TEUCHOS_TEST_FOR_EXCEPT_MSG(nodeSets->find(nodeSetName) != nodeSets->end(), "**** Duplicate node set found: " + nodeSetName + "\n");
      (*nodeSets)[nodeSetName] = vector<int>();
      (*nodeSetIds)[nodeSetName] = exodusNodeSetIds[i];
    }
    for(map<string, int>::iterator it = nodeSetIds->begin() ; it != nodeSetIds->end() ; it++){
      string nodeSetName = it->first;
      int nodeSetId = it->second;
      int numNodesInSet, numDistributionFactorsInSet;
      retval = ex_get_node_set_param(exodusFileId, nodeSetId, &numNodesInSet, &numDistributionFactorsInSet);
      if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_node_set_param");
      if(numNodesInSet > 0){
        vector<int> nodeSetNodeList(numNodesInSet);
        retval = ex_get_node_set(exodusFileId, nodeSetId, &nodeSetNodeList[0]);
        if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_get_node_set");
        vector<int>& nodeSet = (*nodeSets)[nodeSetName];
        for(int i=0 ; i<numNodesInSet ; ++i){
          // The nodes in the Peridigm sphere mesh that are included in the node set are the nodes
          // at the centers of each exodus element (hex/tet) that contain a node in the exodus node set.
          // This works because, in the Peridigm sphere mesh, the node ids match the element ids (one-to-one relationship).
          const vector<int>& nodes = elementsThatNodeBelongsTo[nodeSetNodeList[i] - 1];  // Note the switch from 1-based indexing to 0-based indexing
          for(unsigned int j=0 ; j<nodes.size() ; ++j){
            int nodeLocalId = nodes[j];
            int nodeGlobalId = threeDimensionalMap->GID(nodeLocalId);
            if( find(nodeSet.begin(), nodeSet.end(), nodeGlobalId) == nodeSet.end() )
              nodeSet.push_back(nodeGlobalId);
          }
        }
      }
    }
  }

  // TODO fix Peridigm to use nodeSetIds instead of string parsing nodeSets, which I think it does in a few places

  if(verbose && myPID == 0){
    stringstream ss;
    ss << "\nGenesis file " << fileName << endl;
    ss << "  title " << title << endl;
    ss << "  number of dimensions " << numDim << endl;
    ss << "  number of nodes " << numNodes << endl;
    ss << "  number of elements " << numElem << endl;
    ss << "  number of blocks " << numElemBlocks << endl;
    ss << "  number of node sets " << numNodeSets << endl;
    ss << "  number of side sets (ignored) " << numSideSets << endl;
    cout << ss.str() << endl;
  }

  // Close the genesis file
  retval = ex_close(exodusFileId);
  if (retval != 0) reportExodusError(retval, "ExodusDiscretization::loadData()", "ex_close");
}

void
PeridigmNS::ExodusDiscretization::constructInterfaceData()
{
  int* const neighPtr = neighborhoodData->NeighborhoodList();
  const int listSize = neighborhoodData->NeighborhoodListSize();

  // faces of a cube as stored by exodus ordering:
  int faces[6][4] = {{4,5,6,7},{7,6,2,3},{0,4,7,3},{0,1,2,3},{0,4,5,1},{1,5,6,2}};

  interfaceData = Teuchos::rcp(new PeridigmNS::InterfaceData);

  TEUCHOS_TEST_FOR_EXCEPTION(storeExodusMesh!=true,logic_error," Exodus mesh should have been stored if this is called.");
  vector<int> leftElements;
  vector<int> rightElements;
  vector<int> numNodes;
  vector<vector<int> > interfaceNodesVec;

  int elemIndex = 0;
  for(int i=0 ; i<listSize; i++){
    const int selfGID = oneDimensionalMap->GID(elemIndex);
    int numNeighbors = neighPtr[i];
    const int numNodesPerElem = exodusMeshElementConnectivity->Map().ElementSize(elemIndex);
    const int myIndex = exodusMeshElementConnectivity->Map().FirstPointInElement(elemIndex);
    // get the connectivity for this element:
    vector<int> selfNodeIds(numNodesPerElem);
    for(int n=0;n<numNodesPerElem;++n){
      selfNodeIds[n] = static_cast<int>( (*exodusMeshElementConnectivity)[myIndex+n] );
    }

    for(int j=0;j<numNeighbors;j++){
      i++;
      int numNodesFound = 0;
      const int GID = oneDimensionalOverlapMap->GID(neighPtr[i]);
      const int neighIndex = exodusMeshElementConnectivity->Map().FirstPointInElement(neighPtr[i]);

      bool shareFace = false;
      // the two elements must have a face in common to be an interface pair:
      vector<int> neighNodeIds(numNodesPerElem);
      vector<int> foundNodeIds(4); // 4 works for tris or quads
      for(int n=0;n<numNodesPerElem;++n)
        neighNodeIds[n] = static_cast<int>( (*exodusMeshElementConnectivity)[neighIndex + n] );

      if(numNodesPerElem==8){ // cube (order needed to prevent folding the interfaces by mixing up the node numbering)
        for(int ni=0;ni<6;++ni){
          numNodesFound = 0;
          for(int nj=0;nj<4;++nj){
            for(int nk=0;nk<numNodesPerElem;++nk){
              if(selfNodeIds[faces[ni][nj]]==neighNodeIds[nk]){
                numNodesFound++;
              }
            }
          }
          if(numNodesFound>2){ // this is a shared face
            for(int kn=0;kn<4;++kn)
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
        for(int ni=0;ni<numNodesPerElem;++ni){
          for(int nj=0;nj<numNodesPerElem;++nj){
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
PeridigmNS::ExodusDiscretization::createNeighborhoodData(int neighborListSize, int* neighborList)
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
PeridigmNS::ExodusDiscretization::filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData)
{
  // Set up a block bonding matrix, which defines whether or not bonds should be formed across blocks
  int numBlocks = getNumBlocks();
  vector< vector<bool> > blockBondingMatrix(numBlocks);
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
PeridigmNS::ExodusDiscretization::getGlobalOwnedMap(int d) const
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
                         endl << "ExodusDiscretization::getGlobalOwnedMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::ExodusDiscretization::getGlobalOverlapMap(int d) const
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
                         endl << "ExodusDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::ExodusDiscretization::getGlobalBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::ExodusDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::ExodusDiscretization::getHorizon() const
{
  return horizonForEachPoint;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::ExodusDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::ExodusDiscretization::getBlockID() const
{
  return blockID;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::ExodusDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::ExodusDiscretization::getNumBonds() const
{
  return numBonds;
}

unsigned int
PeridigmNS::ExodusDiscretization::getMaxNumBondsPerElem() const
{
  return maxNumBondsPerElem;
}

void PeridigmNS::ExodusDiscretization::getExodusMeshNodePositions(int globalNodeID,
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

double PeridigmNS::ExodusDiscretization::computeMaxElementDimension()
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
      length = pow( (3.0*volume)/(4.0*pi), 1.0/3.0 );
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

  localMaxElementDimension = sqrt(localMaxElementDimension);

  // Parallel communication to determine overall maximum element dimension
  double globalMaxElementDimension;
  comm->MaxAll(&localMaxElementDimension, &globalMaxElementDimension, 1);

  return globalMaxElementDimension;
}

void PeridigmNS::ExodusDiscretization::removeNonintersectingNeighborsFromNeighborList(Teuchos::RCP<Epetra_Vector> x,
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

void PeridigmNS::ExodusDiscretization::ghostExodusMeshData()
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
  sort(globalNodeIds.begin(), globalNodeIds.end());

  Epetra_BlockMap overlapExodusMeshNodePositionsMap(numGlobalElements, globalNodeIds.size(), &globalNodeIds[0], 3, indexBase, *comm);
  Teuchos::RCP<Epetra_Vector> overlapExodusMeshNodePositions = Teuchos::rcp(new Epetra_Vector(overlapExodusMeshNodePositionsMap));
  Epetra_Import exodusMeshNodePositionsImporter(overlapExodusMeshNodePositionsMap, exodusMeshNodePositions->Map());
  overlapExodusMeshNodePositions->Import(*exodusMeshNodePositions, exodusMeshNodePositionsImporter, Insert);

  // Set the exodus node positions vector and the connectivity vector to the new overlap vectors
  exodusMeshNodePositions = overlapExodusMeshNodePositions;
  exodusMeshElementConnectivity = overlapExodusMeshElementConnectivity;
}

void PeridigmNS::ExodusDiscretization::reportExodusError(int errorCode, const char *methodName, const char*exodusMethodName)
{
  stringstream ss;
  if (errorCode < 0) { // error
    if (numPID > 1) ss << "Error on PID #" << myPID << ": ";
    ss << "PeridigmNS::OutputManager_ExodusII::" << methodName << "() -- Error code: " << errorCode << " (" << exodusMethodName << ")";
    TEUCHOS_TEST_FOR_EXCEPTION(1, invalid_argument, ss.str());
  }  
  else {
    if (numPID > 1) ss << "Warning on PID #" << myPID << ": ";
    ss << "PeridigmNS::OutputManager_ExodusII::" << methodName << "() -- Warning code: " << errorCode << " (" << exodusMethodName << ")";
    cout << ss.str() << endl;
  }
}

// double PeridigmNS::ExodusDiscretization::hexMaxElementDimension(vector<double*>& nodeCoordinates) const
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

// double PeridigmNS::ExodusDiscretization::tetMaxElementDimension(vector<double*>& nodeCoordinates) const
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
