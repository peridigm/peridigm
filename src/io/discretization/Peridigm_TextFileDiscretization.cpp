/*! \file Peridigm_TextFileDiscretization.cpp */

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

#include "Peridigm_TextFileDiscretization.hpp"
#include "Peridigm_HorizonManager.hpp"
#include "pdneigh/NeighborhoodList.h"
#include "pdneigh/PdZoltan.h"

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_RCP.hpp>

#include <sstream>
#include <fstream>

#include <boost/algorithm/string/trim.hpp>

using namespace std;

PeridigmNS::TextFileDiscretization::TextFileDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                           const Teuchos::RCP<Teuchos::ParameterList>& params) :
  minElementRadius(1.0e50),
  maxElementRadius(0.0),
  maxElementDimension(0.0),
  numBonds(0),
  maxNumBondsPerElem(0),
  myPID(epetra_comm->MyPID()),
  numPID(epetra_comm->NumProc()),
  bondFilterCommand("None"),
  comm(epetra_comm)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Text File", "Invalid Type in TextFileDiscretization");

  string meshFileName = params->get<string>("Input Mesh File");
  if(params->isParameter("Omit Bonds Between Blocks"))
    bondFilterCommand = params->get<string>("Omit Bonds Between Blocks");

  // Set up bond filters
  createBondFilters(params);

  QUICKGRID::Data decomp = getDecomp(meshFileName, params);

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

  // 3D only
  TEUCHOS_TEST_FOR_EXCEPT_MSG(decomp.dimension != 3, "Invalid dimension in decomposition (only 3D is supported)");

  // fill the x vector with the current positions (owned positions only)
  initialX = Teuchos::rcp(new Epetra_Vector(Copy, *threeDimensionalMap, decomp.myX.get()));

  // fill cell volumes
  cellVolume = Teuchos::rcp(new Epetra_Vector(Copy,*oneDimensionalMap,decomp.cellVolume.get()) );

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

PeridigmNS::TextFileDiscretization::~TextFileDiscretization() {}


QUICKGRID::Data PeridigmNS::TextFileDiscretization::getDecomp(const string& textFileName,
                                                              const Teuchos::RCP<Teuchos::ParameterList>& params) {

  // Read data from the text file
  vector<double> coordinates;
  vector<double> volumes;
  vector<int> blockIds;

  // Read the text file on the root processor
  if(myPID == 0){
    ifstream inFile(textFileName.c_str());
    TEUCHOS_TEST_FOR_EXCEPT_MSG(!inFile.is_open(), "**** Error opening discretization text file.\n");
    while(inFile.good()){
      string str;
      getline(inFile, str);
      boost::trim(str);
      // Ignore comment lines, otherwise parse
      if( !(str[0] == '#' || str[0] == '/' || str[0] == '*' || str.size() == 0) ){
        istringstream iss(str);
        vector<double> data;
        copy(istream_iterator<double>(iss),
             istream_iterator<double>(),
             back_inserter<vector<double> >(data));
        // Check for obvious problems with the data
        if(data.size() != 5){
          string msg = "\n**** Error parsing text file, invalid line: " + str + "\n";
          TEUCHOS_TEST_FOR_EXCEPT_MSG(data.size() != 5, msg);
        }
        // Store the coordinates, block id, and volumes
        coordinates.push_back(data[0]);
        coordinates.push_back(data[1]);
        coordinates.push_back(data[2]);
        blockIds.push_back(static_cast<int>(data[3]));
        volumes.push_back(data[4]);
      }
    }
    inFile.close();
  }

  int numElements = static_cast<int>(blockIds.size());
  TEUCHOS_TEST_FOR_EXCEPT_MSG(myPID == 0 && numElements < 1, "**** Error reading discretization text file, no data found.\n");

  // Record the block ids on the root processor
  set<int> uniqueBlockIds;
  if(myPID == 0){
    for(unsigned int i=0 ; i<blockIds.size() ; ++i)
      uniqueBlockIds.insert(blockIds[i]);
  }

  // Broadcast necessary data from root processor
  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  int numGlobalElements;
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, 1, &numElements, &numGlobalElements);

  // Broadcast the unique block ids so that all processors are aware of the full block list
  // This is necessary because if a processor does not have any elements for a given block, it will be unaware the
  // given block exists, which causes problems downstream
  int numLocalUniqueBlockIds = static_cast<int>( uniqueBlockIds.size() );
  int numGlobalUniqueBlockIds;
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, 1, &numLocalUniqueBlockIds, &numGlobalUniqueBlockIds);
  vector<int> uniqueLocalBlockIds(numGlobalUniqueBlockIds, 0);
  int index = 0;
  for(set<int>::const_iterator it = uniqueBlockIds.begin() ; it != uniqueBlockIds.end() ; it++)
    uniqueLocalBlockIds[index++] = *it;
  vector<int> uniqueGlobalBlockIds(numGlobalUniqueBlockIds);  
  reduceAll(*teuchosComm, Teuchos::REDUCE_SUM, numGlobalUniqueBlockIds, &uniqueLocalBlockIds[0], &uniqueGlobalBlockIds[0]);

  // Create list of global ids
  vector<int> globalIds(numElements);
  for(unsigned int i=0 ; i<globalIds.size() ; ++i)
    globalIds[i] = i;

  // Copy data into a decomp object
  int dimension = 3;
  QUICKGRID::Data decomp = QUICKGRID::allocatePdGridData(numElements, dimension);
  decomp.globalNumPoints = numGlobalElements;
  memcpy(decomp.myGlobalIDs.get(), &globalIds[0], numElements*sizeof(int)); 
  memcpy(decomp.cellVolume.get(), &volumes[0], numElements*sizeof(double)); 
  memcpy(decomp.myX.get(), &coordinates[0], 3*numElements*sizeof(double));

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
  for(unsigned int i=0 ; i<blockIds.size() ; ++i)
    tempBlockIDPtr[i] = blockIds[i];

  // call the rebalance function on the current-configuration decomp
  decomp = PDNEIGH::getLoadBalancedDiscretization(decomp);

  // create a (throw-away) one-dimensional owned map in the rebalanced configuration
  Epetra_BlockMap rebalancedMap(decomp.globalNumPoints, decomp.numPoints, decomp.myGlobalIDs.get(), 1, 0, *comm);

  // Create a (throw-away) blockID vector corresponding to the load balanced decomposition
  Epetra_Vector rebalancedBlockID(rebalancedMap);
  Epetra_Import rebalancedImporter(rebalancedBlockID.Map(), tempBlockID.Map());
  rebalancedBlockID.Import(tempBlockID, rebalancedImporter, Insert);

  // Initialize the element list for each block
  // Force blocks with no on-processor elements to have an entry in the elementBlocks map
  for(unsigned int i=0 ; i<uniqueGlobalBlockIds.size() ; i++){
    stringstream blockName;
    blockName << "block_" << uniqueGlobalBlockIds[i];
    (*elementBlocks)[blockName.str()] = std::vector<int>();
  }

  // Create the element list for each block
  for(int i=0 ; i<rebalancedBlockID.MyLength() ; ++i){
    stringstream blockName;
    blockName << "block_" << rebalancedBlockID[i];
    TEUCHOS_TEST_FOR_EXCEPT_MSG(elementBlocks->find(blockName.str()) == elementBlocks->end(),
                                "\n**** Error in TextFileDiscretization::getDecomp(), invalid block id.\n");
    int globalID = rebalancedBlockID.Map().GID(i);
    (*elementBlocks)[blockName.str()].push_back(globalID);
  }

  // Record the horizon for each point
  PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();
  Teuchos::RCP<Epetra_Vector> rebalancedHorizonForEachPoint = Teuchos::rcp(new Epetra_Vector(rebalancedMap));
  double* rebalancedX = decomp.myX.get();
  for(map<string, vector<int> >::const_iterator it = elementBlocks->begin() ; it != elementBlocks->end() ; it++){
    const string& blockName = it->first;
    const vector<int>& globalIds = it->second;

    bool hasConstantHorizon = horizonManager.blockHasConstantHorizon(blockName);
    double constantHorizonValue(0.0);
    if(hasConstantHorizon)
      constantHorizonValue = horizonManager.getBlockConstantHorizonValue(blockName);

    for(unsigned int i=0 ; i<globalIds.size() ; ++i){
      int localId = rebalancedMap.LID(globalIds[i]);
      if(hasConstantHorizon){
        (*rebalancedHorizonForEachPoint)[localId] = constantHorizonValue;
      }
      else{
        double x = rebalancedX[localId*3];
        double y = rebalancedX[localId*3 + 1];
        double z = rebalancedX[localId*3 + 2];
        double horizon = horizonManager.evaluateHorizon(blockName, x, y, z);
        (*rebalancedHorizonForEachPoint)[localId] = horizon;
      }
    }
  }

  // execute neighbor search and update the decomp to include resulting ghosts
  std::tr1::shared_ptr<const Epetra_Comm> commSp(comm.getRawPtr(), NonDeleter<const Epetra_Comm>());
  Teuchos::RCP<PDNEIGH::NeighborhoodList> list;
  TEUCHOS_TEST_FOR_EXCEPT_MSG(bondFilters.size() > 1, "\n****Error:  Multiple bond filters currently unsupported.\n");

  if(bondFilters.size() == 0)
    list = Teuchos::rcp(new PDNEIGH::NeighborhoodList(commSp,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,rebalancedHorizonForEachPoint));
  else if(bondFilters.size() == 1)
    list = Teuchos::rcp(new PDNEIGH::NeighborhoodList(commSp,decomp.zoltanPtr.get(),decomp.numPoints,decomp.myGlobalIDs,decomp.myX,rebalancedHorizonForEachPoint,bondFilters[0]));
  decomp.neighborhood=list->get_neighborhood();
  decomp.sizeNeighborhoodList=list->get_size_neighborhood_list();
  decomp.neighborhoodPtr=list->get_neighborhood_ptr();

  // Create all the maps.
  createMaps(decomp);

  // Create the blockID vector corresponding to the load balanced decomposition
  blockID = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  Epetra_Import tempImporter(blockID->Map(), tempBlockID.Map());
  blockID->Import(tempBlockID, tempImporter, Insert);

  // Create the horizonForEachPonit vector corresponding to the load balanced decomposition
  horizonForEachPoint = Teuchos::rcp(new Epetra_Vector(*oneDimensionalMap));
  Epetra_Import horizonImporter(horizonForEachPoint->Map(), rebalancedHorizonForEachPoint->Map());
  horizonForEachPoint->Import(*rebalancedHorizonForEachPoint, horizonImporter, Insert);

  return decomp;
}

void
PeridigmNS::TextFileDiscretization::createMaps(const QUICKGRID::Data& decomp)
{
  int dimension;

  // oneDimensionalMap
  // used for global IDs and scalar data
  dimension = 1;
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOwnedMap(*comm, decomp, dimension)));

  // oneDimensionalOverlapMap
  // used for global IDs and scalar data, includes ghosts
  dimension = 1;
  oneDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOverlapMap(*comm, decomp, dimension)));

  // threeDimensionalMap
  // used for R3 vector data, e.g., u, v, etc.
  dimension = 3;
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOwnedMap(*comm, decomp, dimension)));

  // threeDimensionalOverlapMap
  // used for R3 vector data, e.g., u, v, etc.,  includes ghosts
  dimension = 3;
  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(Discretization::getOverlapMap(*comm, decomp, dimension)));
}

void
PeridigmNS::TextFileDiscretization::createNeighborhoodData(const QUICKGRID::Data& decomp)
{
   neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
   neighborhoodData->SetNumOwned(decomp.numPoints);
   memcpy(neighborhoodData->OwnedIDs(), 
 		 Discretization::getLocalOwnedIds(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.numPoints*sizeof(int));
   memcpy(neighborhoodData->NeighborhoodPtr(), 
 		 decomp.neighborhoodPtr.get(),
 		 decomp.numPoints*sizeof(int));
   neighborhoodData->SetNeighborhoodListSize(decomp.sizeNeighborhoodList);
   memcpy(neighborhoodData->NeighborhoodList(),
 		 Discretization::getLocalNeighborList(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.sizeNeighborhoodList*sizeof(int));
   neighborhoodData = filterBonds(neighborhoodData);
}

Teuchos::RCP<PeridigmNS::NeighborhoodData>
PeridigmNS::TextFileDiscretization::filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData)
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
PeridigmNS::TextFileDiscretization::getGlobalOwnedMap(int d) const
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
                         std::endl << "TextFileDiscretization::getGlobalOwnedMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::TextFileDiscretization::getGlobalOverlapMap(int d) const
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
                         std::endl << "TextFileDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::TextFileDiscretization::getGlobalBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getHorizon() const
{
  return horizonForEachPoint;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::TextFileDiscretization::getBlockID() const
{
  return blockID;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::TextFileDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::TextFileDiscretization::getNumBonds() const
{
  return numBonds;
}

unsigned int
PeridigmNS::TextFileDiscretization::getMaxNumBondsPerElem() const
{
  return maxNumBondsPerElem;
}

