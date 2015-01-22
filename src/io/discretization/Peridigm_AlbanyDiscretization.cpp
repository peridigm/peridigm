/*! \file Peridigm_AlbanyDiscretization.cpp */

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

#include "Peridigm_AlbanyDiscretization.hpp"
#include "Peridigm_ProximitySearch.hpp"
#include "Peridigm_HorizonManager.hpp"

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

PeridigmNS::AlbanyDiscretization::AlbanyDiscretization(const MPI_Comm& mpiComm,
                                                       const Teuchos::RCP<Teuchos::ParameterList>& params,
						       int numGlobalIds,
						       const int* globalIds,
						       const double* refCoord,
						       const double* volume,
						       const int* blockId) :
  minElementRadius(1.0e50),
  maxElementRadius(0.0),
  maxElementDimension(0.0),
  numBonds(0),
  maxNumBondsPerElem(0),
  bondFilterCommand("None")
{
  comm = Teuchos::rcp(new Epetra_MpiComm(mpiComm));

  Teuchos::RCP<Teuchos::ParameterList> discretizationParams = Teuchos::rcpFromRef(params->sublist("Discretization", true));
  Teuchos::RCP<Teuchos::ParameterList> blockParams = Teuchos::rcpFromRef(params->sublist("Blocks", true));

  TEUCHOS_TEST_FOR_EXCEPT_MSG(discretizationParams->get<string>("Type") != "Albany", "Invalid Type in AlbanyDiscretization");

  // Create the owned maps
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(-1, numGlobalIds, globalIds, 1, 0, *comm));
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(-1, numGlobalIds, globalIds, 3, 0, *comm));

  // Create Epetra_Vectors and fill them with the provided data
  cellVolume = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*oneDimensionalMap));
  blockID = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*oneDimensionalMap));
  for(int i=0 ; i<numGlobalIds ; ++i){
    (*cellVolume)[i] = volume[i];
    (*blockID)[i] = blockId[i];
  }
  initialX = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*threeDimensionalMap));
  for(int i=0 ; i<3*numGlobalIds ; ++i){
    (*initialX)[i] = refCoord[i];
  }

  // Set up bond filters
  if(discretizationParams->isParameter("Omit Bonds Between Blocks"))
    bondFilterCommand = discretizationParams->get<string>("Omit Bonds Between Blocks");
  createBondFilters(discretizationParams);

  // Create a list of on-processor elements for each block
  createBlockElementLists();

  // Assign the correct horizon to each node
  PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();
  horizonManager.loadHorizonInformationFromBlockParameters(*blockParams);
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

  // Perform the proximity search to identify neighbors
  int neighborListSize;
  int* neighborList;
  ProximitySearch::GlobalProximitySearch(initialX, horizonForEachPoint, oneDimensionalOverlapMap, neighborListSize, neighborList, bondFilters);

  createNeighborhoodData(neighborListSize, neighborList);

  // Create the three-dimensional overlap map based on the one-dimensional overlap map
  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(-1, 
                                                                oneDimensionalOverlapMap->NumMyElements(),
                                                                oneDimensionalOverlapMap->MyGlobalElements(),
                                                                3,
                                                                0,
                                                                oneDimensionalOverlapMap->Comm()));

  // \todo Move this functionality to base class, it's currently duplicated in multiple discretizations.
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
  comm->MinAll(&localMin[0], &globalMin[0], 1);
  minElementRadius = globalMin[0];
  localMin[0] = maxElementRadius;
  comm->MaxAll(&localMin[0], &globalMin[0], 1);
  maxElementRadius = globalMin[0];
}

PeridigmNS::AlbanyDiscretization::~AlbanyDiscretization() {}

void PeridigmNS::AlbanyDiscretization::createBlockElementLists() {

  // Each processor needs the compete list of block ids.
  // Parallel communication is needed because if a processor has zero nodes/elements for a given block,
  // it will not know that block exists.

  // Find the unique block ids on processor
  set<int> uniqueBlockIds;
  for(int i=0 ; i<blockID->MyLength() ; ++i)
    uniqueBlockIds.insert( static_cast<int>( (*blockID)[i] ) );

  // Broadcast the unique block ids so that all processors are aware of the full block list
  Teuchos::RCP<const Teuchos::Comm<int> > teuchosComm = Teuchos::createMpiComm<int>(Teuchos::opaqueWrapper<MPI_Comm>(MPI_COMM_WORLD));
  int numLocalUniqueBlockIds = static_cast<int>( uniqueBlockIds.size() );
  int maxNumberOfUniqueBlockIds;
  reduceAll(*teuchosComm, Teuchos::REDUCE_MAX, 1, &numLocalUniqueBlockIds, &maxNumberOfUniqueBlockIds);
  int arraySize = maxNumberOfUniqueBlockIds * comm->NumProc();
  vector<int> uniqueBlockIdsLocal(arraySize, -1);
  int index = maxNumberOfUniqueBlockIds * comm->MyPID();
  for(set<int>::const_iterator it = uniqueBlockIds.begin() ; it != uniqueBlockIds.end() ; it++)
    uniqueBlockIdsLocal[index++] = *it;
  vector<int> uniqueBlockIdsGlobal(arraySize);  
  reduceAll(*teuchosComm, Teuchos::REDUCE_MAX, arraySize, &uniqueBlockIdsLocal[0], &uniqueBlockIdsGlobal[0]);

  // Insert off-processor block ids into the set of unique block ids
  for(int i=0 ; i<arraySize ; ++i){
    if(uniqueBlockIdsGlobal[i] != -1)
      uniqueBlockIds.insert(uniqueBlockIdsGlobal[i]);
  }

  // Initialize the element list for each block
  // Force blocks with no on-processor elements to have an entry in the elementBlocks map
  for(set<int>::iterator it=uniqueBlockIds.begin() ; it!=uniqueBlockIds.end() ; it++){
    stringstream blockName;
    blockName << "block_" << *it;
    (*elementBlocks)[blockName.str()] = std::vector<int>();
  }

  // Create the element list for each block
  for(int i=0 ; i<blockID->MyLength() ; ++i){
    stringstream blockName;
    blockName << "block_" << (*blockID)[i];
    TEUCHOS_TEST_FOR_EXCEPT_MSG(elementBlocks->find(blockName.str()) == elementBlocks->end(),
                                "\n**** Error in AlbanyDiscretization::loadData(), invalid block id.\n");
    int globalID = blockID->Map().GID(i);
    (*elementBlocks)[blockName.str()].push_back(globalID);
  }
}

void
PeridigmNS::AlbanyDiscretization::createNeighborhoodData(int neighborListSize, int* neighborList)
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
PeridigmNS::AlbanyDiscretization::filterBonds(Teuchos::RCP<PeridigmNS::NeighborhoodData> unfilteredNeighborhoodData)
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
PeridigmNS::AlbanyDiscretization::getGlobalOwnedMap(int d) const
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
                         std::endl << "AlbanyDiscretization::getGlobalOwnedMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::AlbanyDiscretization::getGlobalOverlapMap(int d) const
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
                         std::endl << "AlbanyDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::AlbanyDiscretization::getGlobalBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::AlbanyDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::AlbanyDiscretization::getHorizon() const
{
  return horizonForEachPoint;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::AlbanyDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::AlbanyDiscretization::getBlockID() const
{
  return blockID;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::AlbanyDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::AlbanyDiscretization::getNumBonds() const
{
  return numBonds;
}

unsigned int
PeridigmNS::AlbanyDiscretization::getMaxNumBondsPerElem() const
{
  return maxNumBondsPerElem;
}

