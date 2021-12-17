/*! \file Peridigm_Discretization.cpp */

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

#include "Peridigm_Discretization.hpp"
#include "Peridigm_GenesisToPolygons.hpp"
#include <sstream>
#include <set>

using std::set;
using std::string;
using std::stringstream;

Epetra_BlockMap PeridigmNS::Discretization::getOverlap(int ndf, int numShared, int*shared, int numOwned,const  int* owned, const Epetra_Comm& comm){

  int numPoints = numShared+numOwned;
  UTILITIES::Array<int> ids(numPoints);
  int *ptr = ids.get();

  for(int j=0;j<numOwned;j++,ptr++)
    *ptr=owned[j];

  for(int j=0;j<numShared;j++,ptr++)
    *ptr=shared[j];

  return Epetra_BlockMap(-1,numPoints, ids.get(),ndf, 0,comm);
}

UTILITIES::Array<int> PeridigmNS::Discretization::getSharedGlobalIds(const QUICKGRID::Data& gridData){
  set<int> ownedIds(gridData.myGlobalIDs.get(),gridData.myGlobalIDs.get()+gridData.numPoints);
  set<int> shared;
  int *neighPtr = gridData.neighborhoodPtr.get();
  int *neigh = gridData.neighborhood.get();
  set<int>::const_iterator ownedIdsEnd = ownedIds.end();
  for(size_t p=0;p<gridData.numPoints;p++){
    int ptr = neighPtr[p];
    int numNeigh = neigh[ptr];
    for(int n=1;n<=numNeigh;n++){
      int id = neigh[ptr+n];
      /*
       * look for id in owned points
       */
      if(ownedIdsEnd == ownedIds.find(id)){
        /*
         * add this point to shared
         */
        shared.insert(id);
      }
    }
  }

  // Copy set into shared ptr
  UTILITIES::Array<int> sharedGlobalIds(shared.size());
  int *sharedPtr = sharedGlobalIds.get();
  set<int>::iterator it;
  for ( it=shared.begin() ; it != shared.end(); it++, sharedPtr++ )
    *sharedPtr = *it;

  return sharedGlobalIds;
}

std::shared_ptr<int> PeridigmNS::Discretization::getLocalOwnedIds(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap){
  UTILITIES::Array<int> localIds(gridData.numPoints);
  int *lIds = localIds.get();
  int *end = localIds.get()+gridData.numPoints;
  int *gIds = gridData.myGlobalIDs.get();
  for(; lIds != end;lIds++, gIds++)
    *lIds = overlapMap.LID(*gIds);
  return localIds.get_shared_ptr();
}

std::shared_ptr<int> PeridigmNS::Discretization::getLocalNeighborList(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap){
  UTILITIES::Array<int> localNeighborList(gridData.sizeNeighborhoodList);
  int *localNeig = localNeighborList.get();
  int *neighPtr = gridData.neighborhoodPtr.get();
  int *neigh = gridData.neighborhood.get();
  for(size_t p=0;p<gridData.numPoints;p++){
    int ptr = neighPtr[p];
    int numNeigh = neigh[ptr];
    localNeig[ptr]=numNeigh;
    for(int n=1;n<=numNeigh;n++){
      int gid = neigh[ptr+n];
      int localId = overlapMap.LID(gid);
      localNeig[ptr+n] = localId;
    }
  }
  return localNeighborList.get_shared_ptr();
}

Epetra_BlockMap PeridigmNS::Discretization::getOwnedMap(const Epetra_Comm& comm,const QUICKGRID::Data& gridData, int ndf) {
  int numShared=0;
  int *sharedPtr=NULL;
  int numOwned = gridData.numPoints;
  const int *ownedPtr = gridData.myGlobalIDs.get();
  return getOverlap(ndf, numShared,sharedPtr,numOwned,ownedPtr,comm);
}

Epetra_BlockMap PeridigmNS::Discretization::getOverlapMap(const Epetra_Comm& comm,const QUICKGRID::Data& gridData, int ndf) {
  UTILITIES::Array<int> sharedGIDS = getSharedGlobalIds(gridData);
  std::shared_ptr<int> sharedPtr = sharedGIDS.get_shared_ptr();
  int numShared = sharedGIDS.get_size();
  int *shared = sharedPtr.get();
  int *owned = gridData.myGlobalIDs.get();
  int numOwned = gridData.numPoints;
  return getOverlap(ndf,numShared,shared,numOwned,owned,comm);
}

void PeridigmNS::Discretization::createBondFilters(const Teuchos::RCP<Teuchos::ParameterList>& params){
  if(params->isSublist("Bond Filters")){
    Teuchos::RCP<Teuchos::ParameterList> bondFilterParameters = sublist(params, "Bond Filters");
    for (Teuchos::ParameterList::ConstIterator it = bondFilterParameters->begin(); it != bondFilterParameters->end(); ++it) {
      string parameterListName = it->first;
      Teuchos::ParameterList params = bondFilterParameters->sublist(parameterListName);
      string type = params.get<string>("Type");
      if(type == "Rectangular_Plane"){
        double normal[3], lowerLeftCorner[3], bottomUnitVector[3], bottomLength, sideLength;
        normal[0] = params.get<double>("Normal_X");
        normal[1] = params.get<double>("Normal_Y");
        normal[2] = params.get<double>("Normal_Z");
        lowerLeftCorner[0] = params.get<double>("Lower_Left_Corner_X");
        lowerLeftCorner[1] = params.get<double>("Lower_Left_Corner_Y");
        lowerLeftCorner[2] = params.get<double>("Lower_Left_Corner_Z");
        bottomUnitVector[0] = params.get<double>("Bottom_Unit_Vector_X");
        bottomUnitVector[1] = params.get<double>("Bottom_Unit_Vector_Y");
        bottomUnitVector[2] = params.get<double>("Bottom_Unit_Vector_Z");
        bottomLength = params.get<double>("Bottom_Length");
        sideLength = params.get<double>("Side_Length");
        PdBondFilter::FinitePlane finitePlane(normal, lowerLeftCorner, bottomUnitVector, bottomLength, sideLength);
        std::shared_ptr<PdBondFilter::BondFilter> bondFilter(new PdBondFilter::FinitePlaneFilter(finitePlane));
        bondFilters.push_back(bondFilter);
      }
      else if(type == "Disk"){
        double center[3], normal[3], radius;
        center[0] = params.get<double>("Center_X");
        center[1] = params.get<double>("Center_Y");
        center[2] = params.get<double>("Center_Z");
        normal[0] = params.get<double>("Normal_X");
        normal[1] = params.get<double>("Normal_Y");
        normal[2] = params.get<double>("Normal_Z");
        radius = params.get<double>("Radius");
        std::shared_ptr<PdBondFilter::BondFilter> bondFilter(new PdBondFilter::DiskFilter(center, normal, radius));
        bondFilters.push_back(bondFilter);
      }
      else if(type == "Exodus Mesh"){
        std::string file_name = params.get<string>("File Name");
        std::vector< std::vector< std::vector<double> > > polygons;
        GenesisToPolygons(file_name, polygons);
        std::shared_ptr<PdBondFilter::BondFilter> bondFilter(nullptr);
        for(unsigned int i=0 ; i<polygons.size() ; i++){
          if(polygons[i].size() == 3)
            bondFilter = std::make_shared<PdBondFilter::TriangleFilter>(PdBondFilter::TriangleFilter(polygons[i][0].data(), polygons[i][1].data(), polygons[i][2].data()));
          else if(polygons[i].size() == 4)
            bondFilter = std::make_shared<PdBondFilter::QuadFilter>(PdBondFilter::QuadFilter(polygons[i][0].data(), polygons[i][1].data(), polygons[i][2].data(), polygons[i][3].data()));
          bondFilters.push_back(bondFilter);
        }
      }
      else{
        string msg  = "\n**** Error, invalid bond filter type:  " + type;
               msg += "\n**** Allowable types are:  Rectangular_Plane, Disk, and Exodus Mesh.\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
      }
    }
  }
}

int PeridigmNS::Discretization::blockNameToBlockId(string blockName) const {

  int bID = -1;

  size_t underscore_loc = blockName.find_last_of('_');

  // Dream3D naming convention is GRAIN1_SET
  size_t prefix_loc = blockName.find_first_of("GRAIN");
  size_t prefix_len = 5;
  size_t suffix_loc = blockName.find_last_of("_SET");
  size_t suffix_len = 4;

  if (prefix_loc != string::npos && suffix_loc != string::npos) {
    // Dream3D naming convention
    size_t block_num_len = suffix_loc - suffix_len - prefix_loc - prefix_len + 1;
    stringstream blockIDSS(blockName.substr(prefix_loc+prefix_len, block_num_len));
    blockIDSS >> bID;
  }
  else if (underscore_loc != string::npos) {
    // Cubit naming convention
    stringstream blockIDSS(blockName.substr(underscore_loc+1, blockName.size()));
    blockIDSS >> bID;
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(underscore_loc == string::npos, "\n**** Parse error, invalid block name: " + blockName + "\n");
  }

  return bID;
}

void PeridigmNS::Discretization::createBondMapAndCheckForZeroNeighbors(Teuchos::RCP<Epetra_BlockMap>& bondMap,
                                                                       const Teuchos::RCP<Epetra_BlockMap> oneDimensionalMap,
                                                                       const Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData,
                                                                       unsigned int & numBonds,
                                                                       unsigned int & maxNumBondsPerElem) const {
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
      /*
       * Print warning message when no neighbor exist for this point.
       * This can cause in problems downstream e.g. in computation of the
       * jacobian matrix in quasi-static or implicit solver.
       */
      int globalID = oneDimensionalMapGlobalElements[i];
      int localID  = oneDimensionalMap->LID(globalID);
      Teuchos::RCP<Epetra_Vector> xPtr = getInitialX();
      std::cout << "\n *** WARNING: PeridigmNS::Discretization::createBondMapAndCheckForZeroNeighbors(..)\n";
      std::cout << "\tNumber of neighbors is less or equal 1!\n\tThis is probably a problem for the evaluation of a jacobian.\n";
      std::cout << "\tGlobal point id = " << globalID << "\n"
                << "\tx,y,z = " << (*xPtr)[3*localID] << ", " << (*xPtr)[3*localID+1] << ", " << (*xPtr)[3*localID+2] << "\n" << std::endl;
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
}