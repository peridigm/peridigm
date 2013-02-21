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

using std::tr1::shared_ptr;

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
	std::set<int> ownedIds(gridData.myGlobalIDs.get(),gridData.myGlobalIDs.get()+gridData.numPoints);
	std::set<int> shared;
	int *neighPtr = gridData.neighborhoodPtr.get();
	int *neigh = gridData.neighborhood.get();
	std::set<int>::const_iterator ownedIdsEnd = ownedIds.end();
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
    std::set<int>::iterator it;
	for ( it=shared.begin() ; it != shared.end(); it++, sharedPtr++ )
		*sharedPtr = *it;

	return sharedGlobalIds;
}

shared_ptr<int> PeridigmNS::Discretization::getLocalOwnedIds(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap){
	UTILITIES::Array<int> localIds(gridData.numPoints);
	int *lIds = localIds.get();
	int *end = localIds.get()+gridData.numPoints;
	int *gIds = gridData.myGlobalIDs.get();
	for(; lIds != end;lIds++, gIds++)
		*lIds = overlapMap.LID(*gIds);
	return localIds.get_shared_ptr();
}

shared_ptr<int> PeridigmNS::Discretization::getLocalNeighborList(const QUICKGRID::Data& gridData, const Epetra_BlockMap& overlapMap){
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
	shared_ptr<int> sharedPtr = sharedGIDS.get_shared_ptr();
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
        std::tr1::shared_ptr<PdBondFilter::BondFilter> bondFilter(new PdBondFilter::FinitePlaneFilter(finitePlane));
        bondFilters.push_back(bondFilter);        
      }
      else{
        string msg = "\n**** Error, invalid bond filter type:  " + type;
        msg += "\n**** Allowable types are:  Rectangular_Plane\n";
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
      }
    }
  }
}



