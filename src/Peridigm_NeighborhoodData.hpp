/*! \file Peridigm_NeighborhoodData.hpp */

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

#ifndef PERIDIGM_NEIGHBORHOODDATA_HPP
#define PERIDIGM_NEIGHBORHOODDATA_HPP

namespace PeridigmNS {

class NeighborhoodData {

public:

  NeighborhoodData() 
    : numOwnedPoints(0), ownedIDs(0), neighborhoodListSize(0), neighborhoodList(0), neighborhoodPtr(0) {}

  NeighborhoodData(const NeighborhoodData& other)
    : numOwnedPoints(0), ownedIDs(0), neighborhoodListSize(0), neighborhoodList(0), neighborhoodPtr(0)
  {
    SetNumOwned(other.NumOwnedPoints());
    SetNeighborhoodListSize(other.NeighborhoodListSize());
    memcpy(ownedIDs, other.ownedIDs, numOwnedPoints*sizeof(int));
    memcpy(neighborhoodPtr, other.neighborhoodPtr, numOwnedPoints*sizeof(int));
    memcpy(neighborhoodList, other.neighborhoodList, neighborhoodListSize*sizeof(int));
  }

  ~NeighborhoodData(){
	if(ownedIDs != 0)
	  delete[] ownedIDs;
	if(neighborhoodList != 0)
	  delete[] neighborhoodList;
    if(neighborhoodPtr != 0)
      delete[] neighborhoodPtr;
  }

  void SetNumOwned(int numOwned){
	numOwnedPoints = numOwned;
	if(ownedIDs != 0)
	  delete[] ownedIDs;
	ownedIDs = new int[numOwned];
    if(neighborhoodPtr != 0)
      delete[] neighborhoodPtr;
    neighborhoodPtr = new int[numOwned];
  }

  void SetNeighborhoodListSize(int neighborhoodSize){
	neighborhoodListSize = neighborhoodSize;
	if(neighborhoodList != 0)
	  delete[] neighborhoodList;
	neighborhoodList = new int[neighborhoodListSize];
  }

  const int NumOwnedPoints() const{
	return numOwnedPoints;
  }

  int* const OwnedIDs() const{
	return ownedIDs;
  }

  int* const NeighborhoodPtr() const{
	return neighborhoodPtr;
  }

  const int NeighborhoodListSize() const{
	return neighborhoodListSize;
  }

  int* const NeighborhoodList() const{
	return neighborhoodList;
  }

protected:
  int numOwnedPoints;
  int* ownedIDs;
  int neighborhoodListSize;
  int* neighborhoodList;
  int* neighborhoodPtr;
};

}

#endif // PERIDIGM_NEIGHBORHOODDATA_HPP
