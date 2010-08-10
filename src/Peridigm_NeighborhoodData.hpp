/*! \file Peridigm_NeighborhoodData.hpp */

// ***********************************************************************
//
//                             Peridigm
//                 Copyright (2009) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// 
// Questions? 
// David J. Littlewood   djlittl@sandia.gov 
// John A. Mitchell      jamitch@sandia.gov
// Michael L. Parks      mlparks@sandia.gov
// Stewart A. Silling    sasilli@sandia.gov
//
// ***********************************************************************

#ifndef PERIDIGM_NEIGHBORHOODDATA_HPP
#define PERIDIGM_NEIGHBORHOODDATA_HPP

namespace Peridigm {

class NeighborhoodData {

public:

  NeighborhoodData() 
    : numOwnedPoints(0), ownedIDs(0), neighborhoodListSize(0), neighborhoodList(0), neighborhoodPtr(0) {}

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
