/*! \file Peridigm_Compute_Neighborhood_Volume.cpp */

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

#include <vector>

#include "Peridigm_Compute_Number_Of_Neighbors.hpp"
#include "Peridigm_Field.hpp"

PeridigmNS::Compute_Number_Of_Neighbors::Compute_Number_Of_Neighbors(Teuchos::RCP<const Epetra_Comm> epetraComm_)
  : Compute(epetraComm_), m_partialVolumeFieldId(-1), m_numberOfNeighborsFieldId(-1)
{
  FieldManager& fieldManager = FieldManager::self();
  m_numberOfNeighborsFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Number_Of_Neighbors");
  m_fieldIds.push_back(m_numberOfNeighborsFieldId);
  if(fieldManager.hasField("Partial_Volume")){
    m_partialVolumeFieldId = fieldManager.getFieldId("Partial_Volume");
    m_fieldIds.push_back(m_partialVolumeFieldId);
  }
}

PeridigmNS::Compute_Number_Of_Neighbors::~Compute_Number_Of_Neighbors(){}

void PeridigmNS::Compute_Number_Of_Neighbors::initialize( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {

  std::vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    double *numberOfNeighbors;
    blockIt->getData(m_numberOfNeighborsFieldId, PeridigmField::STEP_NONE)->ExtractView(&numberOfNeighbors);

    double *partialVolume = 0;
    if( m_partialVolumeFieldId != -1 )
      blockIt->getData(m_partialVolumeFieldId, PeridigmField::STEP_NONE)->ExtractView(&partialVolume);

    int neighborhoodListIndex = 0;
    int bondIndex = 0;
    for(int iID=0 ; iID<numOwnedPoints ; ++iID){
      int numNeighbors = neighborhoodList[neighborhoodListIndex++];

      if(partialVolume == 0){
      
        // The analysis does not utilize partial volumes,
        // all neighbors are counted.

        numberOfNeighbors[iID] = numNeighbors;
      }
      else {

        // Partial volumes are present.

        // Allow for the case in which an element is within the neighborhood list,
        // but actually does not have any volume within the neighborhood.
        // This may or may not be possible, depending on the details of the neighborhood
        // search, which must extend beyond the neighborhood when partial volumes are
        // considered.

        numberOfNeighbors[iID] = 0;
        for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
          double fraction = partialVolume[bondIndex++];
          if(fraction > 0.0)
            numberOfNeighbors[iID] += 1;
        }
      }
      neighborhoodListIndex += numNeighbors;
    }
  }
}

int PeridigmNS::Compute_Number_Of_Neighbors::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const {
  return 0;
}
