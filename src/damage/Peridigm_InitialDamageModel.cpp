/*! \file Peridigm_InitialDamageModel.cpp */

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

#include "Peridigm_InitialDamageModel.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_GenesisToTriangles.hpp"
#include "BondFilter.h"

PeridigmNS::InitialDamageModel::InitialDamageModel(const Teuchos::ParameterList& params)
  : DamageModel(params),  m_modelCoordinatesFieldId(-1), m_damageFieldId(-1), m_bondDamageFieldId(-1)
{
  Teuchos::ParameterList temp = params;
  m_bondFilterParameters = temp.sublist("Bond Filters", true);

  PeridigmNS::FieldManager& fieldManager = PeridigmNS::FieldManager::self();
  m_modelCoordinatesFieldId = fieldManager.getFieldId("Model_Coordinates");
  m_damageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::ELEMENT, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Damage");
  m_bondDamageFieldId = fieldManager.getFieldId(PeridigmNS::PeridigmField::BOND, PeridigmNS::PeridigmField::SCALAR, PeridigmNS::PeridigmField::TWO_STEP, "Bond_Damage");

  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_bondDamageFieldId);
}

PeridigmNS::InitialDamageModel::~InitialDamageModel()
{}

void
PeridigmNS::InitialDamageModel::initialize(const double dt,
                                           const int numOwnedPoints,
                                           const int* ownedIDs,
                                           const int* neighborhoodList,
                                           PeridigmNS::DataManager& dataManager) const
{
  // Construct the bond filters
  std::vector< std::shared_ptr<PdBondFilter::BondFilter> > bondFilters;
  for (Teuchos::ParameterList::ConstIterator it = m_bondFilterParameters.begin(); it != m_bondFilterParameters.end(); ++it) {
    std::string parameterListName = it->first;
    Teuchos::ParameterList params = m_bondFilterParameters.sublist(parameterListName);
    std::string type = params.get<std::string>("Type");
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
    else if(type == "Exodus Mesh"){
      std::string file_name = params.get<std::string>("File Name");
      std::vector< std::vector< std::vector<double> > > triangles;
      GenesisToTriangles(file_name, triangles);
      for(unsigned int i=0 ; i<triangles.size() ; i++){
        std::shared_ptr<PdBondFilter::BondFilter> bondFilter(new PdBondFilter::TriangleFilter(triangles[i][0].data(), triangles[i][1].data(), triangles[i][2].data()));
        bondFilters.push_back(bondFilter);
      }
    }
    else{
      std::string msg = "\n**** Error, invalid bond filter type:  " + type;
      msg += "\n**** Allowable types are:  Rectangular_Plane, Exodus Mesh\n";
      TEUCHOS_TEST_FOR_EXCEPT_MSG(true, msg);
    }
  }

  double *x, *damageN, *damageNP1, *bondDamageN, *bondDamageNP1;
  dataManager.getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE)->ExtractView(&x);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_N)->ExtractView(&damageN);
  dataManager.getData(m_damageFieldId, PeridigmField::STEP_NP1)->ExtractView(&damageNP1);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_N)->ExtractView(&bondDamageN);
  dataManager.getData(m_bondDamageFieldId, PeridigmField::STEP_NP1)->ExtractView(&bondDamageNP1);

  // Apply the bond filters
  int neighborhoodListIndex = 0;
  int bondIndex = 0;
  std::vector<int> treeList(1,0);
  std::vector<double> xOverlap(3);
  int ptLocalID = -1;
  bool bondFlag;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    int nodeID = ownedIDs[iID];
    double* pt = &x[3*nodeID];
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      int neighborID = neighborhoodList[neighborhoodListIndex++];
      xOverlap[0] = x[3*neighborID];
      xOverlap[1] = x[3*neighborID+1];
      xOverlap[2] = x[3*neighborID+2];
      bondDamageNP1[bondIndex] = 0.0;
      for (std::vector< std::shared_ptr<PdBondFilter::BondFilter> >::iterator it = bondFilters.begin(); it != bondFilters.end() ; it++) {
        bondFlag = false;
        (*it)->filterBonds(treeList, pt, ptLocalID, xOverlap.data(), &bondFlag);
        if (bondFlag)
          bondDamageNP1[bondIndex] = 1.0;
      }
      bondIndex += 1;
    }
  }

  //  Update the element damage (percent of bonds broken)
  neighborhoodListIndex = 0;
  bondIndex = 0;
  for(int iID=0 ; iID<numOwnedPoints ; ++iID){
    int nodeID = ownedIDs[iID];
    int numNeighbors = neighborhoodList[neighborhoodListIndex++];
    neighborhoodListIndex += numNeighbors;
    double totalDamage = 0.0;
    for(int iNID=0 ; iNID<numNeighbors ; ++iNID){
      bondDamageN[bondIndex] = bondDamageNP1[bondIndex]; // set bond damage for both step N and NP1
      totalDamage += bondDamageNP1[bondIndex++];
    }
    if(numNeighbors > 0)
      totalDamage /= numNeighbors;
    else
      totalDamage = 0.0;
    damageN[nodeID] = totalDamage;
    damageNP1[nodeID] = totalDamage;
  }
}

void
PeridigmNS::InitialDamageModel::computeDamage(const double dt,
                                              const int numOwnedPoints,
                                              const int* ownedIDs,
                                              const int* neighborhoodList,
                                              PeridigmNS::DataManager& dataManager) const
{}
