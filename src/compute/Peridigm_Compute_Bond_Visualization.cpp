/*! \file Peridigm_Compute_Bond_Visualization.cpp */

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

#include "Peridigm_Compute_Bond_Visualization.hpp"
#include <sstream>
#include <fstream>

using namespace std;

PeridigmNS::Compute_Bond_Visualization::Compute_Bond_Visualization(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                                   Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                                   Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_), m_modelCoordinatesFieldId(-1), m_bondVisualizationFieldId(-1), m_specifiedNodesOnly(false)
{
  // Name of the output file
  string baseName = "BondVisualization";
  if(!params.is_null() && params->isParameter("File Name")){
    baseName = params->get<string>("File Name");
  }
  stringstream ss;
  ss << baseName << "_Proc" << epetraComm_->MyPID() << ".vtk";
  m_fileName = ss.str();

  // By default, display all the bonds in the model
  // Alternatively, show bonds for only a specified node
  if(!params.is_null() && params->isParameter("Node Id")){
    m_specifiedNodesOnly = true;
    int nodeId = params->get<int>("Node Id") - 1; // Switch from 1-based indexing (exodus) to 0-based indexing (epetra)
    m_specifiedNodes.push_back(nodeId);
  }

  FieldManager& fieldManager = FieldManager::self();
  m_modelCoordinatesFieldId  = fieldManager.getFieldId("Model_Coordinates");
  // This is a dummy field, it's just here so that the user can trigger the creation of this compute class
  // by requesting "Bond_Visualization" in the Output section of the input deck
  m_bondVisualizationFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, "Bond_Visualization");
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_bondVisualizationFieldId);
}

PeridigmNS::Compute_Bond_Visualization::~Compute_Bond_Visualization(){}

void PeridigmNS::Compute_Bond_Visualization::initialize( Teuchos::RCP< vector<PeridigmNS::Block> > blocks ) {

  std::string filename;
  std::set<int> globalNodeIdSet;
  std::vector<int> globalNodeIds;
  std::map< int, std::vector<double> > coordinates;
  std::vector< std::pair<int, int> > connectivity;
  int globalId(-1);
  bool storeData;

  // Record the global id and model coordinates for all the nodes on this processor
  vector<PeridigmNS::Block>::iterator blockIt;
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    Teuchos::RCP<Epetra_Vector> modelCoordinates = blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE);
    const Epetra_BlockMap& map = modelCoordinates->Map();

    vector<double> coord(3);
    for(int iNode=0 ; iNode<map.NumMyElements() ; ++iNode){
      globalId = map.GID(iNode);
      coord[0] = (*modelCoordinates)[3*iNode];
      coord[1] = (*modelCoordinates)[3*iNode+1];
      coord[2] = (*modelCoordinates)[3*iNode+2];
      coordinates[globalId] = coord;
    }
  }

  // Record the bonds as bar elements (pair of global ids)
  for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++){
    Teuchos::RCP<PeridigmNS::NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
    const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
    const int* neighborhoodList = neighborhoodData->NeighborhoodList();
    Teuchos::RCP<const Epetra_BlockMap> map = blockIt->getOverlapVectorPointMap();

    int iID, iNID, localNeighborId, globalId, globalNeighborId, numNeighbors(0), neighborhoodListIndex(0);
    for(iID=0 ; iID<numOwnedPoints ; ++iID){
      globalId = map->GID(iID);

      storeData = true;
      if(m_specifiedNodesOnly){
	storeData = false;
	if( find(m_specifiedNodes.begin(), m_specifiedNodes.end(), globalId) != m_specifiedNodes.end() )
	  storeData = true;
      }

      numNeighbors = neighborhoodList[neighborhoodListIndex++];
      for(iNID=0 ; iNID<numNeighbors ; ++iNID){
        localNeighborId = neighborhoodList[neighborhoodListIndex++];
        globalNeighborId = map->GID(localNeighborId);
	if(storeData){
	  connectivity.push_back( std::pair<int, int>(globalId, globalNeighborId) );
	  globalNodeIdSet.insert(globalId);
	  globalNodeIdSet.insert(globalNeighborId);
	}
      }
    }
  }

  for(std::set<int>::iterator it = globalNodeIdSet.begin() ; it != globalNodeIdSet.end() ; it++)
    globalNodeIds.push_back(*it);
  sort(globalNodeIds.begin(), globalNodeIds.end());

  writeVTK(m_fileName, globalNodeIds, coordinates, connectivity);
}

int PeridigmNS::Compute_Bond_Visualization::compute( Teuchos::RCP< vector<PeridigmNS::Block> > blocks ) const {
  return 0;
}

void PeridigmNS::Compute_Bond_Visualization::writeVTK(std::string fileName,
                                                      std::vector<int>& globalNodeIds,
                                                      std::map< int, std::vector<double> >& coordinates,
                                                      std::vector< std::pair<int, int> >& connectivity)
{
  int globalId;
  ofstream visFile;
  visFile.open(fileName.c_str());

  // file version and identifier
  visFile << "# vtk DataFile Version 3.0" << endl;

  // header
  visFile << "This file was generated by the Peridigm Bond Visualization Compute Class" << endl;

  // file format ASCII | BINARY
  visFile << "ASCII" << endl;

  // dataset structure STRUCTURED_POINTS | STRUCTURED_GRID | UNSTRUCTURED_GRID | POLYDATA | RECTILINEAR_GRID | FIELD
  visFile << "DATASET UNSTRUCTURED_GRID" << endl;

  std::map<int, int> globalIdToVtkId;

  // data
  visFile << "POINTS " << globalNodeIds.size() << " float" << endl;
  for(unsigned int i=0 ; i<globalNodeIds.size() ; ++i){
    globalId = globalNodeIds[i];
    const vector<double>& coord = coordinates.at(globalId);
    visFile << coord[0] << " " << coord[1] << " " << coord[2] << " " << endl;
    globalIdToVtkId[globalId] = i;
  }

  int id_1, id_2;

  visFile << "CELLS " << connectivity.size() << " " << 3*connectivity.size() << endl;
  for(vector< pair<int, int> >::iterator it = connectivity.begin() ; it != connectivity.end() ; it++){
    id_1 = globalIdToVtkId.at(it->first);
    id_2 = globalIdToVtkId.at(it->second);
    visFile << "2 " << id_1 << " " << id_2 << endl;
  }

  // note that cell type 3 is VTK_LINE
  visFile << "CELL_TYPES " << connectivity.size() << endl;
  for(unsigned int i=0 ; i<connectivity.size() ; ++i)
    visFile << 3 << endl;

  visFile.close();
}
