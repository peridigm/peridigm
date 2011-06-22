/*! \file Peridigm_STKDiscretization.cpp */

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

#include "Peridigm_STKDiscretization.hpp"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdZoltan.h"

#include <Epetra_MpiComm.h>
#include <Teuchos_MPIComm.hpp>
#include <Ionit_Initializer.h>
#include <stk_io/util/UseCase_mesh.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <vector>
#include <list>
#include <sstream>

using namespace std;

PeridigmNS::STKDiscretization::STKDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) :
  comm(epetra_comm),
  numBonds(0),
  myPID(comm->MyPID()),
  numPID(comm->NumProc())
{
  TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "Exodus", "Invalid Type in STKDiscretization");

  string meshFileName = params->get<string>("Input Mesh File");
  double horizon = params->get<double>("Horizon");
  PdGridData decomp = getDecomp(meshFileName, horizon);

  createMaps(decomp);
  createNeighborhoodData(decomp);

  // \todo Move this functionality to base class, it's currently duplicated in PdQuickGridDiscretization.
  // Create the bondMap, a local map used for constitutive data stored on bonds.
  // Due to Epetra_BlockMap restrictions, there can not be any entries with length zero.
  // This means that points with no neighbors can not appear in the bondMap.
  int numMyElementsUpperBound = oneDimensionalMap->NumMyElements();
  int numGlobalElements = -1; 
  int numMyElements = 0;
  int* oneDimensionalMapGlobalElements = oneDimensionalMap->MyGlobalElements();
  int* myGlobalElements = new int[numMyElementsUpperBound];
  int* elementSizeList = new int[numMyElementsUpperBound];
  int* neighborhood = decomp.neighborhood.get();
  int neighborhoodIndex = 0;
  int numPointsWithZeroNeighbors = 0;
  for(int i=0 ; i<decomp.numPoints ; ++i){
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
    neighborhoodIndex += 1 + numNeighbors;
  }
  int indexBase = 0;
  bondMap = Teuchos::rcp(new Epetra_BlockMap(numGlobalElements, numMyElements, myGlobalElements, elementSizeList, indexBase, *comm));
  delete[] myGlobalElements;
  delete[] elementSizeList;

  // 3D only
  TEST_FOR_EXCEPT_MSG(decomp.dimension != 3, "Invalid dimension in decomposition (only 3D is supported)");

  // fill the x vector with the current positions (owned positions only)
  initialX = Teuchos::rcp(new Epetra_Vector(Copy,*threeDimensionalMap,decomp.myX.get()) );

  // fill cell volumes
  cellVolume = Teuchos::rcp(new Epetra_Vector(Copy,*oneDimensionalMap,decomp.cellVolume.get()) );

//   cout << "DECOMP:" << endl;
//   cout << "  dimension " << decomp.dimension << endl;
//   cout << "  globalNumPoints " << decomp.globalNumPoints << endl;
//   cout << "  numPoints " << decomp.numPoints << endl;
//   cout << "  sizeNeighborhoodList " << decomp.sizeNeighborhoodList << endl;
//   cout << "  numExport " << decomp.numExport << endl;
//   cout << "  unPack " << decomp.unPack << endl;
}

PeridigmNS::STKDiscretization::~STKDiscretization() {}


PdGridData PeridigmNS::STKDiscretization::getDecomp(const string& meshFileName,
                                                    double horizon) {

  string meshType = "exodusii";
  string workingDirectory = "";
  Teuchos::RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp_dynamic_cast<const Epetra_MpiComm>(comm, true);
  metaData = Teuchos::rcp(new stk::mesh::fem::FEMMetaData);
  meshData = Teuchos::rcp(new stk::io::util::MeshData);
  Ioss::Init::Initializer io;
  stk::io::util::create_input_mesh(meshType,
                                   meshFileName,
                                   workingDirectory,
                                   mpiComm->Comm(),
                                   *metaData,
                                   *meshData);
  
  int numberOfDimensions = metaData->spatial_dimension();
  TEST_FOR_EXCEPTION(numberOfDimensions != 3, std::invalid_argument, "Peridigm operates only on three-dimensional meshes.");

  // This assigns a null Ioss::GroupingEntity attribute to the universal part
  stk::io::put_io_part_attribute(metaData->universal_part());

  // Loop over the parts and store the element parts, side parts, and node parts.
  const stk::mesh::PartVector& parts = metaData->get_parts();
  stk::mesh::PartVector elementBlocks;
//   stk::mesh::PartVector sideSets;
  stk::mesh::PartVector nodeSets;
  for(stk::mesh::PartVector::const_iterator it = parts.begin(); it != parts.end(); ++it){
    stk::mesh::Part* const part = *it;
    if(part->name()[0] == '{')
      continue;
    if(part->primary_entity_rank() == metaData->element_rank())
      elementBlocks.push_back(part);
//     else if(part->primary_entity_rank() == metaData->side_rank())
//       sideSets.push_back(part);
    else if(part->primary_entity_rank() == metaData->node_rank())
      nodeSets.push_back(part);
    else
      if(myPID == 0)
        cout << "Warning, unknown part type for part " << part->name() << endl;
  }
  if(myPID == 0){
    stringstream ss;
    ss << "Converting input file " << meshFileName << " to sphere mesh:" << endl;
    ss << "  Element blocks:";
    for(stk::mesh::PartVector::const_iterator it = elementBlocks.begin(); it != elementBlocks.end(); ++it)
      ss << " " << (*it)->name();
    ss << endl;
//     ss << "  Side sets:";
//     for(stk::mesh::PartVector::const_iterator it = sideSets.begin(); it != sideSets.end(); ++it)
//       ss << " " << (*it)->name();
//     ss << endl;
    ss << "  Node sets:";
    for(stk::mesh::PartVector::const_iterator it = nodeSets.begin(); it != nodeSets.end(); ++it)
      ss << " " << (*it)->name();
    ss << endl;
    cout << ss.str() << endl;
  }

  if (!metaData->is_FEM_initialized())
    metaData->FEM_initialize(numberOfDimensions);

  stk::mesh::BulkData bulkData(stk::mesh::fem::FEMMetaData::get_meta_data(*metaData), mpiComm->Comm());

  metaData->commit();
  stk::io::util::populate_bulk_data(bulkData, *meshData, "exodusii");
  bulkData.modification_end();

  stk::mesh::Field<double, stk::mesh::Cartesian>* coordinatesField = 
    metaData->get_field< stk::mesh::Field<double, stk::mesh::Cartesian> >("coordinates");

  // Create a selector to select everything in the universal part that is either locally owned or globally shared
  stk::mesh::Selector selector = 
    stk::mesh::Selector( metaData->universal_part() ) & ( stk::mesh::Selector( metaData->locally_owned_part() ) | stk::mesh::Selector( metaData->globally_shared_part() ) );

  // Select the mesh entities that match the selector
  std::vector<stk::mesh::Entity*> elements;
  stk::mesh::get_selected_entities(selector, bulkData.buckets(metaData->element_rank()), elements);

  // Determine the total number of elements in the model
  // \todo There must be a cleaner way to determine the number of elements in a model.
  Teuchos::MPIComm teuchosComm;
  int localElemCount = elements.size();
  int globalElemCount(0);
  teuchosComm.allReduce(&localElemCount, &globalElemCount, 1, Teuchos::MPIComm::INT, Teuchos::MPIComm::SUM);

  // Copy data from stk into a decomp object
  int myNumElements = elements.size();
  int dimension = 3;
  PdGridData decomp = PdQuickGrid::allocatePdGridData(myNumElements, dimension);
  decomp.globalNumPoints = globalElemCount;
  int* globalIds = decomp.myGlobalIDs.get();
  double* volumes = decomp.cellVolume.get();
  double* centroids = decomp.myX.get();

  // loop over the elements and fill the decomp data structure
  for(unsigned int iElem=0 ; iElem<elements.size() ; ++iElem){
    stk::mesh::PairIterRelation nodeRelations = elements[iElem]->node_relations();
    centroids[iElem*3] = 0.0;
    centroids[iElem*3+1] = 0.0;
    centroids[iElem*3+2] = 0.0;
    std::vector<double*> nodeCoordinates(nodeRelations.size());
    int iNode = 0;
    for(stk::mesh::PairIterRelation::iterator it=nodeRelations.begin() ; it!=nodeRelations.end() ; ++it){
      stk::mesh::Entity* node = it->entity();
      double* coordinates = stk::mesh::field_data(*coordinatesField, *node);
      centroids[iElem*3] += coordinates[0];
      centroids[iElem*3+1] += coordinates[1];
      centroids[iElem*3+2] += coordinates[2];
      nodeCoordinates[iNode++] = coordinates;
    }
    centroids[iElem*3] /= nodeRelations.size();
    centroids[iElem*3+1] /= nodeRelations.size();
    centroids[iElem*3+2] /= nodeRelations.size();
    volumes[iElem] = hexVolume(nodeCoordinates);
    globalIds[iElem] = elements[iElem]->identifier() - 1; 
  }

  // loop over the node sets
  for(unsigned int iNodeSet=0 ; iNodeSet<nodeSets.size() ; iNodeSet++){

    const std::string nodeSetName = nodeSets[iNodeSet]->name();
    TEST_FOR_EXCEPT_MSG(sphereMeshNodeSets->find(nodeSetName) != sphereMeshNodeSets->end(), "**** Duplicate node set found: " + nodeSetName + "\n");
    (*sphereMeshNodeSets)[nodeSetName] = std::vector<int>();
    std::vector<int>& sphereMeshNodeSet = (*sphereMeshNodeSets)[nodeSetName];

    // Create a selector for all nodes in the node set
    stk::mesh::Selector selector = 
      stk::mesh::Selector( *nodeSets[iNodeSet] ) & ( stk::mesh::Selector( metaData->locally_owned_part() ) | stk::mesh::Selector( metaData->globally_shared_part() ) );

    // Select the mesh entities that match the selector
    std::vector<stk::mesh::Entity*> nodesInNodeSet;
    stk::mesh::get_selected_entities(selector, bulkData.buckets(metaData->node_rank()), nodesInNodeSet);
    
    // Loop over the nodes in this node set
    std::set<int> elementGlobalIds;
    for(unsigned int iNode=0 ; iNode<nodesInNodeSet.size() ; iNode++){

      // Get all the elements relations for this node and record the element ID in the node set
      // This works because the original element ID is the same as the node ID in the sphere mesh
      stk::mesh::PairIterRelation elementRelations = nodesInNodeSet[iNode]->relations(metaData->element_rank()); 
      for(stk::mesh::PairIterRelation::iterator it=elementRelations.begin() ; it!=elementRelations.end() ; ++it){
        stk::mesh::Entity* element = it->entity();
        int globalId = element->identifier() - 1;
        elementGlobalIds.insert(globalId);
      }
    }

    // Load the node set into the sphereMeshNodeSets container
    for(std::set<int>::iterator it=elementGlobalIds.begin() ; it!=elementGlobalIds.end() ; it++)
      sphereMeshNodeSet.push_back(*it);
  }

  // free the meshData
  meshData = Teuchos::RCP<stk::io::util::MeshData>();

  // call the rebalance function on the current-configuration decomp
  decomp = getLoadBalancedDiscretization(decomp);
  
  // execute neighbor search and update the decomp to include resulting ghosts
  decomp = createAndAddNeighborhood(decomp, horizon);

  return decomp;
}


void
PeridigmNS::STKDiscretization::createMaps(const PdGridData& decomp)
{
  int dimension;

  // oneDimensionalMap
  // used for global IDs and scalar data
  dimension = 1;
  oneDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGrid::getOwnedMap(*comm, decomp, dimension)));

  // oneDimensionalOverlapMap
  // used for global IDs and scalar data, includes ghosts
  dimension = 1;
  oneDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGrid::getOverlapMap(*comm, decomp, dimension)));

  // threeDimensionalMap
  // used for R3 vector data, e.g., u, v, etc.
  dimension = 3;
  threeDimensionalMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGrid::getOwnedMap(*comm, decomp, dimension)));

  // threeDimensionalOverlapMap
  // used for R3 vector data, e.g., u, v, etc.,  includes ghosts
  dimension = 3;
  threeDimensionalOverlapMap = Teuchos::rcp(new Epetra_BlockMap(PdQuickGrid::getOverlapMap(*comm, decomp, dimension)));

}

void
PeridigmNS::STKDiscretization::createNeighborhoodData(const PdGridData& decomp)
{
   neighborhoodData = Teuchos::rcp(new PeridigmNS::NeighborhoodData);
   neighborhoodData->SetNumOwned(decomp.numPoints);
   memcpy(neighborhoodData->OwnedIDs(), 
 		 PdQuickGrid::getLocalOwnedIds(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.numPoints*sizeof(int));
   memcpy(neighborhoodData->NeighborhoodPtr(), 
 		 decomp.neighborhoodPtr.get(),
 		 decomp.numPoints*sizeof(int));
   neighborhoodData->SetNeighborhoodListSize(decomp.sizeNeighborhoodList);
   memcpy(neighborhoodData->NeighborhoodList(),
 		 PdQuickGrid::getLocalNeighborList(decomp, *oneDimensionalOverlapMap).get(),
 		 decomp.sizeNeighborhoodList*sizeof(int));
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::STKDiscretization::getMap(int d) const
{
  switch (d) {
    case 1:
      return oneDimensionalMap;
      break;
    case 3:
      return threeDimensionalMap;
      break;
    default:
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         std::endl << "STKDiscretization::getMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::STKDiscretization::getOverlapMap(int d) const
{
  switch (d) {
    case 1:
      return oneDimensionalOverlapMap;
      break;
    case 3:
      return threeDimensionalOverlapMap;
      break;
    default:
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter, 
                         std::endl << "STKDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::STKDiscretization::getBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::STKDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::STKDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::STKDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::STKDiscretization::getNumBonds() const
{
  return numBonds;
}

double PeridigmNS::STKDiscretization::scalarTripleProduct(std::vector<double>& a,
                                                          std::vector<double>& b,
                                                          std::vector<double>& c) const
{
  double tripleProduct = 
    a[0]*(b[1] - c[2])*(b[2] - c[1]) - a[1]*(b[0] - c[2])*(b[2] - c[0]) + a[2]*(b[0] - c[1])*(b[1] - c[0]);

  return tripleProduct;
}

double PeridigmNS::STKDiscretization::hexVolume(std::vector<double*>& nodeCoordinates) const
{
  int map[] = {0,1,3,2,4,5,7,6};
  std::vector<double> x17(3), x27(3), x47(3), x06(3), x05(3), x03(3), A(3), B(3), C(3);  
  for(int dof=0 ; dof<3 ; ++dof){
    x17[dof] = nodeCoordinates[map[7]][dof] - nodeCoordinates[map[1]][dof];
    x27[dof] = nodeCoordinates[map[7]][dof] - nodeCoordinates[map[2]][dof];
    x47[dof] = nodeCoordinates[map[7]][dof] - nodeCoordinates[map[4]][dof];
    x06[dof] = nodeCoordinates[map[6]][dof] - nodeCoordinates[map[0]][dof];
    x05[dof] = nodeCoordinates[map[5]][dof] - nodeCoordinates[map[0]][dof];
    x03[dof] = nodeCoordinates[map[3]][dof] - nodeCoordinates[map[0]][dof];
  }

  for(int dof=0 ; dof<3 ; ++dof){
    A[dof] = x17[dof] + x06[dof];
    B[dof] = x27[dof] + x05[dof];
    C[dof] = x47[dof] + x03[dof];
  }

  double v1 = fabs( scalarTripleProduct(A,x27,x03) );
  double v2 = fabs( scalarTripleProduct(x06,B,x47) );
  double v3 = fabs( scalarTripleProduct(x17,x05,C) );

  return (v1+v2+v3)/12.0;
}
