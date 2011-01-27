/*! \file Peridigm_PdQuickGridDiscretization.cpp */

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

#include "Peridigm_PdQuickGridDiscretization.hpp"
#include "PdQuickGrid.h"
#include "PdQuickGridParallel.h"
#include "PdZoltan.h"
#include <vector>
#include <list>
#include <sstream>

using namespace std;

PeridigmNS::PdQuickGridDiscretization::PdQuickGridDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                                 const Teuchos::RCP<Teuchos::ParameterList>& params) :
  comm(epetra_comm),
  numBonds(0),
  myPID(comm->MyPID()),
  numPID(comm->NumProc())
{
  TEST_FOR_EXCEPT_MSG(params->get<string>("Type") != "PdQuickGrid", "Invalid Type in PdQuickGridDiscretization");
  PdGridData decomp = getDiscretization(params);

  createMaps(decomp);
  createNeighborhoodData(decomp);

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
}


PeridigmNS::PdQuickGridDiscretization::PdQuickGridDiscretization(const Teuchos::RCP<const Epetra_Comm>& epetra_comm,
                                                                 const Teuchos::RCP<const PdGridData>& decomp) :
  comm(epetra_comm),
  numBonds(0),
  myPID(comm->MyPID()),
  numPID(comm->NumProc())
{
  createMaps(*decomp);
  createNeighborhoodData(*decomp);

  // Create the bondMap, a local map used for constitutive data stored on bonds.
  // Due to Epetra_BlockMap restrictions, there can not be any entries with length zero.
  // This means that points with no neighbors can not appear in the bondMap.
  int numMyElementsUpperBound = oneDimensionalMap->NumMyElements();
  int numGlobalElements = -1; 
  int numMyElements = 0;
  int* oneDimensionalMapGlobalElements = oneDimensionalMap->MyGlobalElements();
  int* myGlobalElements = new int[numMyElementsUpperBound];
  int* elementSizeList = new int[numMyElementsUpperBound];
  int* neighborhood = decomp->neighborhood.get();
  int neighborhoodIndex = 0;
  int numPointsWithZeroNeighbors = 0;
  for(int i=0 ; i<decomp->numPoints ; ++i){
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
  TEST_FOR_EXCEPT_MSG(decomp->dimension != 3, "Invalid dimension in decomposition (only 3D is supported)");

  // fill the x vector with the current positions (owned positions only)
  initialX = Teuchos::rcp(new Epetra_Vector(Copy, *threeDimensionalMap, decomp->myX.get()) );

  // fill cell volumes
  cellVolume = Teuchos::rcp(new Epetra_Vector(Copy, *oneDimensionalMap, decomp->cellVolume.get()) );
}

PeridigmNS::PdQuickGridDiscretization::~PdQuickGridDiscretization() {}

PdGridData PeridigmNS::PdQuickGridDiscretization::getDiscretization(const Teuchos::RCP<Teuchos::ParameterList>& params) {

  // This is the type of norm used to create neighborhood lists
  PdQuickGrid::NormFunctionPointer neighborhoodType = PdQuickGrid::NoOpNorm;
  Teuchos::ParameterEntry* normTypeEntry=params->getEntryPtr("NeighborhoodType");
  if(NULL!=normTypeEntry){
    std::string normType = params->get<string>("NeighborhoodType");
    if(normType=="Spherical") neighborhoodType = PdQuickGrid::SphericalNorm;
  }

  // Get the horizion
  horizon = params->get<double>("Horizon");

  // param list should have a "sublist" with different types that we switch on here
  PdGridData decomp;
  if (params->isSublist("TensorProduct3DMeshGenerator")){
    Teuchos::RCP<Teuchos::ParameterList> pdQuickGridParamList = Teuchos::rcp(&(params->sublist("TensorProduct3DMeshGenerator")), false);
    double xStart = pdQuickGridParamList->get<double>("X Origin");
    double yStart = pdQuickGridParamList->get<double>("Y Origin");
    double zStart = pdQuickGridParamList->get<double>("Z Origin");
    double xLength = pdQuickGridParamList->get<double>("X Length");
    double yLength = pdQuickGridParamList->get<double>("Y Length");
    double zLength = pdQuickGridParamList->get<double>("Z Length");
    int nx = pdQuickGridParamList->get<int>("Number Points X");
    int ny = pdQuickGridParamList->get<int>("Number Points Y");
    int nz = pdQuickGridParamList->get<int>("Number Points Z");

    const PdQuickGrid::PdQPointSet1d xSpec(nx,xStart,xLength);
    const PdQuickGrid::PdQPointSet1d ySpec(ny,yStart,yLength);
    const PdQuickGrid::PdQPointSet1d zSpec(nz,zStart,zLength);

    // Create abstract decomposition iterator
    PdQuickGrid::TensorProduct3DMeshGenerator cellPerProcIter(numPID,horizon,xSpec,ySpec,zSpec,neighborhoodType);
    decomp =  PdQuickGrid::getDiscretization(myPID, cellPerProcIter);
    // Load balance and write new decomposition
    #ifdef HAVE_MPI
      decomp = getLoadBalancedDiscretization(decomp);
    #endif
  } 
  else if (params->isSublist("TensorProductCylinderMeshGenerator")){
    Teuchos::RCP<Teuchos::ParameterList> pdQuickGridParamList = Teuchos::rcp(&(params->sublist("TensorProductCylinderMeshGenerator")), false);
    double innerRadius    = pdQuickGridParamList->get<double>("Inner Radius");
    double outerRadius    = pdQuickGridParamList->get<double>("Outer Radius");
    double cylinderLength = pdQuickGridParamList->get<double>("Cylinder Length");
    int numRings          = pdQuickGridParamList->get<int>("Number Points Radius");
    double xC             = pdQuickGridParamList->get<double>("Ring Center x");
    double yC             = pdQuickGridParamList->get<double>("Ring Center y");
    double zStart         = pdQuickGridParamList->get<double>("Z Origin");

    // Create 2d Ring
    std::valarray<double> center(0.0,3);
    center[0] = xC;
    center[1] = yC;
    center[2] = 0;

    // Note that zStart is used for the 1D spec along cylinder axis
    PdQuickGrid::PdQRing2d ring2dSpec(center,innerRadius,outerRadius,numRings);

    // Create 1d Spec along cylinder axis
    // Compute number of cells along length of cylinder so that aspect ratio
    // is cells is approximately 1.
    // Cell sizes along axis are not exactly "cellSize" since last cell
    // would be a fraction of a cellSize -- so 1 is added to numCellsAlongAxis.
    // Actual cell sizes are slightly smaller than "cellSize" because of this.
    double cellSize = ring2dSpec.getRaySpec().getCellSize();
    int numCellsAxis = (int)(cylinderLength/cellSize)+1;
    PdQuickGrid::PdQPointSet1d axisSpec(numCellsAxis,zStart,cylinderLength);

    // Create abstract decomposition iterator
    PdQuickGrid::TensorProductCylinderMeshGenerator cellPerProcIter(numPID, horizon,ring2dSpec, axisSpec,neighborhoodType);
    decomp =  PdQuickGrid::getDiscretization(myPID, cellPerProcIter);
    // Load balance and write new decomposition
    #ifdef HAVE_MPI
      decomp = getLoadBalancedDiscretization(decomp);
    #endif
  } 
  else { // ERROR
    TEST_FOR_EXCEPT_MSG(true, "Invalid Type in PdQuickGridDiscretization");
  }

  return decomp;
}

void
PeridigmNS::PdQuickGridDiscretization::createMaps(const PdGridData& decomp)
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
PeridigmNS::PdQuickGridDiscretization::createNeighborhoodData(const PdGridData& decomp)
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
PeridigmNS::PdQuickGridDiscretization::getMap(int d) const
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
                         std::endl << "PdQuickGridDiscretization::getMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::PdQuickGridDiscretization::getOverlapMap(int d) const
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
                         std::endl << "PdQuickGridDiscretization::getOverlapMap(int d) only supports dimensions d=1 or d=3. Supplied dimension d=" << d << std::endl); 
    }
}

Teuchos::RCP<const Epetra_BlockMap>
PeridigmNS::PdQuickGridDiscretization::getBondMap() const
{
  return bondMap;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::PdQuickGridDiscretization::getInitialX() const
{
  return initialX;
}

Teuchos::RCP<Epetra_Vector>
PeridigmNS::PdQuickGridDiscretization::getCellVolume() const
{
  return cellVolume;
}

Teuchos::RCP<PeridigmNS::NeighborhoodData> 
PeridigmNS::PdQuickGridDiscretization::getNeighborhoodData() const
{
  return neighborhoodData;
}

unsigned int
PeridigmNS::PdQuickGridDiscretization::getNumBonds() const
{
  return numBonds;
}

