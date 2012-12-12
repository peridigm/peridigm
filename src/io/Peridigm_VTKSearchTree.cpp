/*! \file Peridigm_VTKSearchTree.cpp */
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

#include "Peridigm_VTKSearchTree.hpp"
#include <vtkDoubleArray.h>
#include <vtkIdList.h>

using namespace std;

PeridigmNS::VTKSearchTree::VTKSearchTree(int numPoints, double* coordinates) : SearchTree(numPoints, coordinates)
{
  vtkSmartPointer<vtkUnstructuredGrid> overlapGrid = getGrid(coordinates, numPoints);
  vtkKdTreePointLocator* kdTree = vtkKdTreePointLocator::New();
  kdTree->SetDataSet(overlapGrid);
}

PeridigmNS::VTKSearchTree::~VTKSearchTree()
{
  kdTree->Delete();
}

void PeridigmNS::VTKSearchTree::FindPointsWithinRadius(const double* point, double searchRadius, vector<int>& neighborList)
{
  vtkIdList* kdTreeList = vtkIdList::New();

  // Note that list returned includes this point.
  kdTree->FindPointsWithinRadius(searchRadius, point, kdTreeList);

  neighborList.clear();
  for(int i=0; i<kdTreeList->GetNumberOfIds(); ++i){
    neighborList.push_back( kdTreeList->GetId(i) );
  }
}

vtkSmartPointer<vtkCellArray> PeridigmNS::VTKSearchTree::getCellArray(vtkIdType numCells){

	// This is a really simply link array
	// It assumes that the cell id # is same as point id #
	// By construction, the points generated here satisfy this condition

	vtkSmartPointer<vtkIdTypeArray>  links = vtkSmartPointer<vtkIdTypeArray>::New();
	vtkIdType onePointPerCell(1);
	links->SetNumberOfValues(2*numCells);
	for(int cell=0; cell<numCells; cell++){
      links->SetValue(2*cell, onePointPerCell);
      links->SetValue(2*cell+1, cell);
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->SetCells(numCells, links);

	return cells;
}

vtkSmartPointer<vtkUnstructuredGrid> PeridigmNS::VTKSearchTree::getGrid(double *points, int numPoints){

  // NOTE: pointer to data is used directly here;
  // This method lets the user specify data to be held by the array.  The
  // array argument is a pointer to the data.  Set save to 1 to keep the class
  // from deleting the array when it cleans up or reallocates memory.
  // The class uses the actual array provided; it does not copy the data
  // from the suppled array.
  // NOTE: Also -- the supplied coordinates 'y' should be a vector with 3 components
  // per point

  vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkDoubleArray> ptsData = vtkSmartPointer<vtkDoubleArray>::New();
  int numComponents = 3;
  ptsData->SetNumberOfComponents(numComponents);
  int save = 1;
  ptsData->SetArray(points, numPoints*numComponents, save);
  pts->SetData(ptsData);

  // Creates a grid of "points" and "cells"
  vtkSmartPointer<vtkCellArray> cells = getCellArray(numPoints);
  vtkSmartPointer<vtkUnstructuredGrid> grid = getGrid(pts, cells);

  return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> PeridigmNS::VTKSearchTree::getGrid(const vtkSmartPointer<vtkPoints>& points, const vtkSmartPointer<vtkCellArray>& cells, VTKCellType type){

  // Cells constructed are "vertex" type cells
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  grid->SetPoints(points);
  grid->SetCells(type, cells);

  return grid;
}
