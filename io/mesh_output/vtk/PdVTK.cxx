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

#include "PdVTK.h"
#include "Vector.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVertex.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include <vtkDoubleArray.h>

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <sstream>


namespace PdVTK {

/*
 * These are private
 */
vtkSmartPointer<vtkUnstructuredGrid>  getGrid(const vtkSmartPointer<vtkPoints>& x);
vtkSmartPointer<vtkCellArray> getCellArray(vtkIdType numCells);

CollectionWriter::CollectionWriter(const char* _fileName, int numProcs, int rank, VTK_FILE_TYPE type)
: fileName(_fileName), times(), writer(getWriter(_fileName,numProcs,rank,type)){}

/**
 * This call writes a single pvtu file associated with the grid at this time step
 * @param t -- time step added to collection
 * @param grid -- unstructured grid associated with time step
 */
void CollectionWriter::writeTimeStep(double t, vtkSmartPointer<vtkUnstructuredGrid> grid){
	times.push_back(t);
	std::size_t index = times.size() - 1;
	string name(getPVTU_fileName(index,fileName));
	writer->SetFileName(name.c_str());
	write(writer,grid);
}

string CollectionWriter::getPVTU_fileName(int index, const char* _fileName) const {
	std::stringstream ss;
	ss << _fileName << "_t" << index << ".pvtu";
	return ss.str();
}

/**
 **This function writes the collection file that allows paraview to open and read all time steps
 * @param comment -- Comment added to xml header; should include "\n"
 */
void CollectionWriter::close(const string& comment) {
	string outFileName(fileName);
	outFileName += ".pvd";
	std::fstream fStream(outFileName.c_str(), fstream::out);

	fStream << "<?xml version=\"1.0\"?>\n"
			<< "<!--\n"
			<< comment
			<<	"-->\n"
		    << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
		    << "\t<Collection>\n";

	std::size_t numPVTU = times.size();

	for(std::size_t f=0;f<numPVTU;f++){
		double t = times[f];
		fStream << "\t\t<DataSet timestep=\"" << t << "\" file=\"" << getPVTU_fileName(f,fileName) << "\"/>\n";
	}

	fStream << "\t</Collection>\n"
		    << "</VTKFile>\n";

	fStream.close();
}

vtkSmartPointer<vtkXMLPUnstructuredGridWriter > getWriter(const char* _fileName, int numProcs, int rank, VTK_FILE_TYPE type){
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter > w = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();
	w->SetNumberOfPieces(numProcs);
	w->SetStartPiece(rank);
	w->SetEndPiece(rank);
	w->SetFileName(_fileName);
	/*
	 * If incoming file type is ascii -- write this type of file
	 */
	if(vtkASCII==type) w->SetDataModeToAscii();

	return w;
}


void write(vtkSmartPointer<vtkXMLPUnstructuredGridWriter> w, vtkSmartPointer<vtkUnstructuredGrid> g){
	w->SetInput(g);
	w->Write();
}


vtkSmartPointer<vtkUnstructuredGrid> getGrid(shared_ptr<double>& y, int numPoints){
	return getGrid(y.get(),numPoints);
}

vtkSmartPointer<vtkPoints> createVTK_Points(double *yPtr, int numPoints){
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	/*
	 * Add coordinates to grid
	 * This directly uses the pointer to data provided  -- this is the part that
	 * the note above refers to.
	 */
	vtkSmartPointer<vtkDoubleArray> ptsData = vtkSmartPointer<vtkDoubleArray>::New();
	int numComponents=3;
	ptsData->SetNumberOfComponents(numComponents);
	int save=1;
	ptsData->SetArray(yPtr,numPoints*numComponents,save);
	pts->SetData(ptsData);
	return pts;
}

/*
 * NOTE: pointer to data is used directly here;
 * This method lets the user specify data to be held by the array.  The
 * array argument is a pointer to the data.  Set save to 1 to keep the class
 * from deleting the array when it cleans up or reallocates memory.
 * The class uses the actual array provided; it does not copy the data
 * from the suppled array.
 * NOTE: Also -- the supplied coordinates 'y' should be a vector with 3 components
 * per point
 */
vtkSmartPointer<vtkUnstructuredGrid> getGrid(double *yPtr, int numPoints){

	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

	/*
	 * Add coordinates to grid
	 * This directly uses the pointer to data provided  -- this is the part that
	 * the note above refers to.
	 */
	vtkSmartPointer<vtkDoubleArray> ptsData = vtkSmartPointer<vtkDoubleArray>::New();
	int numComponents=3;
	ptsData->SetNumberOfComponents(numComponents);
	int save=1;
	ptsData->SetArray(yPtr,numPoints*numComponents,save);
	pts->SetData(ptsData);

	/*
	 * Creates a grid of "points" and "cells"
	 */
	vtkSmartPointer<vtkCellArray> cells = getCellArray(numPoints);
	vtkSmartPointer<vtkUnstructuredGrid> grid = getGrid(pts,cells);

	return grid;
}

vtkSmartPointer<vtkCellArray> createVTK_quadCells(size_t* vLinks, int numCells){
	vtkSmartPointer<vtkIdTypeArray>  links = vtkSmartPointer<vtkIdTypeArray>::New();
	vtkIdType npe(4);
	links->SetNumberOfValues((npe+1)*numCells);
	for(int cell=0;cell<numCells;cell++){
		links->SetValue((npe+1)*cell,npe);
		for(int n=0;n<npe;n++){
			int node = vLinks[npe*cell+n];
			links->SetValue((npe+1)*cell+(n+1),node);
		}
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->SetCells(numCells,links);
	cells->SetNumberOfCells(numCells);
	return cells;
}

vtkSmartPointer<vtkUnstructuredGrid>  getGrid(const vtkSmartPointer<vtkPoints>& x, const vtkSmartPointer<vtkCellArray>& cells, VTKCellType type){

	/**
	 * Cells constructed are "vertex" type cells
	 */
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(x);
	grid->SetCells(type,cells);

	return grid;
}

vtkSmartPointer<vtkCellArray> getCellArray(vtkIdType numCells){

	// This is a really simply link array
	// It assumes that the cell id # is same as point id #
	// By construction, the points generated here satisfy this condition

	vtkSmartPointer<vtkIdTypeArray>  links = vtkSmartPointer<vtkIdTypeArray>::New();
	vtkIdType onePointPerCell(1);
	links->SetNumberOfValues(2*numCells);
	for(int cell=0;cell<numCells;cell++){
		links->SetValue(2*cell,onePointPerCell);
		links->SetValue(2*cell+1,cell);
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->SetCells(numCells,links);

	return cells;
}

/*
 * Initial hack for post processing step; call this from output manager
 * after all data has been set on 'grid'
 */
void expandRingPostProcess(double current_time, vtkSmartPointer<vtkUnstructuredGrid> grid, int myRank){

	vtkKdTreePointLocator* kdTree = vtkKdTreePointLocator::New();
	kdTree->SetDataSet(grid);
	double ri(16.0), ro(17.0), h(1.0);
	double horizon = 1.01/2;
	double rho(2.71e-3);
	/*
	 * Center point of search: Find points within sphere of radius = horizon
	 */
	double xC[] = {(ri+ro)/2.0,0.0,h/2.0};

	/*
	 * Search
	 */
	vtkIdList* kdTreeList = vtkIdList::New();
	kdTree->FindPointsWithinRadius(horizon, xC, kdTreeList);

	/*
	 * number of points to average over
	 */
	size_t numPoints = kdTreeList->GetNumberOfIds();

	/*
	 * Compute 'volume' weighted sum over points found
	 */
	double sum[]={0.0,0.0,0.0,0.0};
	vtkPoints* points = grid->GetPoints();
	vtkPointData* pointData = grid->GetPointData();
	vtkDataArray* vData = pointData->GetVectors("Velocity");
	vtkDataArray* uData = pointData->GetVectors("Displacement");
	vtkDataArray* volData = pointData->GetScalars("Volume");
	vtkDataArray* lambdaData = pointData->GetScalars("Lambda");
	/*
	 * Weighted volume averages
	 */
	for(size_t n=0;n<numPoints;n++){
		vtkIdType localId = kdTreeList->GetId(n);
		double *xyz = points->GetPoint(localId);
		double x = xyz[0];
		double y = xyz[1];
		double theta = atan2(y,x);
		double *u = uData->GetTuple3(localId);
		double ur = u[0]*cos(theta) + u[1]*sin(theta);
		double volume = volData->GetTuple1(localId);
		double lambda = lambdaData->GetTuple1(localId);
		sum[0]+= volume;
		sum[1]+= ur*volume;
		sum[2]+= lambda*volume;
	}
	kdTree->Delete();
	kdTreeList->Delete();

	/*
	 * Calculate kinetic energy
	 * HERE WE LOOP OVER ALL OWNED POINTS on processor
	 */
	for(vtkIdType n=0;n<grid->GetNumberOfPoints();n++){
		vtkIdType localId = n;
		double *v = vData->GetTuple3(localId);
		double volume = volData->GetTuple1(localId);
		sum[3] += (v[0]*v[0]+v[1]*v[1]+v[2]*v[2]) * volume;
	}


	/*
	 * Sum across all processors
	 */
	int four(4);
	double totalSum[]={0.0,0.0,0.0,0.0};
	MPI_Allreduce(sum,totalSum,four,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	/*
	 * Compute averages
	 */
	double urAvg = totalSum[1]/totalSum[0];
	double lambdaAvg = totalSum[2]/totalSum[0];
	double KE = 0.5 * rho * totalSum[3];

	/*
	 * write out result on processor '0'
	 */

	if(0==myRank){
		std::ofstream out;
		out.open("meanValues.dat",std::ios_base::app);
		out.precision(9);
		out << current_time << ", " << urAvg << ", " << lambdaAvg << ", " << KE << std::endl;
		out.close();
	}

}

/**
 * Computes volume of a single HEX8: Smokin Fast!
 * Perhaps this calculation is done with very nearly the least
 * number of flops.
 * Note that there can only be 8 incoming points
 * Order of points must be according to VTK convention (also same as exodus)  */
double compute_hex8_volume(vtkPoints* points){
	using namespace UTILITIES;
	/*
	 * swap points 2 & 3
	 * swap points 6 & 7
	 */
	int map[] = {0,1,3,2,4,5,7,6};
	double xyz[3];
	Vector3D x[8];
	for(int p=0;p<8;p++){
		points->GetPoint(map[p],xyz);
		x[p]=Vector3D(xyz);
	}
	Minus m;
	Vector3D x17(m(x[7],x[1]));
	Vector3D x27(m(x[7],x[2]));
	Vector3D x47(m(x[7],x[4]));
	Vector3D x06(m(x[6],x[0]));
	Vector3D x05(m(x[5],x[0]));
	Vector3D x03(m(x[3],x[0]));

	Plus p;
	Vector3D A(p(x17,x06));
	double v1 = UTILITIES::scalar_triple_product(A,x27,x03);
	Vector3D B(p(x27,x05));
	double v2 = UTILITIES::scalar_triple_product(x06,B,x47);
	Vector3D C(p(x47,x03));
	double v3 = UTILITIES::scalar_triple_product(x17,x05,C);
	return (v1+v2+v3)/12;
}


} // PdVTK
