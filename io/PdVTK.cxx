/*
 * PdVTK.cxx
 *
 *  Created on: Nov 14, 2009
 *      Author: awesome
 */
#include "PdVTK.h"
#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkSmartPointer.h"
#include "vtkUnstructuredGrid.h"
#include "vtkVertex.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include <vtkDoubleArray.h>
#include "PdQuickGrid.h"
#include <iostream>
#include <fstream>


namespace PdVTK {


using PdQuickGrid::PdQPointSet1d;
using PdQuickGrid::PdQRing2d;
/*
 * These are private
 */
vtkSmartPointer<vtkUnstructuredGrid> getGrid(const vtkSmartPointer<vtkPoints>& x, const vtkSmartPointer<vtkCellArray>& cells);
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
 * This function writes the collection file that allows paraview to open and read all time steps
 */
void CollectionWriter::close() {
	string outFileName(fileName);
	outFileName += ".pvd";
	std::fstream fStream(outFileName.c_str(), fstream::out);

	fStream << "<?xml version=\"1.0\"?>\n"
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
	/*
	 * Add rank to each point; this should be the rank based upon original construction of writer (see getWriter)
	 */
	int pieceNum = w->GetStartPiece();
	vtkSmartPointer<vtkIdTypeArray> cellRanks = vtkSmartPointer<vtkIdTypeArray>::New();
	cellRanks->SetNumberOfComponents(1);
	vtkCellData *cellData = g->GetCellData();
	vtkIdType numCells = g->GetNumberOfCells();
	cellRanks->SetName("rank");
	for(int p=0;p<numCells;p++)
		cellRanks->InsertNextValue(pieceNum);

	cellData->AddArray(cellRanks);
	w->SetInput(g);
	w->Write();
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
vtkSmartPointer<vtkUnstructuredGrid> getGrid(shared_ptr<double>& y, int numPoints){
	// Set points and cells
	// note number of points is same as number of cells
	vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
	int numCells = numPoints;
	/*
	 * Add coordinates to grid
	 * This directly uses the pointer to data provided  -- this is the part that
	 * the note above refers to.
	 */
	vtkSmartPointer<vtkDoubleArray> ptsData = vtkSmartPointer<vtkDoubleArray>::New();
	int numComponents=3;
	ptsData->SetNumberOfComponents(numComponents);
	int save=1;
	ptsData->SetArray(y.get(),numCells*numComponents,save);
	pts->SetData(ptsData);

	vtkSmartPointer<vtkCellArray> cells = getCellArray(numCells);
	vtkSmartPointer<vtkUnstructuredGrid> grid = getGrid(pts,cells);

	return grid;
}


vtkSmartPointer<vtkUnstructuredGrid>  getGrid(const vtkSmartPointer<vtkPoints>& x, const vtkSmartPointer<vtkCellArray>& cells){

	/**
	 * Cells constructed are "vertex" type cells
	 */
	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(x);
	grid->SetCells(VTK_VERTEX,cells);

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

} // PdVTK
