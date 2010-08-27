/*
 * PdVTK.h
 *
 *  Created on: Nov 14, 2009
 *      Author: awesome
 */

#ifndef PDVTK_H_
#define PDVTK_H_
#include "vtkUnstructuredGrid.h"
#include "vtkSmartPointer.h"
#include <vtkDoubleArray.h>
#include "vtkCellData.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include "PdGridData.h"
#include <tr1/memory>
#include "Field.h"
#include <deque>
#include <string>

using std::deque;
using std::string;

namespace PdVTK {

enum VTK_FILE_TYPE { vtkASCII=0, vtkBINARY };

vtkSmartPointer<vtkUnstructuredGrid> getGrid(shared_ptr<double>& y, int numPoints);
vtkSmartPointer<vtkXMLPUnstructuredGridWriter> getWriter(const char* _fileName, int numProcs, int rank, VTK_FILE_TYPE type=vtkBINARY);
void write(vtkSmartPointer<vtkXMLPUnstructuredGridWriter> w, vtkSmartPointer<vtkUnstructuredGrid> g);

class CollectionWriter {
public:
	CollectionWriter(const char* _fileName, int numProcs, int rank, VTK_FILE_TYPE t=vtkBINARY);
	void writeTimeStep(double t, vtkSmartPointer<vtkUnstructuredGrid> grid);
	void close();
private:
	string getPVTU_fileName(int index, const char* _fileName) const;
	const char* fileName;
	deque<double> times;
	vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writer;
};

/*
 * DEMO On how to write a field
 * NOTE: field passed in must go out of scope AFTER the vtkUnstructuredGrid grid "g" otherwise
 * vtk writer may try to access the data stored on "g" but it would have been deleted
 */
template<class T>
void writeField(vtkSmartPointer<vtkUnstructuredGrid> g, Field_NS::Field<T> field) {
	vtkCellData *cellData = g->GetCellData();
//	vtkIdType numCells = g->GetNumberOfCells();
	/*
	 * Add field to grid
	 */
	vtkSmartPointer<vtkDoubleArray> cellField = vtkSmartPointer<vtkDoubleArray>::New();
	cellField->SetNumberOfComponents(field.getLength());
	cellField->SetName(field.getLabel().c_str());
	int save=1;
//	cellField->SetArray(field.get(),numCells*field.getLength(),save);
	Pd_shared_ptr_Array<T> fieldArray = field.getArray();
	cellField->SetArray(fieldArray.get(),fieldArray.getSize(),save);
	cellData->AddArray(cellField);
}

} // PdVTK

#endif /* PDVTK_H_ */
