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
vtkSmartPointer<vtkUnstructuredGrid> getGrid(double *y, int numPoints);

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
void writeField
(
		vtkSmartPointer<vtkUnstructuredGrid>& g,
		const char* name_null_terminated,
		std::size_t degree,
		std::size_t size,
		T* data
) {
	vtkCellData *cellData = g->GetCellData();
	/*
	 * Add field to grid
	 */
	vtkSmartPointer<vtkDoubleArray> cellField = vtkSmartPointer<vtkDoubleArray>::New();
	cellField->SetName(name_null_terminated);
	cellField->SetNumberOfComponents(degree);
	int save=1;
	cellField->SetArray(data,size,save);
	cellData->AddArray(cellField);
}

template<class T>
void writeField(vtkSmartPointer<vtkUnstructuredGrid>& g, Field_NS::Field<T> field) {
	const char* name = field.getLabel().c_str();
	std::size_t degree = field.getLength();
	std::size_t size = field.getArray().getSize();
	T* data = field.getArray().get();
	writeField(g,name,degree,size,data);
}

template<class T>
void writeField(vtkSmartPointer<vtkUnstructuredGrid>& g, const Field_NS::FieldSpec& spec, T *data) {
	const char* name = spec.getLabel().c_str();
	std::size_t degree =spec.getLength();
	std::size_t size = degree*g->GetNumberOfCells();
	writeField(g,name,degree,size,data);
}

} // PdVTK

#endif /* PDVTK_H_ */
