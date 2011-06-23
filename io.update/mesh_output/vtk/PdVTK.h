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
#include <vtkIntArray.h>
#include <vtkCharArray.h>
#include <vtkCellType.h>
#include "vtkPointData.h"
#include "vtkIdList.h"
#include "vtkKdTree.h"
#include "vtkKdTreePointLocator.h"
#include "vtkXMLPUnstructuredGridWriter.h"
#include <tr1/memory>
#include "Field.h"
#include <deque>
#include <string>

using std::deque;
using std::string;
using std::tr1::shared_ptr;

class vtkPoints;

namespace PdVTK {

enum VTK_FILE_TYPE { vtkASCII=0, vtkBINARY };

template<typename T>
struct vtk_trait {};

template<>
struct vtk_trait<char> {
   typedef vtkSmartPointer<vtkCharArray> vtk_type;
};

template<>
struct vtk_trait<int> {
   typedef vtkSmartPointer<vtkIntArray> vtk_type;
};

template<>
struct vtk_trait<double> {
   typedef vtkSmartPointer<vtkDoubleArray> vtk_type;
};

vtkSmartPointer<vtkUnstructuredGrid> getGrid(shared_ptr<double>& y, int numPoints);

/*
 * TEMPORARY
 */
void expandRingPostProcess(double current_time, vtkSmartPointer<vtkUnstructuredGrid> grid, int myRank);
/*
 * END TEMPORARY
 */
double compute_hex8_volume(vtkPoints* points);
vtkSmartPointer<vtkPoints> createVTK_Points(double *yPtr, int numPoints);
vtkSmartPointer<vtkCellArray> createVTK_quadCells(size_t* vLinks, int numCells);
vtkSmartPointer<vtkUnstructuredGrid> getGrid(shared_ptr<double>& y, int numPoints);
vtkSmartPointer<vtkUnstructuredGrid> getGrid(double *y, int numPoints);
vtkSmartPointer<vtkUnstructuredGrid> getGrid(const vtkSmartPointer<vtkPoints>& x, const vtkSmartPointer<vtkCellArray>& cells, VTKCellType type=VTK_VERTEX);
vtkSmartPointer<vtkXMLPUnstructuredGridWriter> getWriter(const char* _fileName, int numProcs, int rank, VTK_FILE_TYPE type=vtkBINARY);
void write(vtkSmartPointer<vtkXMLPUnstructuredGridWriter> w, vtkSmartPointer<vtkUnstructuredGrid> g);

class CollectionWriter {
public:
	CollectionWriter(const char* _fileName, int numProcs, int rank, VTK_FILE_TYPE t=vtkBINARY);
	void writeTimeStep(double t, vtkSmartPointer<vtkUnstructuredGrid> grid);
	void close(const string& comment="");
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
template<typename T>
void writeField
(
vtkSmartPointer<vtkUnstructuredGrid>& g,
const char* name_null_terminated,
std::size_t degree,
T *data
) {

	/*
	 * Write "point data"
	 */
	vtkPointData *pointData = g->GetPointData();
	/*
	 * Add field to grid
	 */
	typename vtk_trait<T>::vtk_type field;
	field = vtk_trait<T>::vtk_type::New();

	field->SetName(name_null_terminated);
	field->SetNumberOfComponents(degree);
	int save=1;
	std::size_t size = degree*g->GetNumberOfPoints();
	field->SetArray(data,size,save);
	pointData->AddArray(field);

}


template<class T>
void writeField(vtkSmartPointer<vtkUnstructuredGrid>& g, PdITI::Field_NS::Field<T> field) {
	const char* name = field.getLabel().c_str();
	std::size_t degree = field.getLength();
	T* data = field.get();
	writeField(g,name,degree,data);
}

template<class T>
void writeField(vtkSmartPointer<vtkUnstructuredGrid>& g, const PdITI::Field_NS::FieldSpec& spec, T *data) {
	const char* name = spec.getLabel().c_str();
	std::size_t degree =spec.getLength();
	writeField(g,name,degree,data);
}

} // PdVTK

#endif /* PDVTK_H_ */
