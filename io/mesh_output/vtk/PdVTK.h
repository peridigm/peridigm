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
void writeField(vtkSmartPointer<vtkUnstructuredGrid>& g, Field_NS::Field<T> field) {
	const char* name = field.getLabel().c_str();
	std::size_t degree = field.getLength();
	T* data = field.get();
	writeField(g,name,degree,data);
}

template<class T>
void writeField(vtkSmartPointer<vtkUnstructuredGrid>& g, const Field_NS::FieldSpec& spec, T *data) {
	const char* name = spec.getLabel().c_str();
	std::size_t degree =spec.getLength();
	writeField(g,name,degree,data);
}

} // PdVTK

#endif /* PDVTK_H_ */
