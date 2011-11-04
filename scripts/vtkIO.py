'''
Created on Aug 27, 2010

@author: jamitch
'''
import vtk.io as io
from vtk import vtkKdTreePointLocator
from xml.dom.minidom import parse
from operator import itemgetter

def GetReader(fileName):
    """Creates and returns grid reader for input file
    
    Input:  string filename -- file type should be 'pvtu' 
    Output: vtkXMLPUnstructuredGridReader
    """
    reader=io.vtkXMLPUnstructuredGridReader()
    reader.SetFileName(fileName)
    reader.Update()
    return reader

def GetGrid(fileName):
	"""Returns vtkUnstructuredGrid
	
	Input:  string filename -- file type should be 'pvtu'
	Output: vtkUntructuredGrid
	"""
	reader=GetReader(fileName)
	grid = reader.GetOutputAsDataSet()
	return grid
    
def GetPointTuples(grid):
	"""Returns point tuples associated with grid
	
	Input:  vtkUnstructuredGrid --> vtkPointSet
	Output: vtkDataArray
	"""

	xyz = grid.GetPoints().GetData()
	return xyz

def GetKdTreePointLocator(grid):
	"""Returns a vtk KdTreePointLocator for grid
	
	Input: vtkUnstructuredGrid
	Output: vtkKdTree
	"""
	kdTree = vtkKdTreePointLocator()
	kdTree.SetDataSet(grid);
	return kdTree
	

def GetData(grid, relationType='point'):
	"""Returns cellData associated with grid
	
	Input: vtkUnstructuredGrid --> vtkDataSet
	Input: relationType -- allowed: 'point' or 'cell'
	Output: data -- type is vtkPointData OR vtkCellData
	"""
	if 'point'==relationType:
		data=grid.GetPointData()
	elif 'cell'==relationType:
		data=grid.GetCellData()
	else:
		raise TypeError("relationType must be: \'point\' or \'cell\'; input value = "+repr(relationType))

	return data

def GetFieldTuple(fieldName,data):
	"""Returns field data tuples on cellData 
	
	Input:  string fieldName
	Input:  data -- type is vtkPointData OR vtkCellData
	Output: vtkDoubleArray
	"""
	return data.GetVectors(fieldName)

def GetDataFieldNames(data):
	"""List of field names on data
	
	Input:  data -- type is vtkPointData OR vtkCellData 
	Output: python list
	"""
	return [data.GetArrayName(i) for i in range(data.GetNumberOfArrays())]

def GetDataTuplesDictionary(data):
	"""Returns a dictionary : fieldName,fieldTuple
	
	Input:  data -- type is vtkPointData OR vtkCellData
	Output: dictionary (keys=field names, values=vtkDoubleArray)
	"""
	return dict([(data.GetArrayName(i),data.GetArray(i)) for i in range(data.GetNumberOfArrays())])
	
def GetTimeCollection(filename):
	"""Opens and reads "pvd" file; Returns List of Tuples: [(timestep, filename)] 
	
	Input: string filename of "pvd" file to read
	Preconditions: file must be an xml paraview collection file
	Output: return list is sorted by time 
	"""

	dom=parse(filename)
	vtkFileElement=dom.documentElement
	collectionNodeList=vtkFileElement.getElementsByTagName("Collection")
	collectionElement=collectionNodeList.item(0)
	dataSetNodeList=collectionElement.getElementsByTagName("DataSet")
	d=[(float(dataSetNodeList.item(i).getAttribute("timestep")),dataSetNodeList.item(i).getAttribute("file")) for i in range(dataSetNodeList.length)]
	return sorted(d,key=itemgetter(0))
	
def GetPieceFileNames(filename):
	"""Parses input file name and extract piece part names;  
	
	Input: string filename of ".pvtu" file to read
	Preconditions: file must be an xml paraview .pvtu file
	Output: return list containing piece part 'Source' filenames
	"""
	dom=parse(filename)
	vtkFileElement=dom.documentElement
	gridList=vtkFileElement.getElementsByTagName("PUnstructuredGrid")
	grid=gridList.item(0)
	pieceList=grid.getElementsByTagName("Piece")
	return [pieceList.item(i).getAttribute('Source') for i in range(pieceList.length)]
	
	
	

	

