'''
Created on Aug 27, 2010

@author: jamitch
'''
import vtk.io as io
from xml.dom.minidom import parse

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
	
def GetCellData(grid):
	"""Returns cellData associated with grid
	
	Input:  vtkUnstructuredGrid --> vtkDataSet
	Output: vtkCellData
	"""
	cellData=grid.GetCellData()
	return cellData

def GetFieldTuple(fieldName,cellData):
	"""Returns field data tuples on cellData 
	
	Input:  string fieldName
	Input:  vtkCellData cellData
	Output: vtkDoubleArray
	"""
	return cellData.GetVectors(fieldName)

def GetCellDataFieldNames(cellData):
	"""List of field names on cellData
	
	Input:  vtkCellData cellData
	Output: python list
	"""
	return [cellData.GetArrayName(i) for i in range(cellData.GetNumberOfArrays())]

def GetCellDataTuplesDictionary(cellData):
	"""Returns a dictionary : fieldName,fieldTuple
	
	Input:  vtkCellData cellData
	Output: dictionary (keys=field names, values=vtkDoubleArray)
	"""
	return dict([(cellData.GetArrayName(i),cellData.GetArray(i)) for i in range(cellData.GetNumberOfArrays())])

def GetTimeCollection(filename):
	"""Opens and reads "pvd" file; Returns dictionary (filenames, timestep) 
	
	Input: string filename of "pvd" file to read
	Preconditions: file must be an xml paraview collection file
	Output: dictionary with key = "pvtu" file, value=timestep
	"""
	dom=parse(filename)
	vtkFileElement=dom.documentElement
	collectionNodeList=vtkFileElement.getElementsByTagName("Collection")
	collectionElement=collectionNodeList.item(0)
	dataSetNodeList=collectionElement.getElementsByTagName("DataSet")
	d=dict([(dataSetNodeList.item(i).getAttribute("file"),dataSetNodeList.item(i).getAttribute("timestep")) for i in range(dataSetNodeList.length)])
	return d
	
	

