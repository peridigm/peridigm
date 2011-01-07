import sys
sys.path.append('/home/jamitch/pimp_projects/pimp/operators/unitTest/eigen1D')
sys.path.append('/home/jamitch/peridigm_projects/peridigm.votd/scripts')

from vtkIO import *
from math import fabs, pi, sqrt, cos, sin
import matplotlib.pyplot as plot
from operator import itemgetter

class Demo:
	"""
	
	"""
	f='peridigmEigen1D.pvd'
	numPeriods=1
	periodNumber=1
	cylinderLength = 15.0
	E = 68.9e3
	rho = 2.7e-3

	def GetReader(self,fileName):
	    """Creates and returns grid reader for input file
	    
	    Input:  string filename -- file type should be 'vtu' 
	    Output: vtkXMLUnstructuredGridReader
	    """
	    reader=io.vtkXMLUnstructuredGridReader()
	    reader.SetFileName(fileName)
	    reader.Update()
	    return reader

	def GetGrid(self,fileName):
		"""Returns vtkUnstructuredGrid
		
		Input:  string filename -- file type should be 'vtu'
		Output: vtkUntructuredGrid
		"""
		reader=self.GetReader(fileName)
		grid = reader.GetOutputAsDataSet()
		return grid

	def __init__(self):
		self.files = GetTimeCollection(Demo.f)
		self.beta=self.periodNumber*pi/self.cylinderLength
		self.omega = self.beta*sqrt(Demo.E/Demo.rho)
		self.period = 2.0*pi/self.omega
		self.v0 = 1.0
		# default plot
		self.makeAnalytical=self.MakeAnalyticalDisplacement;
	
	def GetFiles(self,start,stop,increment):
		return [self.files[i][1] for i in range(start,stop,increment)]

	def GetFileTimes(self,start,stop,increment):
		return [self.files[i][0] for i in range(start,stop,increment)]
	
	def MakeAnalyticalVelocity(self,time):
		c1=self.v0/self.omega
		return lambda z: c1 * self.omega*cos(self.omega*time) * cos(self.beta*z)

	def MakeAnalyticalDisplacement(self,time):
		c1=self.v0/self.omega
		return lambda z: c1 * sin(self.omega*time) * cos(self.beta*z)
			
	def Plot(self,start, stop, increment):
		files=self.GetFiles(start,stop,increment)
		times=self.GetFileTimes(start,stop,increment)
		file0=files[0]
		grid0=grid=self.GetGrid(file0)
		xyz0=GetPointTuples(grid0)
		for i,f in enumerate(files):
			[z,uz]=self.GetAxisData(f,xyz0)
			line,=plot.plot(z,uz)
			t=times[i]
			u=self.makeAnalytical(t)
			zSkip=[z[j] for j in range(0,len(z),50)]
			ua=[u(jZ) for jZ in zSkip]
			plot.plot(zSkip,ua,mfc=line.get_color(),marker='o', linewidth=0.0)
		plot.xlim(0,15.0)
	
	def GetAxisXYZ(self):
	# any file in time will do
		file=self.GetFiles(0,1,1)[0]
    # This calls local GetGrid() -- NOT the vtkIO.GetGrid()
		grid=self.GetGrid(file)
		xyz=GetPointTuples(grid)
		numPoints = xyz.GetNumberOfTuples()
		axis=[]
		for i in range(numPoints):
			t = xyz.GetTuple3(i)
			x = t[0]
			y = t[1]
			if fabs(x) < 1.0e-5 and fabs(y) < 1.0e-5:
				axis.append(t)
		return sorted(axis,key=itemgetter(2))

	def GetAxisData(self,file,xyz):
		p=[]
		grid=self.GetGrid(file)
		pointData = grid.GetPointData()
		coordZ = pointData.GetVectors("Z_Position")
		numPoints = xyz.GetNumberOfTuples()
		for i in range(numPoints):
			t = xyz.GetTuple3(i)
       # current position along z-axis
			yZ = coordZ.GetTuple1(i)
			x = t[0]
			y = t[1]
			z = t[2]
       # displacement along z-axis
			uZ = yZ - z
			if fabs(x) < 1.0e-5 and fabs(y) < 1.0e-5:
				p.append((z,uZ))
		q=sorted(p,key=itemgetter(0))
		numPoints=len(q)
		px=[q[i][0] for i in range(numPoints)]
		py=[q[i][1] for i in range(numPoints)]
		return [px,py]


