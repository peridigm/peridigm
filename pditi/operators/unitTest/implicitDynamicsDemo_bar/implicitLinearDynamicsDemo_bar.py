import sys
sys.path.append('/home/jamitch/pditi_projects/pditi/operators/unitTest')
sys.path.append('/home/jamitch/peridigm_projects/peridigm.build/scripts')

from vtkIO import *
from math import fabs, pi, sqrt, cos, sin
import matplotlib.pyplot as plot
from operator import itemgetter

class Demo:
	"""
	
	"""
	f='utPdImp_implicitLinearDynamicsDemo_bar.pvd'
	numStepsPerPeriod=50
	numPeriods=1
	periodNumber=1
#	cylinderLength = 14.823529411881999
	cylinderLength = 15.0
	E = 68.9e3
	rho = 2.7e-3
	
	def __init__(self):
		self.files = GetTimeCollection(Demo.f)
		self.beta=self.periodNumber*pi/self.cylinderLength
		self.omega = self.beta*sqrt(Demo.E/Demo.rho)
		self.period = 2.0*pi/self.omega
		self.dt = self.period/Demo.numStepsPerPeriod
		self.v0 = 1.0
		# default plot
		self.makeAnalytical=self.MakeAnalyticalDisplacement;
		self.field="Displacement"
	
	def GetFiles(self,start,stop,increment):
		return [self.files[i][1] for i in range(start,stop,increment)]

	def GetFileTimes(self,start,stop,increment):
		return [self.files[i][0] for i in range(start,stop,increment)]

	def GetTimes(self,start,stop,increment):
		return [i*self.dt for i in range(start,stop,increment)]
	
	def MakeAnalyticalVelocity(self,time):
		c1=self.v0/self.omega
		return lambda z: c1 * self.omega*cos(self.omega*time) * cos(self.beta*z)

	def MakeAnalyticalDisplacement(self,time):
		c1=self.v0/self.omega
		return lambda z: c1 * sin(self.omega*time) * cos(self.beta*z)
		
	def PlotD(self,start, stop, increment):
		self.makeAnalytical=self.MakeAnalyticalDisplacement;
		self.field="Displacement"
		self.Plot(start,stop,increment)

	def PlotV(self,start, stop, increment):
		self.makeAnalytical=self.MakeAnalyticalVelocity;
		self.field="v"
		self.Plot(start,stop,increment)
	
	def Plot(self,start, stop, increment):
		files=self.GetFiles(start,stop,increment)
		times=self.GetTimes(start,stop,increment)
		fTimes=self.GetFileTimes(start,stop,increment)
		for i,f in enumerate(files):
			[z,uz]=self.GetAxisData(f,self.field)
			line,=plot.plot(z,uz)
			t=times[i]
			u=self.makeAnalytical(t)
			zSkip=[z[j] for j in range(0,len(z),6)]
			ua=[u(jZ) for jZ in zSkip]
		#	plot.plot(zSkip,ua, marker='o',linewidth=0.0)
			plot.plot(zSkip,ua,mfc=line.get_color(),marker='o', linewidth=0.0)
		plot.xlim(0,15.0)
	
	def GetAxisXYZ(self):
	# any file in time will do
		file=self.GetFiles(0,1,1)[0]
		grid=GetGrid(file)
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

	def GetAxisData(self,file,field):
		p=[]
		grid=GetGrid(file)
		xyz=GetPointTuples(grid)
		cellData = GetCellData(grid)
		u = cellData.GetVectors(field)
		numPoints = xyz.GetNumberOfTuples()
		for i in range(numPoints):
			t = xyz.GetTuple3(i)
			d = u.GetTuple3(i)
			x = t[0]
			y = t[1]
			z = t[2]
			if fabs(x) < 1.0e-5 and fabs(y) < 1.0e-5:
				p.append((z,d[2]))
		q=sorted(p,key=itemgetter(0))
		numPoints=len(q)
		px=[q[i][0] for i in range(numPoints)]
		py=[q[i][1] for i in range(numPoints)]
		return [px,py]

# ipython -pylab
# files = GetCollection()
# pFiles = [files[i][1] for i in range(20,25,1)]
# figure()
# hold
# for f in pFiles:
#     [p1,p2]=getAxisData(f)
#     plot(p1,p2)
# savefig("fig.svg",dpi=600,format='svg')


