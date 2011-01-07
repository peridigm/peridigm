import sys
#sys.path.append('/home/jamitch/pimp_projects/pimp/operators/unitTest/eigen1D')
sys.path.append('/home/jamitch/peridigm_projects/peridigm/scripts')

from vtkIO import *
from math import fabs, pi, sqrt, cos, sin
import matplotlib.pyplot as plot
import matplotlib
from operator import itemgetter
from matplotlib.font_manager import FontProperties
import scipy


class Demo:
	"""
	
	"""
	f='utPdImp_eigen1D.pvd'
	numPeriods=1
	periodNumber=1
	cylinderLength = 15.0
	E = 68.9e3
	rho = 2.7e-3
	numCells=150
	
	def __init__(self):
		self.files = GetTimeCollection(Demo.f)
		self.beta=self.periodNumber*pi/self.cylinderLength
		self.omega = (self.beta)*sqrt(Demo.E/Demo.rho)
		self.period = 2.0*pi/self.omega
		self.v0 = 1.0
		# default plot
		self.makeAnalytical=self.MakeAnalyticalDisplacement;
		self.field="Displacement"

	def plot(self):
		f=plot.figure(figsize=(7.0,7.5))
		gs=matplotlib.gridspec.GridSpec(4,4,hspace=0.5)
		top=plot.subplot(gs[0,:])
		bottom=plot.subplot(gs[1:,:])
		self.pTop(top)
		self.pBottom(bottom)
		f.subplots_adjust(left=.14, right=.93)
		
	def pTop(self,ax):
		#f=plot.figure(figsize=(7,1.5))
       	# top plot of initial condition
		#ax=f.add_subplot(111,autoscale_on=False,xlim=(0,15),ylim=(-1,1))
		ax.set_title('Cylinder Initial condition: $v_0=V_0\cos(\pi z/L)$',fontsize=20)
		ax.annotate('Cylinder axis',xy=(1,-.25),
			xycoords='data',horizontalalignment='left',verticalalignment='top',fontsize=15,color='green')

		ax.set_xlim(0,15)
		ax.set_ylim(-1,1)
		#ax.set_yticklabels([r"$-1.0$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"],fontsize=15,visible=False)
		ax.set_xticks(range(1,15))
		ax.set_xticklabels(["" for i in range(1,15)],fontsize=15)
		# turn off all ticks associated with y-axis
		for tick in ax.yaxis.get_major_ticks():
			tick.label1On=False
			tick.label2On=False
			tick.tick1On=False
			tick.tick2On=False
		# turn off all ticks associated with x-axis
		for tick in ax.xaxis.get_major_ticks():
			tick.label1On=False
			tick.label2On=False
			tick.tick1On=False
			tick.tick2On=False

		#ax.set_xticklabels(["$"+repr(i)+"$" for i in range(1,15)],fontsize=15)
		c=lambda x: cos(self.beta*x)
		z=scipy.arange(0,self.cylinderLength,self.cylinderLength/self.numCells)
		ax.plot([z[i] for i in range(self.numCells)],[c(z[i]) for i in range(self.numCells)])
		ax.plot([0,self.cylinderLength],[0,0])
		# left arrow that points right
		ax.annotate("",
		xytext=(1,0.5),xycoords='data',
		fontsize=30,xy=(4,0.5),
		arrowprops=dict(facecolor='black',
		shrink=0.05,width=1, headwidth=4),
		horizontalalignment='left', verticalalignment='bottom',)
		# right arrow that points left
		ax.annotate("",
		xytext=(14,-0.5),xycoords='data',
		fontsize=30,xy=(11,-0.5),
		arrowprops=dict(facecolor='black',
		shrink=0.05,width=1, headwidth=4),
		horizontalalignment='left', verticalalignment='bottom',)


		#f.subplots_adjust(left=.14, right=.94,bottom=.19)

	def pBottom(self,ax):
		#f=plot.figure(figsize=(7,6))
		#ax=f.add_subplot(111,autoscale_on=False,xlim=(0,15),ylim=(-.001,.001))
		ax.set_xlim(0,15)
		ax.set_ylim(-.001,.001)
		self.Plot((0,6,12,30,36),ax)
		ax.annotate('Symbols: exact local',xy=(0.525,.62),
			xycoords='figure fraction',horizontalalignment='center',fontsize=20)
		ax.set_title('$1$D verification',fontsize=20)
		ax.annotate(r'$\times 10^{-3}$',xy=(0,1),xycoords='axes fraction',fontsize=15)
		ax.set_yticklabels([r"$-1.0$",r"$-0.5$",r"$0$",r"$0.5$",r"$1$"],fontsize=15)
		ax.set_xticks(range(1,15))
		ax.set_xticklabels(["$"+repr(i)+"$" for i in range(1,15)],fontsize=15)
		lines=ax.lines
		lines[0].set_label('$\omega t=0$')
		lines[2].set_label('$\omega t=\pi/4$')
		lines[4].set_label('$\omega t=\pi/2$')
		lines[6].set_label('$\omega t=5\pi/4$')
		lines[8].set_label('$\omega t=3\pi/2$')
		font=FontProperties(size=16)
		plot.legend(loc='center',bbox_to_anchor=(0.9,0.5),prop=font)
		plot.xlabel('Axial coordinate $z$ (mm)',fontsize=20)
		plot.ylabel('Axial displacement (mm)',fontsize=20)
		#f.subplots_adjust(top=.93,left=.14, right=.94,bottom=.11)

	def GetFiles(self,tuplesIndex):
		return [self.files[i][1] for i in tuplesIndex]

	def GetFileTimes(self,tuplesIndex):
		return [self.files[i][0] for i in tuplesIndex]
	
	def MakeAnalyticalVelocity(self,time):
		c1=self.v0/self.omega
		return lambda z: c1 * self.omega*cos(self.omega*time) * cos(self.beta*z)

	def MakeAnalyticalDisplacement(self,time):
		c1=self.v0/self.omega
		return lambda z: c1 * sin(self.omega*time) * cos(self.beta*z)
		
	def PlotD(self,tuplesIndex):
		self.makeAnalytical=self.MakeAnalyticalDisplacement;
		self.field="Displacement"
		self.Plot(tuplesIndex)

	def PlotV(self,tuplesIndex):
		self.makeAnalytical=self.MakeAnalyticalVelocity;
		self.field="v"
		self.Plot(tuplesIndex)
	
	def Plot(self,tuplesIndex, ax):
		files=self.GetFiles(tuplesIndex)
		times=self.GetFileTimes(tuplesIndex)
		for i,f in enumerate(files):
			skip=8
			[z,uz]=self.GetAxisData(f,self.field)
			line,=ax.plot(z,uz)
			t=times[i]
			u=self.makeAnalytical(t)
			zSkip=[z[j] for j in range(2,len(z),skip)]
			ua=[u(jZ) for jZ in zSkip]
			ax.plot(zSkip,ua,mfc=line.get_color(),marker='o', linewidth=0.0)
		ax.legend()
#		plot.xlim(0,15.0)
	
	def GetAxisXYZ(self):
	# any file in time will do
		file=self.GetFiles((0,1,1))[0]
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
		cellData = GetData(grid)
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


