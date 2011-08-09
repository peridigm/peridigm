#! /usr/bin/env python

from pylab import *

def GetTickLabels(tickValues):
    lab=["$"+str(d)+"$" for d in tickValues]
    return lab
    
# load computed data
d=loadtxt("utPeridigm_twoPointMaxwellRelaxation.dat")

t=[d[i][0] for i in range(len(d))]
u=[d[i][1] for i in range(len(d))]
f=[d[i][2] for i in range(len(d))]

f1=figure(figsize=(11,6))
ax1=f1.add_subplot(121,autoscale_on=False,xlim=(0,4.5),ylim=(0,.3))
xTicks=linspace(0,4,5)
ax1.set_xticks(xTicks)
xTickLabels=GetTickLabels(xTicks)
ax1.set_xticklabels(xTickLabels,fontsize=20)
xminorticks=MultipleLocator(0.5)
xminorticks.view_limits(0.5,4.0)
ax1.xaxis.set_minor_locator(xminorticks)
yTicks=linspace(0,.3,7)
ax1.set_yticks(yTicks)
yTickLabels=GetTickLabels(yTicks)
yminorticks=MultipleLocator(0.025)
ax1.yaxis.set_minor_locator(yminorticks)
#ax1.grid(which='minor')
ax1.set_yticklabels(yTickLabels,fontsize=20)
line1=ax1.plot(t,f,linewidth=2.0, color='b')
ax1.plot([t[0]],[f[0]],'o',color='g')
ax1.plot([t[-1]],[f[-1]],'o',color='g')


ax2=f1.add_subplot(122,autoscale_on=False,xlim=(0,4.5),ylim=(0,1.1e-6))
ax2.set_xticks(xTicks)
xTickLabels=GetTickLabels(xTicks)
ax2.set_xticklabels(xTickLabels,fontsize=20)
xminorticks=MultipleLocator(0.5)
xminorticks.view_limits(0.5,4.0)
ax2.xaxis.set_minor_locator(xminorticks)
yTicks=linspace(0,11*pow(10,-7),12)
ax2.set_yticks(yTicks)
yTickLabels=GetTickLabels(yTicks*pow(10,6))
ax2.set_yticklabels(yTickLabels,fontsize=20)
ax2.annotate(r'$\times 10^{-6}$',xy=(0,1),xycoords='axes fraction',fontsize=20)
line2=ax2.plot(t,u,linewidth=2.0, color='b')

