#! /usr/bin/env python

from pylab import *

def twoPoint_Analytical(alpha, tau, tau_b, ed0, time_steps):
    tau_star=1/tau_b - 1/tau
    td0=alpha*ed0
    if tau == tau_b:
        td_long=0.0
    else:
        td_long=alpha * ed0 * tau_b / tau_star
    print "tau_star = ", tau_star
    analytical=lambda t: (td0-td_long) * exp(-t/tau_b) + td_long
    td=[analytical(t) for t in time_steps]
    return td
    
def GetTickLabels(tickValues):
    lab=["$"+str(d)+"$" for d in tickValues]
    return lab

# this is may or may not be required to export 'eps' files
# this is not required for 'svg'
# Also check or matplotlibrc file -->~/.matplotlib/matplotlibrc
#from matplotlib import rc
#rc('text',usetex=True)
   
# load computed data
d=loadtxt("utPeridigm_twoPointMaxwellRelaxation.dat")
d2=loadtxt("utPeridigm_twoPoint_SLS_Relaxation.dat")

t=[d[i][0] for i in range(len(d))]
u=[d[i][1] for i in range(len(d))]
f=[d[i][2] for i in range(len(d))]
f2=[d2[i][2] for i in range(len(d))]
ta=[t[i] for i in range(0,len(t),5)]

# Analytical solution for single bond relaxation
E=68.9e3
nu=.33
G = E / (1+nu) / 2
# weighted volume
m=2 
# input step -- shear
my_gamma=1.0e-6
alpha = 15.0 * G / m
ed0=my_gamma/sqrt(2)
tau_b=2.0

# analytical solution for maxwell case
tau=tau_b
fa=twoPoint_Analytical(alpha,tau,tau_b,ed0,ta);



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
line12=ax1.plot(t,f2,linewidth=2.0, color='g')
ax1.plot([t[0]],[f[0]],'o',color='g')
ax1.plot([t[-1]],[f[-1]],'o',color='g')
xlabel("Time (miliseconds)",fontsize=20)
ylabel("Scalar force state: "+ r'$t^d$',fontsize=20)


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
xlabel("Time (miliseconds)",fontsize=20)
ylabel("Step displacement input: "+r'$\delta$',fontsize=20)

f1.subplots_adjust(bottom=.13, wspace=.34)