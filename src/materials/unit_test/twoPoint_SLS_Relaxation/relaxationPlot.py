#! /usr/bin/env python

from pylab import *
from relaxation import standard_linear_solid

files=[(r'$\lambda_i=.01$',"twoPoint_SLS_Elastic.dat"),
       (r'$\lambda_i=.5$',"twoPoint_SLS_Relaxation.dat"),
       (r'$\lambda_i=.99$',"twoPoint_Maxwell_Relaxation.dat")
       ]
    
def GetTickLabels(tickValues):
    lab=["$"+str(d)+"$" for d in tickValues]
    return lab

   
# load computed data
d0=loadtxt(files[0][1])
d1=loadtxt(files[1][1])
d2=loadtxt(files[2][1])

t=[d0[i][0] for i in range(len(d0))]
u=[d0[i][1] for i in range(len(d0))]
f=[d0[i][2] for i in range(len(d0))]
f2=[d1[i][2] for i in range(len(d1))]
f1000=[d2[i][2] for i in range(len(d2))]
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

# make model elastic
modulus=.01
maxwell=standard_linear_solid(alpha,tau_b,modulus)
fa_1=maxwell.relax(ed0,ta)
modulus=.5
# analytical value for standard linear solid
sls=standard_linear_solid(alpha,tau_b,modulus)
fa_2=sls.relax(ed0,ta)

# analytical solution for maxwell case
modulus=.99
sls_elastic=standard_linear_solid(alpha,tau_b,modulus)
fa_3=sls_elastic.relax(ed0,ta)


f1=figure(figsize=(8,7))
ax1=f1.add_subplot(111,autoscale_on=False,xlim=(0,6.5),ylim=(0,.3))
xTicks=linspace(0,6,5)
ax1.set_xticks(xTicks)
xTickLabels=GetTickLabels(xTicks)
ax1.set_xticklabels(xTickLabels,fontsize=20)
xminorticks=MultipleLocator(0.5)
xminorticks.view_limits(0.5,6.0)
ax1.xaxis.set_minor_locator(xminorticks)
yTicks=linspace(0,.3,7)
ax1.set_yticks(yTicks)
yTickLabels=GetTickLabels(yTicks)
yminorticks=MultipleLocator(0.025)
ax1.yaxis.set_minor_locator(yminorticks)

ax1.set_yticklabels(yTickLabels,fontsize=20)
line1=ax1.plot(t,f,linewidth=2.0, color='g')
line12=ax1.plot(t,f2,linewidth=2.0, color='b')
line13=ax1.plot(t,f1000,linewidth=2.0, color='#FF9955')



# plot analytical
# We need to double this force: computing the 
#  force density
for i in range(len(fa_2)):
    ax1.plot([ta[i]],[2*fa_1[i]],'o',color=(4.0/255.0,250.0/255.0,33.0/255.0))
    ax1.plot([ta[i]],[2*fa_2[i]],'o',color=(0.0/255.0,220.0/255.0,251.0/255.0))
    ax1.plot([ta[i]],[2*fa_3[i]],'o',color='#D45500')

ax1.legend([r"$\lambda_i=.01$",r'$\lambda_i=.5$', r'$\lambda_i=.99$','Analytical', 'Analytical', 'Analytical'],loc="upper right",bbox_to_anchor=(0.96, 0.89),shadow=True)

xlabel("Time (miliseconds)",fontsize=20)
ylabel("Force Density Magnitude",fontsize=20)
f1.subplots_adjust(left=.18,top=.95,right=.96)


#ax2=f1.add_subplot(122,autoscale_on=False,xlim=(0,6.5),ylim=(0,1.1e-6))
#ax2.set_xticks(xTicks)
#xTickLabels=GetTickLabels(xTicks)
#ax2.set_xticklabels(xTickLabels,fontsize=20)
#xminorticks=MultipleLocator(0.5)
#xminorticks.view_limits(0.5,6.0)
#ax2.xaxis.set_minor_locator(xminorticks)
#yTicks=linspace(0,11*pow(10,-7),12)
#ax2.set_yticks(yTicks)
#yTickLabels=GetTickLabels(yTicks*pow(10,6))
#ax2.set_yticklabels(yTickLabels,fontsize=20)
#ax2.annotate(r'$\times 10^{-6}$',xy=(0,1),xycoords='axes fraction',fontsize=20)
#line2=ax2.plot(t,u,linewidth=2.0, color='b')
#xlabel("Time (miliseconds)",fontsize=20)
#ylabel("Step displacement input: "+r'$\delta$',fontsize=20)
#f1.subplots_adjust(bottom=.13, wspace=.34)


