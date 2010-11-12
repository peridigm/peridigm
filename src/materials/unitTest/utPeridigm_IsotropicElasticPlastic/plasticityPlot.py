#! /usr/bin/env python

from pylab import *


# load computed data
d=loadtxt("ep.dat")

# time
t=[d[i][0] for i in range(len(d))]
# applied displacement/shear
u=[d[i][1] for i in range(len(d))]
# consistency parameter
c=[d[i][2] for i in range(len(d))]
# force
f=[d[i][3] for i in range(len(d))]

#f1=figure(figsize=(6,6))
f1=figure(figsize=(12,6))
ax1=f1.add_subplot(121,autoscale_on=False,xlim=(-0.05,3),ylim=(-.002,.006))
# this adjusts sub plot margins so that in this case axis labels are not cut off
f1.subplots_adjust(bottom=.12,left=.18)
line1=ax1.plot(t,u,linewidth=2.0, color='b')
ax1.plot([t[0]],[u[0]],'o',color='r')
ax1.plot([t[50]],[u[50]],'o',color='g')
ax1.plot([t[100]],[u[100]],'o',color='y')

# LOADING ANNOTATION
# text
ax1.annotate("Loading",
xy=(1.0,0.0),xycoords='data',
fontsize=20,
horizontalalignment='center', verticalalignment='bottom',
)

# left arrow
ax1.annotate("",
xytext=(1.0,0.0007),xycoords='data',
fontsize=30,xy=(.4,.002),
arrowprops=dict(facecolor='black',
shrink=0.05,width=1, headwidth=4),
horizontalalignment='center', verticalalignment='bottom',
)

# right arrow
ax1.annotate("",
xytext=(1.55,0.00035),xycoords='data',
fontsize=30,xy=(2.5,-.0001),
arrowprops=dict(facecolor='black',
shrink=0.05,width=.9, headwidth=3),
horizontalalignment='right', verticalalignment='center',
)

# UNLOADING ANNOTATION
# arrow
ax1.annotate("",
xy=(1.4,0.003),xycoords='data',
fontsize=30,xytext=(1.8,.0045),
arrowprops=dict(facecolor='black',
shrink=0.05,width=1, headwidth=4),
horizontalalignment='left', verticalalignment='bottom',
)
# text
ax1.annotate("Unloading",
xy=(1.4,.0045),xycoords='data',
fontsize=20,
)

ax1.annotate(r'$\times 10^{-3}$',xy=(0,1),xycoords='axes fraction',fontsize=20)
ax1.set_xticks([0,1,2,3])
ax1.set_xticklabels([r"$0$",r"$1$",r"$2$",r"$3$"],fontsize=20)
ax1.set_yticks([-.002,0,.002,.004,.006])
ax1.set_yticklabels([r"$-2$",r"$0$",r"$2$",r"$4$",r"$6$"],fontsize=20)
#title("Prescribed Kinematics",fontsize=20)
xlabel("Load step parameter",fontsize=20)
ylabel("Prescribed Shear "+r'$\gamma$',fontsize=20)
figName="/home/awesome/Desktop/Docs/peridynamics/figures/twoPointPlasticityLoading.eps"
#f1.savefig(figName,dpi=600,pad_inches=0.0)

#f2=figure(figsize=(6,6))
f2=f1
ax2=f2.add_subplot(122,autoscale_on=False,xlim=(-.0015,.0055),ylim=(-1000,1000))
f2.subplots_adjust(bottom=.12,left=.18)
ax2.annotate(r'$\times 10^{3}$',xy=(0,1),xycoords='axes fraction',fontsize=20)
ax2.annotate(r'$\times 10^{-3}$',xy=(0.95,-0.1),xycoords='axes fraction',fontsize=20)
# add background grid lines on "both" major and minor ticks
ax2.grid(True,which='both')
line2=ax2.plot(u,f,linewidth=2.0, color='b')
ax2.plot([u[0]],[f[0]],'o',color='r')
ax2.plot([u[50]],[f[50]],'o',color='g')
ax2.plot([u[100]],[f[100]],'o',color='y')

#title("Time Integration of Plasticity Model",fontsize=20)
ax2.set_xticks([-.001,-.0005,0,.0005,.001,.0015,.002,.0025,.003,.0035,.004,.0045,.005])
ax2.set_xticklabels([r"$-1$","",r"$0$","",r"$1$","",r"$2$","",r"$3$","",r"$4$","",r"$5$"],fontsize=20)
ax2.set_yticks([-1e3,-.75e3,-.5e3,-.25e3,0,.25e3,.5e3,.75e3,1.0e3])
ax2.set_yticklabels([r"$-1$","",r"$-.5$","",r"$0$","",r"$.5$","",r"$1$"],fontsize=20)

# This places y-axis label on right
#ax2.yaxis.label_position='right'

# Since I can't get the yaxis label positioned correctly, the following does it directly
text(6.0e-3,0,"Internal Force Magnitude",rotation="vertical",fontsize=20,va='center')

#ylabel("Internal Force Magnitude",fontsize=20)
xlabel("Shear Strain  "+r'$\gamma$',fontsize=20)
#figName="/home/awesome/Desktop/Docs/peridynamics/figures/twoPointPlasticityTimeIntegration.eps"
#figName="/home/awesome/Desktop/Docs/peridynamics/figures/twoPointPlasticityCombined.eps"
figName="twoPointPlasticityCombined.eps"
f2.savefig(figName,dpi=600,pad_inches=0.0)

