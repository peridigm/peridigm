#! /usr/bin/env python

import sys
import os
from subprocess import Popen

if __name__ == "__main__":

    executable = sys.argv[-1]
    base_name = executable.split('/')[-1]
    logfile = open(base_name + ".log", 'w')

    if not os.path.exists(executable):
        logfile.write("\nError:  " + executable + " not found!\n\n")
        sys.exit(1)

    # look for existing table and data files
    files=["./table_1.tex", "./ut_bondVolumeConvergenceStudy.dat"]
    for f in files:
        if os.path.exists(f):
            logfile.write("\nRemoving existing output file "+f+"\n")
            os.remove(f)

    command = []
    for i in range(len(sys.argv)-1):
        command.append(sys.argv[1+i])

    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if 0<return_code:
        logfile.write("\nExecution of unit test failed\n\n")
        sys.exit(1)
        
    # load computed and gold data
    dFile=open('ut_bondVolumeConvergenceStudy.dat','r')
    dGFile=open('ut_bondVolumeConvergenceStudy.gold.dat','r')
    dLines=dFile.readlines()
    dGLines=dGFile.readlines()
    if len(dLines) != len(dGLines):
        logfile.write("\nError: Files have different number of lines!")
        sys.exit(1)
    dStr=[line.split(" ") for line in dLines]
    #        n           h           m                ed2
    d=[(int(dd[0]),float(dd[1]),float(dd[2]), float(dd[3])) for dd in dStr]
    dGStr=[line.split(" ") for line in dGLines]
    dG=[(int(dd[0]),float(dd[1]),float(dd[2]), float(dd[3])) for dd in dGStr]
    
    # diff files
    logfile.write("\nDiffing Files ...")
    from math import pi;
    from math import pow;
    delta=1.0
    gamma=1.0e-6
    mMag=4.0*pi*pow(delta,5)/5.0
    ed2Mag=4.0*pi*pow(gamma,2)*pow(delta,5)/75.0
    tol=1.0e-15;
    for i in range(len(d)):
        dn=abs(d[i][0]-dG[i][0])
        dh=abs(d[i][1]-dG[i][1])
        dm=abs(d[i][2]-dG[i][2])
        ded2=abs(d[i][3]-dG[i][3])
        if dn != 0 or dh > tol or dm/mMag > tol or ded2/ed2Mag > tol:
            logfile.write("\nut_bondVolumeConvergenceStudy.dat DIFFERs with ut_bondVolumeConvergenceStudy.gold.dat\n\n")
            sys.exit(1)
            
    logfile.write("\nDiffing Files ...DONE\n")
    sys.exit(0)
