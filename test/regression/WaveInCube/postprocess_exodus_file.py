#!/usr/bin/env python

# This script requires exodus.py
# It can be tricky to get exodus.py to work because it requires libnetcdf.so and libexodus.so
# Here's one approach that worked on the CEE LAN:
# 1) Append sys.path, as shown below, to include the bin subdirectory of your Trilinos install directory.
#    The file exodus.py is in this directory.
# 2) Edit exodus.py as follows (approximately line 71):
#    accessPth = "/projects/seacas/linux_rhel6/current"
#    The path above is valid on the CEE LAN.  On other systems, you need to provide a path to a SEACAS build
#    that includes shared libraries.
import sys
sys.path.append('/Users/djlittl/Software/seacas/GCC_4.9.4/lib')
import exodus

import string

if __name__ == "__main__":

    inFileName = "WaveInCube.e"
    inFile = exodus.exodus(inFileName, mode='r')

    outFileLabel = string.splitfields(inFileName, '.')[0]

    # Print database parameters from inFile
    print " "
    print "Database version:         " + str(round(inFile.version.value,2))
    print "Database title:           " + inFile.title()
    print "Database dimensions:      " + str(inFile.num_dimensions())
    print "Number of nodes:          " + str(inFile.num_nodes())
    print "Number of elements:       " + str(inFile.num_elems())
    print "Number of element blocks: " + str(inFile.num_blks())
    print "Number of node sets:      " + str(inFile.num_node_sets())
    print "Number of side sets:      " + str(inFile.num_side_sets())
    print " "

    # Extract nodal displacements and forces
    numNodes = inFile.num_nodes()
    numTimeSteps = inFile.num_times()
    nodeVariableNames = inFile.get_node_variable_names()
    coords = inFile.get_coords()
    num_nodes = len(coords[0])
    print nodeVariableNames
    if 'DisplacementX' not in nodeVariableNames:
        print "\nERROR:  Failed to extract DisplacementX data\n"
        sys.exit(1)

    print "\nProcessing", num_nodes, "nodes...\n"

    target_y = 0.01
    target_z = 0.01
    tol = 1.0e-3

    initial_data = []

    timeStep = 1
    displacement_x = inFile.get_node_variable_values('DisplacementX', timeStep)
    for i in range(num_nodes):
        x = coords[0][i]
        y = coords[1][i]
        z = coords[2][i]
        disp = displacement_x[i]

        if abs(y - target_y) < tol and abs(z - target_z) < tol:
            initial_data.append([x, y, z, disp])

    final_data = []

    timeStep = numTimeSteps
    displacement_x = inFile.get_node_variable_values('DisplacementX', timeStep)
    for i in range(num_nodes):
        x = coords[0][i]
        y = coords[1][i]
        z = coords[2][i]
        disp = displacement_x[i]

        if abs(y - target_y) < tol and abs(z - target_z) < tol:
            final_data.append([x, y, z, disp])

    inFile.close()

    initial_data.sort()
    final_data.sort()

    outFile = open(outFileLabel + ".txt", 'w')
    for i in range(len(initial_data)):
        outFile.write(str(initial_data[i][0]) + " " + str(initial_data[i][3]) + " " + str(final_data[i][3]) + "\n")
    outFile.close()

    print "\nData written to " + outFileLabel + ".txt"
    print
