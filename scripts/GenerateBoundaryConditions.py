#! /usr/bin/env python

import sys
sys.path.append("../test/regression")
from math import pi, sin, cos, atan2
from xml_parser import VTU_XML_Points_Parser
import random

# DEFINE BOUNDARY CONDITION FUNCTIONS HERE

# Each boundary condition needs a global ID (single int), 
# a coordinate ("x", "y", or "z"), and a magnitude (single float)
# boundary conditions are stored as a list of length three,
# for example [5, "x", 1.2]

def Q2CylinderBoundaryConditions(points):

    # find height of cylinder
    min_z = 1.0e10
    max_z = -1.0e10
    for pt in points:
        if pt[2] > max_z:
            max_z = pt[2]
        if pt[2] < min_z:
            min_z = pt[2]
    height = max_z - min_z

    # parameters for determining initial conditions
    v_r0 = 200.0
    v_r1 = 50.0
    v_z0 = 100.0
    a = height/2.0

    # find the velocities at each point and create a 
    # corresponding boundary condition
    bcs = []
    vel = [0.0]*3
    for i in xrange(len(points)):
        pt = points[i]
        z = pt[2] - min_z - a
        v_r = v_r0 - v_r1*(z/a)*(z/a)
        v_z = v_z0*(z/a)
        v_theta = 0.0
        theta = atan2(pt[1], pt[0])
        vel[0] = v_r*cos(theta) - v_theta*sin(theta)
        vel[1] = v_r*sin(theta) + v_theta*cos(theta)
        vel[2] = v_z
        bcs.append([i, "x", vel[0]])
        bcs.append([i, "y", vel[1]])
        bcs.append([i, "z", vel[2]])

    return bcs

# END OF BOUNDARY CONDITION FUNCTIONS

def PerturbInitialConditions(bcs, magnitude):

    for bc in bcs:
        bc[2] += (random.random()*2.0 - 1.0)*magnitude

    return

def PrintBoundaryConditions(bcs):
    
    node_sets = []
    boundary_conditions = []

    for bc in bcs:
        # create a node set
        node_set = "      <Parameter name=\"Node Set " + str(bc[0]) + "\" type=\"string\" value=\"" + str(bc[0]) + "\"/>"
        node_sets.append(node_set)
        # create the boundary condition
        boundary_condition =  "      <ParameterList name=\"Initial Velocity " + str(bc[0]) + " " + bc[1] + "\">\n"
        boundary_condition += "        <Parameter name=\"Type\" type=\"string\" value=\"Initial Velocity\"/>\n"
	boundary_condition += "        <Parameter name=\"Node Set\" type=\"string\" value=\"Node Set " + str(bc[0]) + "\"/>\n"
        boundary_condition += "        <Parameter name=\"Coordinate\" type=\"string\" value=\"" + bc[1] + "\"/>\n"
        boundary_condition += "        <Parameter name=\"Value\" type=\"double\" value=\"" + str(bc[2]) + "\"/>\n"
	boundary_condition += "      </ParameterList>"
        boundary_conditions.append(boundary_condition)

    # remove duplicates from node_set in a rather cryptic but fast way
    # node_sets will end up unsorted
    keys = {}
    for item in node_sets:
        keys[item] = 1
    node_sets = keys.keys()

    print "    <ParameterList name=\"Boundary Conditions\">"
    for node_set in node_sets:
        print node_set
    for boundary_condition in boundary_conditions:
        print boundary_condition
    print "    </ParameterList>"

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print "\n----Peridigm Boundary Condition Generator----\n"
        print "Usage:  PeridigmBoundaryConditions.py <inputfile.xml>\n"
        sys.exit(1)

    infile_name = sys.argv[1]

    # parse the xml output file and get a list of points
    points_parser = VTU_XML_Points_Parser()
    points_parser.Parse(infile_name)
    pts = points_parser.points

    # organize the points data into tuples
    points = []
    for i in xrange(len(pts)/3):
        points.append((pts[i*3], pts[i*3+1], pts[i*3+2]))

    #
    # Now we have the points as a list of (x, y, z) tuples
    # Pass these data to a function that computes boundary conditions
    # Each boundary condition needs a global ID (single int), 
    # a coordinate ("x", "y", or "z"), and a magnitude (single float)
    # boundary conditions are stored as a tuple, for example (5, "x", 1.2)
    #

    # call a function to define the boundary conditions
    bcs = Q2CylinderBoundaryConditions(points)
    
    # perturb the initial conditions
    magnitude = 5.0
    PerturbInitialConditions(bcs, magnitude)

    # print out the boundary conditions
    PrintBoundaryConditions(bcs)
