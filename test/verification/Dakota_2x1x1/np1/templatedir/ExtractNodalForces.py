#! /usr/bin/env python

import sys
import os
import re
import string
from subprocess import Popen

if __name__ == "__main__":

    logfile = open("exotxt.log", 'w')
    command = ["exotxt","peridigm.e","peridigm.txt"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    logfile.close()

    text_file = open("peridigm.txt")
    lines = text_file.readlines()
    text_file.close()

    data = []
    names = []
    namesFlag = False
    dataFlag = False
    for i in range(len(lines)):
        line = lines[i]
        if namesFlag:
            if line[0] == "!":
                namesFlag = False
            else:
                names.append(line)
        if dataFlag:
            if "!" in line:
                dataFlag = False
            else:
                data.append(line)
        if "Variable names" in line:
            names = []
            namesFlag = True
        if "Nodal variables" in line:
            data = []
            dataFlag = True

    var_names = []
    for datum in names:
        vals = string.splitfields(datum)
        for val in vals:
            var_names.append(val)

    num_global_vars = int(var_names[0])
    num_nodal_vars = int(var_names[1])
    num_element_vars = int(var_names[2])

    # get rid of comments read in from exotxt file
    var_names = var_names[8:]

    # find the location of the quantity of interest
    var_of_interest = "FORCE_DENSITYX"
    index = 0
    for i in xrange(len(var_names)):
        if var_names[i] == var_of_interest:
            index = i
            print "\nExtractNodalForces.py found variable", var_of_interest, "at location", index

    values = []
    for datum in data:
        vals = string.splitfields(datum)
        for val in vals:
            values.append(val)

    pt1_vals = values[:num_nodal_vars]
    pt2_vals = values[num_nodal_vars:2*num_nodal_vars]

    pt1_force = float(pt1_vals[index])
    pt2_force = float(pt2_vals[index])

    pt1_truth = -11700000.0
    pt2_truth =  11700000.0

    error1 = pt1_force - pt1_truth
    error2 = pt2_force - pt2_truth

    if len(sys.argv) != 2:
       print "\nError:  No output file specified."
       print "Usage:  ExtractNodalForces.py <output_file>.\n"
       sys.exit(-1)
    results_string = str(error1) + "\n" + str(error2) + "\n"
    dakota_file_name = sys.argv[-1]
    dakota_file = open(dakota_file_name, 'w')
    dakota_file.write(results_string)
    dakota_file.close()
    print "\nResults extracted by ExtractNodalForces.py:\n"
    print results_string
    print "Results written to", dakota_file_name, "\n"
