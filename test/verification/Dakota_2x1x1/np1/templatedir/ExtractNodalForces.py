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
    flag = False
    for i in range(len(lines)):
        line = lines[i]
        if flag:
            if "!" in line:
                flag = False
            else:
                data.append(line)
        if "Nodal variables" in line:
            data = []
            flag = True

    values = []
    for datum in data:
        vals = string.splitfields(datum)
        for val in vals:
            values.append(val)

    num_nodal_variables = 9
    pt1_vals = values[:num_nodal_variables]
    pt2_vals = values[num_nodal_variables:2*num_nodal_variables]

    force_index = 6
    pt1_force = float(pt1_vals[force_index])
    pt2_force = float(pt2_vals[force_index])

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
