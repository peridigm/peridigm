#! /usr/bin/env python

import sys
import os
import string
from test_utils import fequal_tol, read_line

base_name = "PdQuickGrid_3x3_mpi2"

def parse_file(file):

    data = []
    buff = read_line(file)
    while buff != None:
        vals = string.split(buff)
        for val in vals:
            data.append(val)
        buff = read_line(file)

    return data

if __name__ == "__main__":

    verbose = False
    if "-verbose" in sys.argv:
        verbose = True

    log_file_name = base_name + ".log"

    if os.path.exists(log_file_name):
        os.remove(log_file_name)

    os.system("mpiexec -np 2 ../../src/Peridigm PdQuickGrid_3x3.xml &> " + log_file_name)

    if verbose == True:
        os.system("cat " + log_file_name)

    gold_file = open(base_name+".gold")

    gold_results = parse_file(gold_file)
    gold_file.close()
    
    log_file = open(log_file_name)
    log_results = parse_file(log_file)
    log_file.close()

    if verbose == True:
        print "--Test Results--"
        print "\ngold results values", gold_results
        print "test results values", log_results, "\n"

    result = 0
    if len(gold_results) != len(log_results):
        result = 1
    else:
        for i in range(len(gold_results)-1):
            if log_results[i] != gold_results[i]:
                result = 1
        # final value is run time
#        if fequal_tol(float(gold_results[-1]), float(log_results[-1]), 1.0) == False:
#            result = 1

    sys.exit(result)
