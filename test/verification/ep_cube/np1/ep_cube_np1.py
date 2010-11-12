#! /usr/bin/env python

import sys
import os
import re
from subprocess import Popen

test_dir = "ep_cube/np1"
base_name = "ep_cube"

if __name__ == "__main__":

    result = 0

    # log file will be dumped if verbose option is given
    verbose = False
    if "-verbose" in sys.argv:
        verbose = True

    # change to the specified test directory
    os.chdir(test_dir)

    # open log file
    log_file_name = base_name + ".log"
    if os.path.exists(log_file_name):
        os.remove(log_file_name)
    logfile = open(log_file_name, 'w')

    # remove old output files, if any
    # use regular expression module since differentiating
    # between gold files and old output files can be tricky
    files_to_remove = []
    for file in os.listdir(os.getcwd()):
        vals = re.split("\_t[0123456789]*\.vtu", file)
        if len(vals) > 1:
            vtu_basename = re.split("\.p[0123456789]*", vals[0])[0]
            if vtu_basename == base_name:
                files_to_remove.append(file)
        vals = re.split("\_t[0123456789]*\.pvtu", file)
        if len(vals) > 1:
            pvtu_basename = re.split("\.p[0123456789]*", vals[0])[0]
            if pvtu_basename == base_name:
                files_to_remove.append(file)
    for file in files_to_remove:
        os.remove(file)

    # run Peridigm
    command = ["../../../../src/Peridigm", "../"+base_name+".xml"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # compare output files against gold files
    command = ["../../../../scripts/vtkDiff.py", \
                   base_name+".pvd", \
                   base_name+"_gold.pvd", \
                   "--toleranceFile", \
                   "../"+base_name+".comp"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    logfile.close()

    # dump the output if the user requested verbose
    if verbose == True:
        os.system("cat " + log_file_name)

    sys.exit(result)
