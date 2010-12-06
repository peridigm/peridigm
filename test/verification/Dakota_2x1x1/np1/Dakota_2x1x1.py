#! /usr/bin/env python

import sys
import os
import re
from subprocess import Popen

test_dir = "Dakota_2x1x1/np1"
base_name = "Dakota_2x1x1"

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

    # Setup soft link to Peridigm executable so that simulation driver script can find it
    command = ["ln -f -s ../../../../src/Peridigm ."]
    p = Popen(command, stdout=logfile, stderr=logfile, shell=True)
    return_code = p.wait()
    print return_code
    if return_code != 0:
        result = return_code

    # run Dakota
    command = ["dakota -in dakota_peridigm.in > dakota_output.txt"]
    p = Popen(command, stdout=logfile, stderr=logfile, shell=True)
    return_code = p.wait()
    print return_code
    if return_code != 0:
        result = return_code

    # Create results file to parse
    command = ["grep -A 1 -i \"<<<<< Best parameters\" dakota_output.txt | tail -n 1 | sed 's/^[ \t]*//' | cut -d\" \" -f1 > Dakota_2x1x1.dat"]
    p = Popen(command, stdout=logfile, stderr=logfile, shell=True)
    print return_code
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # compare output files against gold files
    command = ["diff "+base_name+".dat "+"../"+base_name+"_gold.dat"]
    print command
    p = Popen(command, stdout=logfile, stderr=logfile, shell=True)
    return_code = p.wait()
    print return_code
    if return_code != 0:
        result = return_code

    logfile.close()

    # dump the output if the user requested verbose
    if verbose == True:
        os.system("cat " + log_file_name)

    sys.exit(result)
