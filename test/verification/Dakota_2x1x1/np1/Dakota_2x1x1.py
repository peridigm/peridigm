#! /usr/bin/env python

import sys
import os
import re
import string
from subprocess import Popen

test_dir = "Dakota_2x1x1/np1"
base_name = "Dakota_2x1x1"

if __name__ == "__main__":

    # This test used Dakota to calibrate a constituive model
    # via optimization.  The know correct value for the bulk
    # moduls is 130.0 GPA
    MODULUS_TRUTH_VALUE = 1.3000000000e+11
    REL_ERROR_TOLERANCE = 1.0e-8

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

    # open log file
    dakota_log_file_name = "Dakota.log"
    if os.path.exists(dakota_log_file_name):
        os.remove(dakota_log_file_name)

    # remove old output files, if any
    files_to_remove = base_name + ".e"
    for file in os.listdir(os.getcwd()):
      if file in files_to_remove:
        os.remove(file)

    # Setup soft link to Peridigm executable so that simulation driver script can find it
    logfile.write("Creating symbolic link to Peridigm in test directory.\n")
    command = ["ln","-f","-s","../../../../src/Peridigm","."]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    logfile.write("  Complete (return code " + str(return_code) + ").\n\n")
    logfile.write("Creating symbolic link to Peridigm in template directory.\n")
    command = ["ln","-f","-s","../../../../../src/Peridigm","./templatedir"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    logfile.write("  Complete (return code " + str(return_code) + ").\n\n")
    logfile.write("Creating symbolic link to exotxt in template directory.\n")
    command = ["ln","-f","-s","../../../../../scripts/exotxt","./templatedir"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    logfile.write("  Complete (return code " + str(return_code) + ").\n\n")

    # run Dakota
    logfile.write("Executing Dakota.\n")
    sys.path.append("../../../../src/Peridigm")
    command = ["dakota","-in","dakota_peridigm.in"]
    dakota_log_file = open(dakota_log_file_name, 'w')
    p = Popen(command, stdout=dakota_log_file, stderr=dakota_log_file)
    return_code = p.wait()
    dakota_log_file.close()
    if return_code != 0:
        result = return_code
    logfile.write("  Complete (return code " + str(return_code) + ").\n\n")

    # extract the modulus that Dakota solved for from the log file
    computed_bulk_modulus = 0.0
    dakota_log_file = open(dakota_log_file_name)
    lines = dakota_log_file.readlines()
    dakota_log_file.close()
    for i in range(len(lines)):
        line = lines[i]
        if "Best parameters" in line:
            vals = string.splitfields(line + lines[i+1])
            for i in range(len(vals)):
                if "BulkModulus" in vals[i]:
                    computed_bulk_modulus = float(vals[i-1])

    logfile.write("Computed bulk modulus: " + str(computed_bulk_modulus) + "\n")
    logfile.write("Truth value: " + str(MODULUS_TRUTH_VALUE) + "\n")
    rel_error = abs((MODULUS_TRUTH_VALUE - computed_bulk_modulus)/MODULUS_TRUTH_VALUE)
    logfile.write("Relative error: " + str(rel_error) + "\n")
    logfile.write("Tolerance: " + str(REL_ERROR_TOLERANCE) + "\n")

    if rel_error > REL_ERROR_TOLERANCE:
       result = -1

    logfile.write("\nTest return code: " + str(result) + "\n\n")
    logfile.close()

    # dump the output if the user requested verbose
    if verbose == True:
        os.system("cat " + log_file_name)

    sys.exit(result)
