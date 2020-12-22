#! /usr/bin/env python

import sys
import os
import glob
import shutil
from subprocess import Popen

test_dir = "Contact_Perforation_With_Restart/np1"
base_name = "Contact_Perforation"

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
    files_to_remove = glob.glob('*.e*')
    for file in os.listdir(os.getcwd()):
      if file in files_to_remove:
        os.remove(file)
    # remove old restart folders, if any
    path = os.getcwd()
    pattern = os.path.join(path, "restart-*")
    for item in glob.glob(pattern):
      if not os.path.isdir(item):
        continue
      shutil.rmtree(item)

    # run Peridigm
    command = ["../../../../src/Peridigm", "../"+base_name+"_Run1.xml"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    command = ["../../../../src/Peridigm", "../"+base_name+"_Run2.xml"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code	

    # compare output files against gold files
    command = ["../../../../scripts/exodiff", \
               "-stat", \
               "-steps",\
               "last",\
               "-f", \
               "../"+base_name+".comp", \
               base_name+"_Run2.e", \
               "../"+base_name+"_gold.e"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    logfile.close()

    # dump the output if the user requested verbose
    if verbose == True:
        os.system("cat " + log_file_name)

    sys.exit(result)
