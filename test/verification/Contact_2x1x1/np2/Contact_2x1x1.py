#! /usr/bin/env python

import sys
import os
import re
import glob
from subprocess import Popen

test_dir = "Contact_2x1x1/np2"
base_name = "Contact_2x1x1"

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

    # run Peridigm
    command = ["mpiexec", "-np", "2", "../../../../src/Peridigm", "../"+base_name+".xml"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # First merge all distributed exodus databases for each time stamp
    files_to_join = ["Contact_2x1x1-s1", "Contact_2x1x1-s2", "Contact_2x1x1-s3", "Contact_2x1x1-s4",
                     "Contact_2x1x1-s5", "Contact_2x1x1-s6", "Contact_2x1x1-s7", "Contact_2x1x1-s8",
                     "Contact_2x1x1-s9", "Contact_2x1x1-s10", "Contact_2x1x1-s11"]
    for file in files_to_join:
      command = ["../../../../scripts/epu", "-p", "2", file]
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
          result = return_code

    # Now combine time series from all databaases
    command = ["../../../../scripts/conjoin", "-output", base_name+".e", 
               "Contact_2x1x1-s1.e", "Contact_2x1x1-s2.e", "Contact_2x1x1-s3.e", "Contact_2x1x1-s4.e",
               "Contact_2x1x1-s5.e", "Contact_2x1x1-s6.e", "Contact_2x1x1-s7.e", "Contact_2x1x1-s8.e",
               "Contact_2x1x1-s9.e", "Contact_2x1x1-s10.e", "Contact_2x1x1-s11.e"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    
    # Now merged output file against gold file
    command = ["../../../../scripts/exodiff", \
                   "-stat", \
                   "-f", \
                   "../"+base_name+".comp", \
                   base_name+".e", \
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
