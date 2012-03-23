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
    files_to_join = ["Contact_2x1x1-s00001", "Contact_2x1x1-s00002", "Contact_2x1x1-s00003", "Contact_2x1x1-s00004",
                     "Contact_2x1x1-s00005", "Contact_2x1x1-s00006", "Contact_2x1x1-s00007", "Contact_2x1x1-s00008",
                     "Contact_2x1x1-s00009", "Contact_2x1x1-s00010", "Contact_2x1x1-s00011"]
    for file in files_to_join:
      command = ["../../../../scripts/epu", "-p", "2", file]
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
          result = return_code

    # Now combine time series from all databaases
    command = ["../../../../scripts/conjoin", "-output", base_name+".e", 
               "Contact_2x1x1-s00001.e", "Contact_2x1x1-s00002.e", "Contact_2x1x1-s00003.e", "Contact_2x1x1-s00004.e",
               "Contact_2x1x1-s00005.e", "Contact_2x1x1-s00006.e", "Contact_2x1x1-s00007.e", "Contact_2x1x1-s00008.e",
               "Contact_2x1x1-s00009.e", "Contact_2x1x1-s00010.e", "Contact_2x1x1-s00011.e"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    
    # Now merged output file against gold file
    command = ["../../../../scripts/exodiff", \
                   base_name+".e", \
                   "../"+base_name+"_gold.e", \
                   "-f", \
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
