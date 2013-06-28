#! /usr/bin/env python

import sys
import os
import re
import glob
from subprocess import Popen

test_dir = "Contact_Perforation/np3"
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

    # run Peridigm
    command = ["mpiexec", "-np", "3", "../../../../src/Peridigm", "../"+base_name+".xml"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # First merge all distributed exodus databases for each time stamp
    files_to_join = ["Contact_Perforation-s1", "Contact_Perforation-s2", "Contact_Perforation-s3", "Contact_Perforation-s4",
                     "Contact_Perforation-s5", "Contact_Perforation-s6", "Contact_Perforation-s7", "Contact_Perforation-s8",
                     "Contact_Perforation-s9", "Contact_Perforation-s10", "Contact_Perforation-s11"]
    for file in files_to_join:
      command = ["../../../../scripts/epu", "-p", "3", file]
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
          result = return_code

    # Now combine time series from all databaases
    command = ["../../../../scripts/conjoin", "-output", base_name+".e", 
               "Contact_Perforation-s1.e", "Contact_Perforation-s2.e", "Contact_Perforation-s3.e", "Contact_Perforation-s4.e",
               "Contact_Perforation-s5.e", "Contact_Perforation-s6.e", "Contact_Perforation-s7.e", "Contact_Perforation-s8.e",
               "Contact_Perforation-s9.e", "Contact_Perforation-s10.e", "Contact_Perforation-s11.e"]

    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # compare output files against gold files
#    command = ["../../../../scripts/epu", "-p", "2", base_name]
#    p = Popen(command, stdout=logfile, stderr=logfile)
#    return_code = p.wait()
#    if return_code != 0:
#        result = return_code
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
