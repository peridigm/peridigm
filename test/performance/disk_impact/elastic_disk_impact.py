#! /usr/bin/env python

import sys
import os
import re
import glob
from subprocess import Popen

test_dir = "elastic_disk_impact"
base_name = "elastic_disk_impact"

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
    files_to_remove.append('*.out')
    files_to_remove.append('*.nem')
    files_to_remove.append('*.pex')
    print "FILES TO REMOVE", files_to_remove
    for file in os.listdir(os.getcwd()):
      if file in files_to_remove:
        os.remove(file)

    # decompose the mesh file
    command = ["../../scripts/decomp", "-p", "4", base_name]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # run Peridigm
    command = ["mpiexec", "-np", "4", "../../../../src/Peridigm", "../"+base_name+".xml"]    
    p = Popen(command, stdout=subprocess.PIPE)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    out, err = p.communicate()
    print "OUT\n", out, "\n"
    print "ERR\n", err, "\n"

    # concatenate output files
    command = ["../../../../scripts/epu", "-p", "4", base_name]
    p = Popen(command, shell=True, stderr=subprocess.PIPE)
    while True:
        out = p.stderr.read(1)
        if out == '' and p.poll() != None:
            break
        if out != '':
            logfile.write
        
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # compare performance statistics against gold statistics
    perf_file_name = base_name + ".perf"





    # compare output against gold file only if the gold file is present
    gold_file_name = base_name + "_gold.e"
    if os.path.exists(gold_file_name):
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
    else:
        logfile.write("\nGold file " + gold_file_name + " not found, skipping exodiff.\n")
            
    


    logfile.close()

    # dump the output if the user requested verbose
    if verbose == True:
        os.system("cat " + log_file_name)

    sys.exit(result)
