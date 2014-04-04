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

    # MLP: Blow away Dakota log file here
    if os.path.exists("dakota_output.txt"):
        os.remove("dakota_output.txt")
    dakota_output = open("dakota_output.txt", 'w')

    # remove old output files, if any
    files_to_remove = base_name + ".e"
    for file in os.listdir(os.getcwd()):
      if file in files_to_remove:
        os.remove(file)

    # Setup soft link to Peridigm executable so that simulation driver script can find it
    command = ["ln","-f","-s","../../../../src/Peridigm","."]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    command = ["ln","-f","-s","../../../../../src/Peridigm","./templatedir"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code
    command = ["ln","-f","-s","../../../../../scripts/exotxt","./templatedir"]
    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # run Dakota
    sys.path.append("../../../../src/Peridigm")
    command = ["dakota","-in","dakota_peridigm.in"]
    p = Popen(command, stdout=dakota_output, stderr=dakota_output)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # Create results file to parse
    command = ["grep -A 1 -i \"<<<<< Best parameters\" dakota_output.txt | tail -n 1 | sed 's/^[ \t]*//' | cut -d\" \" -f1 > Dakota_2x1x1.dat"]
    p = Popen(command, stdout=logfile, stderr=logfile, shell=True)
    return_code = p.wait()
    if return_code != 0:
        result = return_code

    # compare output files against gold files
    filename1 = "../"+base_name+"_gold.dat"
    filename2 = base_name+".dat"
    f1 = open(filename1, 'r')
    f2 = open(filename2, 'r')
    str1 = f1.readline()
    str2 = f2.readline()
    f1.close()
    f2.close()
    val1 = float(str1)
    val2 = float(str2)
    diff = abs((val1 - val2) / val1)
    logfile.write("diff = "+str(diff))
    if diff > 1e-8:
       result = -1
    #command = ["diff "+base_name+".dat "+"../"+base_name+"_gold.dat"]
    #p = Popen(command, stdout=logfile, stderr=logfile, shell=True)
    #return_code = p.wait()
    #if return_code != 0:
    #    result = return_code

    dakota_output.close()
    logfile.close()

    # dump the output if the user requested verbose
    if verbose == True:
        os.system("cat " + log_file_name)

    sys.exit(result)
