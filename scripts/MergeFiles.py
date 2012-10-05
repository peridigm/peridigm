#! /usr/bin/env python

import sys
import os
import re
import glob
from subprocess import Popen

def sort_nicely( l ):
  # Sort the given list in the way that humans expect.
  convert = lambda text: int(text) if text.isdigit() else text
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
  l.sort( key=alphanum_key )

if __name__ == "__main__":

    if len(sys.argv) < 2:
      print "\n----------------Peridigm Exodus Outfile Merger----------------\n"
      print "Usage:  MergeFiles.py   <File Base Name>    <# of Processors>\n"
      sys.exit(1)

    path = sys.argv[0];
    base_name = sys.argv[1]
    num_proc = sys.argv[2]

    result = 0
    
    # Define path to location of MergeFiles to call epu and conjoin later
    start = 0
    end   = path.find('MergeFiles')
    newfile = path[start:end]
    path = newfile

    # remove the output .e file if it exists
    if os.path.exists(base_name+".e"):
      os.remove(base_name+".e")

    # open log file
    log_file_name = base_name + ".log"
    if os.path.exists(log_file_name):
      os.remove(log_file_name)
    logfile = open(log_file_name, 'w')

    # Generate listing of files <basename>-s* for merging
    files = glob.glob('*.e*')
    files_to_join = []
    for file in files:
      start = 0
      end   = file.find('.e')
      newfile = file[start:end];
      files_to_join.append(newfile)

    # Remove duplicates
    files_to_join = list(set(files_to_join))
    files_to_join.sort()

    # First merge all distributed exodus databases for each time stamp
    for file in files_to_join:
      command = ["epu", "-p", num_proc, file]
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
          result = return_code

    # Check for any "-s" files
    check = glob.glob('*-s*')
    if len(check) != 0:
      # Add .e to list of files for conjoin
      files_to_conjoin = []
      for file in files_to_join:
        newfile = file + ".e"
        files_to_conjoin.append(newfile)
      sort_nicely(files_to_conjoin)
      # Now combine time series from all databases
      command = ["conjoin", "-output", base_name+".e"]
      for file in files_to_conjoin:
        command.append(file)
      p = Popen(command, stdout=logfile, stderr=logfile)
      return_code = p.wait()
      if return_code != 0:
        result = return_code
    
      # If error print message to screen
      if result != 0:
        print "\n Error when merging files! Check .log file for error details.\n"

    sys.exit(result)
