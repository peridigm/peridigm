#! /usr/bin/env python

import os
import string
import time

def parse_log_files():

    # find all the log file (may pick up some files unrelated to performance testing)
    log_files = []
    for root, dirs, files in os.walk('./test/performance'):
        for log_file in files:
            if log_file[-4:] == ".log":
                file_with_path = os.path.join(root, log_file) 
                log_files.append(file_with_path)

    results = []

    for log_file in log_files:
        is_valid = False
        test_name = log_file[string.rfind(log_file, '/')+1 : -4]
        wallclock_time = 0.0
        benchmark_value = 0.0
        tolerance = 0.0
        fin = open(log_file, 'r')
        lines = fin.readlines()
        fin.close()
        for line in lines:
            if "wallclock time  =" in line:
                wallclock_time = float( line[string.rfind(line, '= ')+1 : ] )
                is_valid = True
            elif "benchmark value =" in line:
                benchmark_value = float( line[string.rfind(line, '= ')+1 : ] )
            elif "tolerance       =" in line:
                tolerance = float( line[string.rfind(line, '= ')+1 : ] )

        if is_valid:
            results.append([test_name, wallclock_time, benchmark_value, tolerance])

    return results

if __name__ == "__main__":

    results = parse_log_files()

    print "\n------------------------------------------------------------------------------"
    print "\nDate and time:", time.strftime("%m-%d-%Y %H:%M")
    print "\nTest Name                   Wallclock Time     Benchmark Value     Tolerance"
    for result in results:
        print '{0:25}  {1:6}  {2:17}  {3:18}'.format(result[0], result[1], result[2], result[3])
    print
