#! /usr/bin/env python

import os
import sys
import time

def parse_log_files():

    # find all the log file (may pick up some files unrelated to performance testing)
    log_files = []
    for root, dirs, files in os.walk('./'):
        for log_file in files:
            if log_file[-4:] == ".log":
                file_with_path = os.path.join(root, log_file) 
                log_files.append(file_with_path)

    results = []

    for log_file in log_files:
        is_valid = False
        test_name = log_file[log_file.rfind('/')+1 : -4]
        wallclock_time = 0.0
        benchmark_value = 0.0
        tolerance = 0.0
        test_result = "Failed"
        fin = open(log_file, 'r')
        lines = fin.readlines()
        fin.close()
        for line in lines:
            if "wallclock time  =" in line:
                wallclock_time = float( line[line.rfind('= ')+1 : ] )
                is_valid = True
            elif "benchmark value =" in line:
                benchmark_value = float( line[line.rfind('= ')+1 : ] )
            elif "tolerance       =" in line:
                tolerance = float( line[line.rfind('= ')+1 : ] )
            elif "PERFORMANCE TEST PASSED" in line:
                test_result = "Passed"
        if is_valid:
            results.append([test_name, wallclock_time, benchmark_value, tolerance, test_result])

    return results

def summarize_performance_tests():

    results = parse_log_files()

    table_str = ""
    table_str += "\n---------------------------------------------------------------------------------------------\n"
    table_str += "\nDate and time:  " + time.strftime("%m-%d-%Y %H:%M") + "\n"
    table_str += "\nTest Name                   Wallclock Time     Benchmark Value     Tolerance     Test Result\n"
    for result in results:
        table_str += '{0:25}  {1:6}  {2:17}  {3:17}          {4:6}\n'.format(result[0], result[1], result[2], result[3], result[4])
    table_str += "\n"

    return table_str

if __name__ == "__main__":

    log_file_name = "./test_summary/performance_test_summary.txt"
    if os.path.exists(log_file_name):
        os.remove(log_file_name)

    summary_table = summarize_performance_tests()
    log_file = open(log_file_name, 'w')
    log_file.write(summary_table)
    log_file.close()

    sys.exit(0)
