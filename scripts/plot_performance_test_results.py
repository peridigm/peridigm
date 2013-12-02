#! /usr/bin/env python

import sys
import os
import string

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "\nUsage:  process_accumulated_performance_test_results.py <in_file.txt>\n"
        sys.exit(1)

    infile_name = sys.argv[-1]

    infile = open(infile_name, 'r')
    lines = infile.readlines()
    infile.close()

    results = {}

    date = ""
    for line in lines:
        wallclock_time = 0.0
        benchmark_value = 0.0
        tolerance = 0.0
        test_result = "Failed"
        if "Date and time:" in line:
            date = string.splitfields(line)[-2]
        elif "Test Name" not in line and "-----------" not in line:
            vals = string.splitfields(line)
            if len(vals) == 5:
                name = vals[0]
                wallclock_time = float(vals[1])
                benchmark_value = float(vals[2])
                tolerance = float(vals[3])
                test_result = vals[4]
                if name not in results.keys():
                    results[name] = []
                results[name].append([date, wallclock_time, benchmark_value, tolerance])

    data_file_names = []
    for key in results.keys():
        fout_name = key + "_accumulated_performance_data.txt"
        fout = open(fout_name, 'w')
        for result in results[key]:
            fout.write(result[0] + " " + str(result[1]) + " " + str(result[2]) + " " + str(result[3]) + "\n")
        fout.close()
        data_file_names.append(fout_name)

    gnuplot_file_name = "make_accumulated_performance_plots.plotgen"
    pdf_file_names = []
    gfile = open(gnuplot_file_name, 'w')
    for key in results.keys():
        result = results[key]
        title = key.replace('_', '\_') + " Performance Test Data"
        pdf_file_name = key + "_performance.pdf"
        data_file_name = key + "_accumulated_performance_data.txt"
        yrange_min = result[-1][2] - 2.0*result[-1][3]
        yrange_max = result[-1][2] + 2.0*result[-1][3]
        gfile.write("\nset terminal pdf enhanced font \"Times-Roman,24\" size 12in, 8in\n")
        gfile.write("set title \'" + title + "\'\n")
        gfile.write("set output \"" + pdf_file_name + "\"\n")
        gfile.write("set xlabel \"Date\" font \"Times-Roman,32\"\n")
        gfile.write("set xtics rotate by 270 font \"Times-Roman,14\"\n")
        gfile.write("set ylabel \"Run Time (sec.)\" font \"Times-Roman,32\"\n")
        gfile.write("set yrange [" + str(yrange_min) + ":" + str(yrange_max) + "]\n")
        gfile.write("plot \"" + data_file_name + "\" using 2:xticlabel(1) with points pt 7 ps 4 notitle, \\\n")
        gfile.write("     \"" + data_file_name + "\" using 3:xticlabel(1) with lines lt 2 lw 1 notitle, \\\n")
        gfile.write("     \"" + data_file_name + "\" using ($3+$4):xticlabel(1) with lines lt 0 lw 1 notitle, \\\n")
        gfile.write("     \"" + data_file_name + "\" using ($3-$4):xticlabel(1) with lines lt 0 lw 1 notitle\n")
        pdf_file_names.append(pdf_file_name)
    gfile.close()

    print "\nGenerated the following files:"
    for name in data_file_names:
        print "  ", name
    print "\nGnuplot command file written to: ", gnuplot_file_name
    print "\nCommand to generate plots:"
    print "  gnuplot", gnuplot_file_name
    print "\nCommand to concatenate pdf files:"
    cmd_string = "  gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dAutoRotatePages -sOutputFile=performace_plots.pdf"
    for name in pdf_file_names:
        cmd_string += " " + name
    print cmd_string, "\n"
