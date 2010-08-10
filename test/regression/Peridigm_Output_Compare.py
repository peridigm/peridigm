#! /usr/bin/env python

import sys
import os
import re
from xml_parser import Tolerance_XML_Parser, PVTU_XML_Parser, VTU_XML_Info_Parser, VTU_XML_Data_Parser
from test_utils import fequal_tol

def IsInvalid(val):
    # other options exist to check for nan, but none is good
    # 1) use numpy.isnan() and numpy.isinf, which obviously requires numpy
    # 2) use math.isnan() and math.isinf(), which are available only in python versions >= 2.6
    if str(val) == str(1e400*0):
        return True
    return False

class VTUData:
    """Container for vtu file data."""

    def __init__(self, *args):
        filename = args[0]
        self.point_data_names = []
        self.pieces = []
        if filename[-5:] == ".pvtu":
            self.ReadPVTUFile(filename)
        else:
            self.ReadVTUFile(filename)
        self.point_data_names.sort()
        self.pieces.sort()

    def ReadPVTUFile(self, filename):
        parser = PVTU_XML_Parser()
        parser.Parse(filename)
        self.pieces =  parser.pieces[:]
        self.point_data_names = parser.point_data_names[:]

    def ReadVTUFile(self, filename):
        parser = VTU_XML_Info_Parser()
        parser.Parse(filename)
        self.pieces = [filename]
        self.point_data_names = parser.point_data_names[:]

    def __repr__(self):
        text = "VTUData\n"
        text += "  point data:\n"
        for name in self.point_data_names:
            text += "    " + str(name) + "\n"
        text += "  pieces:\n"
        for piece in self.pieces:
            text += "    " + str(piece) + "\n"
        return text

def ParseCommandLineArgs(args):

    tolerance_specification = None
    uniform_tolerance = None
    tolerance_file = None
    verbose = False

    index = 0
    while index < len(args):
        item = args[index]
        if item == "-uniform_tolerance":
            if len(args) < index+2:
                print "ERROR, -uniform_tolerance option must be followed by numeric value for tolerance.\n"
                sys.exit(1)
            tolerance_specification = "Uniform"
            uniform_tolerance = float(args[index+1])
            index += 1
        elif item == "-tolerance_file":
            if len(args) < index+2:
                print "ERROR, -tolerance_file option must be followed by name of xml file.\n"
                sys.exit(1)
            tolerance_specification = "From File"
            tolerance_file = args[index+1]
            index += 1
        elif item == "-verbose":
            verbose = True
        else:
            print "Error, unknown option:", item, "\n"
            sys.exit(1)
        index += 1

    return (tolerance_specification, uniform_tolerance, tolerance_file, verbose)

def Print(text, verbose):
    if verbose == True:
        print text
    return

if __name__ == "__main__":

    print "\n----Peridigm Output Comparison----\n"
    if len(sys.argv) < 5:
        print "Usage:  Peridigm_Output_Compare.py <basename1> <basename2> [options]\n"
        print "Options:"
        print "  -uniform_tolerance <tolerance>    compare all fields with specified tolerance."
        print "  -tolerance_file <file.xml>        compare specific fields using tolerances described in input file"
        print "  -verbose                          prints comparison information for all fields at all time steps.\n"
        print "Either uniform tolerance or tolerance file must be specfied.\n"
        sys.exit(1)

    basename_1 = sys.argv[1]
    basename_2 = sys.argv[2]
    tolerance_specification, uniform_tolerance, tolerance_file, verbose = \
        ParseCommandLineArgs(sys.argv[3:])
    print "Base name 1: ", basename_1, "\nBase name 2: ", basename_2, "\n"

    files_diff = False

    # create a list of all the file names
    vtu_files_1 = []
    pvtu_files_1 = []
    vtu_files_2 = []
    pvtu_files_2 = []
    for file in os.listdir(os.getcwd()):
        # behold the power of regular expressions...
        vals = re.split("\.t[0123456789]*\.vtu", file)
        if len(vals) > 1:
            vtu_basename = re.split("\.p[0123456789]*", vals[0])[0]
            if vtu_basename == basename_1:
                vtu_files_1.append(file)
            elif vtu_basename == basename_2:
                vtu_files_2.append(file)
        vals = re.split("\.t[0123456789]*\.pvtu", file)
        if len(vals) > 1:
            pvtu_basename = re.split("\.p[0123456789]*", vals[0])[0]
            if pvtu_basename == basename_1:
                pvtu_files_1.append(file)
            elif pvtu_basename == basename_2:
                pvtu_files_2.append(file)
    vtu_files_1.sort()
    pvtu_files_1.sort()
    vtu_files_2.sort()
    pvtu_files_2.sort()

    # check to make sure data files have been found
    if len(vtu_files_1) == 0 and len(pvtu_files_1) == 0:
        print "ERROR, no data files found for basename:", basename_1, "\n"
        sys.exit(1)
    if len(vtu_files_2) == 0 and len(pvtu_files_2) == 0:
        print "ERROR, no data files found for basename:", basename_2, "\n"
        sys.exit(1)

    # get basic data from the .pvtu file, if one exists
    # otherwise get basic data from the first .vtu file in the series
    num_time_steps = 0
    num_time_steps_2 = 0
    if len(pvtu_files_1) > 0:
        vtu_data_1 = VTUData(pvtu_files_1[0])
        num_time_steps = len(pvtu_files_1)
    else:
        vtu_data_1 = VTUData(vtu_files_1[0])
        num_time_steps = len(vtu_files_1)
    if len(pvtu_files_2) > 0:
        vtu_data_2 = VTUData(pvtu_files_2[0])
        num_time_steps_2 = len(pvtu_files_2)
    else:
        vtu_data_2 = VTUData(vtu_files_2[0])
        num_time_steps_2 = len(vtu_files_2)

    # check to see that the number of times steps are the same
    if num_time_steps != num_time_steps_2:
        print "WARNING, NUMBERS OF TIME STEPS DO NOT MATCH!"
        print "  Data set #1:", num_time_steps, "steps"
        print "  Data set #2:", num_time_steps_2, "steps\n"

    # check to see that data fields are the same
    if vtu_data_1.point_data_names != vtu_data_2.point_data_names:
        print "WARNING, DATA FIELDS DO NOT MATCH!"
        print "  Data set #1:", vtu_data_1.point_data_names
        print "  Data set #2:", vtu_data_2.point_data_names, "\n"

    tolerances = {}
    if tolerance_specification == "Uniform":
        # assign specifiec tolerance to every field
        for item in vtu_data_1.point_data_names:
            if item != "ID":
                tolerances[item] = uniform_tolerance
    else:
        # read tolerances from file
        parser = Tolerance_XML_Parser()
        parser.Parse(tolerance_file)
        tolerances = parser.tolerances
    tolerance_keys = tolerances.keys()
    tolerance_keys.sort()
        
    # print table of field names and tolerances
    print "Data fields to be checked:"
    print "  Data field               Tolerance"
    for key in tolerance_keys:
        text = "   " + key + " "*(25-len(key)) + str(tolerances[key])
        print text
    print

    print "Number of time steps:", num_time_steps, "\n"

    # loop over the time steps and compare the data
    for time_step in range(num_time_steps):
        vtu_data_text = "Time step " + str(time_step) + "\n"
        if len(pvtu_files_1) > 0:
            vtu_data_1 = VTUData(pvtu_files_1[time_step])
        else:
            vtu_data_1 = VTUData(vtu_files_1[time_step])
        if len(pvtu_files_2) > 0:
            vtu_data_2 = VTUData(pvtu_files_2[time_step])
        else:
            vtu_data_2 = VTUData(vtu_files_2[time_step])
        vtu_data_text += "  comparing file set:  " + vtu_data_1.pieces[0]
        for i in range(len(vtu_data_1.pieces)-1):
            vtu_data_text += ", " + vtu_data_1.pieces[i+1]
        vtu_data_text += "\n  against file set:    " + vtu_data_2.pieces[0]
        for i in range(len(vtu_data_2.pieces)-1):
            vtu_data_text += ", " + vtu_data_2.pieces[i+1]
        Print(vtu_data_text, verbose)

        data_1 = {}
        data_2 = {}
        for data_name in tolerance_keys:
            data_1[data_name] = {}
            data_2[data_name] = {}

        # read every file for this time step
        for filename in vtu_data_1.pieces:
            parser = VTU_XML_Data_Parser(tolerance_keys)
            parser.Parse(filename)
            # check for ID field
            if "ID" not in parser.point_data.keys():
                print "ERROR, required ID field not found in", basename_1, "\n"
                sys.exit(1)
            # put the data into the data array
            for key in tolerance_keys:
                for i in range(len(parser.point_data['ID'])):
                    id = parser.point_data['ID'][i]
                    data_1[key][id] = parser.point_data[key][i]
        for filename in vtu_data_2.pieces:
            parser = VTU_XML_Data_Parser(tolerance_keys)
            parser.Parse(filename)
            if "ID" not in parser.point_data.keys():
                print "ERROR, required ID field not found in", basename_2, "\n"
                sys.exit(1)
            for key in tolerance_keys:
                for i in range(len(parser.point_data['ID'])):
                    id = parser.point_data['ID'][i]
                    data_2[key][id] = parser.point_data[key][i]

        # check data lengths
        length_1 = len(data_1[data_1.keys()[0]])
        length_2 = len(data_2[data_2.keys()[0]])
        if length_1 != length_2:
            text = "ERROR, data field length mismatch "
            text += str(length_1) + " != " + str(length_2) + "\n"
            print text
            sys.exit(1)

        # compare the data
        for key in tolerance_keys:
            data_diff = False
            text = "  comparing " + key + "\n"
            text += "    number of values: " + str(len(data_1[key])) + "\n"
            tolerance = tolerances[key]
            text += "    tolerance: " + str(tolerance) + "\n"
            error_text = ""
            max_diff = 0.0
            for id in data_1[key].keys():
                value_1 = data_1[key][id]
                value_2 = data_2[key][id]
                diff = abs(value_1 - value_2)
                if diff > max_diff:
                    max_diff = diff
                if diff > tolerances[key] or IsInvalid(value_1) or IsInvalid(value_2):
                    error_text += "TOLERANCE EXCEEDED node id: " + str(id)
                    error_text += ", value 1: " + str(value_1) + ", value 2: " + str(value_2)
                    error_text +=  ", difference: " + str(diff) + "\n"
                    data_diff = True
                    files_diff = True
            text += "    maximum difference: " + str(max_diff)
            if data_diff == True:
                Print(vtu_data_text, not verbose)
                print text + "\n" + error_text[:-1]
                Print("", not verbose)
            else:
                Print(text, verbose)

        Print("", verbose)

    return_code = 0
    if files_diff == False:
        print "FILES MATCH\n"
    else:
        print "FILES DIFFER\n"
        return_code = 1

    sys.exit(return_code)
