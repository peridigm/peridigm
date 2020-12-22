#!/usr/bin/env python

""" peridigm_to_yaml.py:  Converts a *.peridigm input file to the *.yaml format." """

__author__ = "David Littlewood (djlittl@sandia.gov)"

# ************************************************************************
#
#
#                             Peridigm
#                 Copyright (2011) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions?
# David J. Littlewood   djlittl@sandia.gov
# John A. Mitchell      jamitch@sandia.gov
# Michael L. Parks      mlparks@sandia.gov
# Stewart A. Silling    sasilli@sandia.gov
#
# ************************************************************************

import sys

def IsComment(line):
    if len(line) > 0 and (line[0] == '#' or line[0] == '\n'):
        return True
    return False

def CheckForAPrePro(line):
    if line[0] == "{":
        print("**** Parsing error, peridigm_to_yaml.py does not support aprepro commands")
        print("**** Remove aprepro commands manually and re-run peridigm_to_yaml.py")
        print("**** Failed to parse line {}".format(line))
        sys.exit(1)
    return

def ParseKeyValue(line):

    key = ""
    value = ""

    quote_positions = [i for i in range(len(line)) if line.startswith('"', i)]

    if len(quote_positions) == 0:
        pos = line.rfind(' ')
        if pos == -1:
            print("**** Parsing error, failed to locate white space separating key and value in line {}".format(line))
            sys.exit(1)
        key = line[:pos]
        value = line[pos:]
    elif len(quote_positions) == 2:
        if line[quote_positions[0]-1] != " ":
            print("**** Parsing error, expected white space to proceed quote in line {}".format(line))
            sys.exit(1)
        key = line[:quote_positions[0]-1]
        value = line[quote_positions[0]-1:]
    else:
        print("**** Parsing error, expected either zero or two quotes in line {}".format(line))
        sys.exit(1)

    return (key, value)

def CorrectBoolEntry(value):

    stripped_value = value.strip()
    if stripped_value == "\"True\"" or stripped_value == "\"true\"":
        value = " true\n"
    elif stripped_value == "\"False\"" or stripped_value == "\"false\"":
        value = " false\n"
    return value

def NumWhiteSpaces(line):
    index = 0
    while index < len(line) and line[index] == ' ':
        index += 1

    return index

if __name__ == "__main__":

    print("\n---- Peridigm to YAML Input Deck Converter\n")

    if len(sys.argv) < 2:
        print("Usage:  peridigm_to_yaml.py <input_deck.peridigm>\n")
        sys.exit(1)

    peridigm_file_name = sys.argv[1]
    if peridigm_file_name[-9:] != ".peridigm":
        print("**** Error:  Expected input file suffix to be \".peridigm\"\n")
        sys.exit(1)

    peridigm_file = open(peridigm_file_name)
    lines = peridigm_file.readlines()
    peridigm_file.close()

    yaml_file_name = peridigm_file_name[:-8] + "yaml"
    yaml_file = open(yaml_file_name, 'w')
    yaml_file.write("Peridigm:\n")

    for i in range(len(lines)):

        line = lines[i]
        next_line = "\n"
        if i+1 < len(lines):
            next_line = lines[i+1]

        CheckForAPrePro(line)

        if IsComment(line) == True:
            yaml_file.write(line)
        else:
            num_white_spaces = NumWhiteSpaces(line)
            num_white_spaces_next_line = NumWhiteSpaces(next_line)
            is_list_header = False
            if num_white_spaces_next_line > num_white_spaces:
                is_list_header = True
            if is_list_header == True:
                line = line[:-1] + ":" + "\n"
                line = "  " + line
                yaml_file.write(line)
            else:
                key, value = ParseKeyValue(line)
                value = CorrectBoolEntry(value)
                line = key + ":" + value
                line = "  " + line
                yaml_file.write(line)

    yaml_file.close()

    print("Original file: {}".format(peridigm_file_name))
    print("New file:      {}\n".format(yaml_file_name))
