#!/usr/bin/env python

import numpy as np
import sys


def clean_data(data):

    strip_comments = [ item.split('!')[0].split() for item in data]

    clean = [item for sublist in strip_comments for item in sublist]

    return clean



with open(sys.argv[-1],'r') as f:
    data = [ line.rstrip() for line in f if (line[0] != '!' and line != '\n')]

num_vars = int(data[25].split()[0])

flat_data = clean_data(data[26:])

var_names = np.array(flat_data[0:num_vars], dtype=np.str)

time_steps = np.array(flat_data[num_vars::num_vars+1], dtype=np.double)

data_vars = [ np.array(flat_data[num_vars+i::num_vars+1], dtype=np.double) for i in range(1,num_vars+1)]


