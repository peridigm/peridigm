#! /usr/bin/env python

import string
import math

if __name__ == "__main__":

    csv_file_name = "NeighborhoodVolume.csv"
    csv_file = open(csv_file_name)
    lines = csv_file.readlines()
    csv_file.close()
    vals = string.splitfields(lines[-2])
    horizon_A = float(string.strip(vals[-6],","))
    volume_A = float(string.strip(vals[-5],","))
    truth_volume_A = 4.0*math.pi*horizon_A*horizon_A*horizon_A/3.0
    horizon_B = float(string.strip(vals[-4],","))
    volume_B = float(string.strip(vals[-3],","))
    truth_volume_B = 4.0*math.pi*horizon_B*horizon_B*horizon_B/3.0
    horizon_C = float(string.strip(vals[-2],","))
    volume_C = float(string.strip(vals[-1],","))
    truth_volume_C = 4.0*math.pi*horizon_C*horizon_C*horizon_C/3.0
    print
    print "      Volume Point A =", volume_A
    print "Truth Volume Point A =", truth_volume_A
    print
    print "      Volume Point B =", volume_B
    print "Truth Volume Point B =", truth_volume_B
    print
    print "      Volume Point C =", volume_C
    print "Truth Volume Point C =", truth_volume_C
    print
