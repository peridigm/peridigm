#! /usr/bin/env python

import math

if __name__ == "__main__":

    csv_file_name = "NeighborhoodVolume.csv"
    csv_file = open(csv_file_name)
    lines = csv_file.readlines()
    csv_file.close()
    vals = lines[-2].split()
    horizon_A = float(vals[-6].strip(","))
    volume_A = float(vals[-5].strip(","))
    truth_volume_A = 4.0*math.pi*horizon_A*horizon_A*horizon_A/3.0
    horizon_B = float(vals[-4].strip(","))
    volume_B = float(vals[-3].strip(","))
    truth_volume_B = 4.0*math.pi*horizon_B*horizon_B*horizon_B/3.0
    horizon_C = float(vals[-2].strip(","))
    volume_C = float(vals[-1].strip(","))
    truth_volume_C = 4.0*math.pi*horizon_C*horizon_C*horizon_C/3.0
    print("")
    print("      Volume Point A = {}".format(volume_A))
    print("Truth Volume Point A = {}".format(truth_volume_A))
    print("")
    print("      Volume Point B = {}".format(volume_B))
    print("Truth Volume Point B = {}".format(truth_volume_B))
    print("")
    print("      Volume Point C = {}".format(volume_C))
    print("Truth Volume Point C = {}".format(truth_volume_C))
    print("")
