#! /usr/bin/env python

import string
import math

if __name__ == "__main__":

    csv_file_name = "WaveInBar_MultiBlock.csv"
    csv_file = open(csv_file_name)
    lines = csv_file.readlines()
    csv_file.close()
    vals = string.splitfields(lines[-2])

    time = float(string.strip(vals[0],","))
    Volume_Block_1 = float(string.strip(vals[1],","))
    Volume_Block_2 = float(string.strip(vals[2],","))
    Volume_Block_3 = float(string.strip(vals[3],","))
    Volume_Block_4 = float(string.strip(vals[4],","))
    Volume_Block_5 = float(string.strip(vals[5],","))
    Volume_Block_6 = float(string.strip(vals[6],","))
    Volume_Block_7 = float(string.strip(vals[7],","))
    Volume_Block_8 = float(string.strip(vals[8],","))
    Volume_Node_Set_11 = float(string.strip(vals[9],","))
    Volume_Node_Set_12 = float(string.strip(vals[10],","))
    Volume_Node_Set_13 = float(string.strip(vals[11],","))
    Volume_Node_Set_14 = float(string.strip(vals[12],","))
    Volume_Node_Set_15 = float(string.strip(vals[13],","))
    Volume_Node_Set_16 = float(string.strip(vals[14],","))
    Volume_Node_Set_17 = float(string.strip(vals[15],","))
    Volume_Node_Set_18 = float(string.strip(vals[16],","))
    Volume_Node_Set_20 = float(string.strip(vals[17],","))
    Volume_Node_Set_30 = float(string.strip(vals[18],","))
    Max_Volume_Node_Set_30 = float(string.strip(vals[19],","))
    Min_Volume_Node_Set_30 = float(string.strip(vals[20],","))

    tol = 1.0e-12
    result = True

    element_vol = 8.0e-9
    num_elem_block_1 = 5
    num_elem_block_2 = 5
    num_elem_block_3 = 125
    num_elem_block_4 = 265
    num_elem_block_5 = 50
    num_elem_block_6 = 27
    num_elem_block_7 = 1
    num_elem_block_8 = 47

    # Each block is also a node set
    # Make sure the Block_Data and Node_Set_Data results match for each Block/Node Set pair

    print "\nChecking consistency of Block_Data and Node_Set_Data values...\n"

    truth_val = num_elem_block_1 * element_vol
    print "Block 1 results:   ", Volume_Block_1, "should equal", Volume_Node_Set_11, "and should equal", truth_val
    if abs(Volume_Block_1 - Volume_Node_Set_11) > tol:
        result = False
    if abs(Volume_Block_1 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_2 * element_vol
    print "Block 2 results:   ", Volume_Block_2, "should equal", Volume_Node_Set_12, "and should equal", truth_val
    if abs(Volume_Block_2 - Volume_Node_Set_12) > tol:
        result = False
    if abs(Volume_Block_2 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_3 * element_vol
    print "Block 3 results:   ", Volume_Block_3, "should equal", Volume_Node_Set_13, "and should equal", truth_val
    if abs(Volume_Block_3 - Volume_Node_Set_13) > tol:
        result = False
    if abs(Volume_Block_3 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_4 * element_vol
    print "Block 4 results:   ", Volume_Block_4, "should equal", Volume_Node_Set_14, "and should equal", truth_val
    if abs(Volume_Block_4 - Volume_Node_Set_14) > tol:
        result = False
    if abs(Volume_Block_4 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_5 * element_vol
    print "Block 5 results:   ", Volume_Block_5, "should equal", Volume_Node_Set_15, "and should equal", truth_val
    if abs(Volume_Block_5 - Volume_Node_Set_15) > tol:
        result = False
    if abs(Volume_Block_5 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_6 * element_vol
    print "Block 6 results:   ", Volume_Block_6, "should equal", Volume_Node_Set_16, "and should equal", truth_val
    if abs(Volume_Block_6 - Volume_Node_Set_16) > tol:
        result = False
    if abs(Volume_Block_6 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_7 * element_vol
    print "Block 7 results:   ", Volume_Block_7, "should equal", Volume_Node_Set_17, "and should equal", truth_val
    if abs(Volume_Block_7 - Volume_Node_Set_17) > tol:
        result = False
    if abs(Volume_Block_7 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_8 * element_vol
    print "Block 8 results:   ", Volume_Block_8, "should equal", Volume_Node_Set_18, "and should equal", truth_val
    if abs(Volume_Block_8 - Volume_Node_Set_18) > tol:
        result = False
    if abs(Volume_Block_8 - truth_val) > tol:
        result = False

    # Node Set 20 should equal Block_1 + Block_2 + Block_3
    block_val = Volume_Block_1 + Volume_Block_2 + Volume_Block_3
    print "\nNode Set 20 results:  ", block_val, "should equal", Volume_Node_Set_20
    if abs(block_val - Volume_Node_Set_20) > tol:
        result = False

    # Node Set 30 is the entire model
    # The volume should be 4.2e-6
    block_val =  Volume_Block_1 + Volume_Block_2 + Volume_Block_3 + Volume_Block_4 + Volume_Block_5 + Volume_Block_6 + Volume_Block_7 + Volume_Block_8
    truth_val = 4.2e-6
    print "Node Set 30 results:  ", block_val, "should equal", Volume_Node_Set_30, "and should equal", truth_val
    if abs(block_val - Volume_Node_Set_30) > tol:
        result = False
    if abs(truth_val - Volume_Node_Set_30) > tol:
        result = False

    # The minimim volume over the entire model should be 8.0e-9
    truth_val = 8.0e-9
    print "\nMin element volume:", Min_Volume_Node_Set_30, "should equal", truth_val
    if abs(truth_val - Min_Volume_Node_Set_30) > tol:
        result = False

    # The maxmimum volume over the entire model should be 8.0e-9
    truth_val = 8.0e-9
    print "Max element volume:", Max_Volume_Node_Set_30, "should equal", truth_val
    if abs(truth_val - Max_Volume_Node_Set_30) > tol:
        result = False

    if result == True:
        print "\nTest Passed.\n"
    else:
        print "\nTest FAILED.\n"
