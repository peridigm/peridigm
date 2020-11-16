#! /usr/bin/env python

import math

if __name__ == "__main__":

    csv_file_name = "WaveInBar_MultiBlock.csv"
    csv_file = open(csv_file_name)
    lines = csv_file.readlines()
    csv_file.close()
    vals = lines[-2].split()

    time = float(vals[0].strip(","))
    Volume_Block_1 = float(vals[1].strip(","))
    Volume_Block_2 = float(vals[2].strip(","))
    Volume_Block_3 = float(vals[3].strip(","))
    Volume_Block_4 = float(vals[4].strip(","))
    Volume_Block_5 = float(vals[5].strip(","))
    Volume_Block_6 = float(vals[6].strip(","))
    Volume_Block_7 = float(vals[7].strip(","))
    Volume_Block_8 = float(vals[8].strip(","))
    Volume_Node_Set_11 = float(vals[9].strip(","))
    Volume_Node_Set_12 = float(vals[10].strip(","))
    Volume_Node_Set_13 = float(vals[11].strip(","))
    Volume_Node_Set_14 = float(vals[12].strip(","))
    Volume_Node_Set_15 = float(vals[13].strip(","))
    Volume_Node_Set_16 = float(vals[14].strip(","))
    Volume_Node_Set_17 = float(vals[15].strip(","))
    Volume_Node_Set_18 = float(vals[16].strip(","))
    Volume_Node_Set_20 = float(vals[17].strip(","))
    Volume_Node_Set_30 = float(vals[18].strip(","))
    Max_Volume_Node_Set_30 = float(vals[19].strip(","))
    Min_Volume_Node_Set_30 = float(vals[20].strip(","))

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

    print("\nChecking consistency of Block_Data and Node_Set_Data values...\n")

    truth_val = num_elem_block_1 * element_vol
    print("Block 1 results:   {} should equal {} and should equal {}".format(Volume_Block_1, Volume_Node_Set_11, truth_val))
    if abs(Volume_Block_1 - Volume_Node_Set_11) > tol:
        result = False
    if abs(Volume_Block_1 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_2 * element_vol
    print("Block 2 results:   {} should equal {} and should equal {}".format(Volume_Block_2, Volume_Node_Set_12, truth_val))
    if abs(Volume_Block_2 - Volume_Node_Set_12) > tol:
        result = False
    if abs(Volume_Block_2 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_3 * element_vol
    print("Block 3 results:   {} should equal {} and should equal {}".format(Volume_Block_3, Volume_Node_Set_13, truth_val))
    if abs(Volume_Block_3 - Volume_Node_Set_13) > tol:
        result = False
    if abs(Volume_Block_3 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_4 * element_vol
    print("Block 4 results:   {} should equal {} and should equal {}".format(Volume_Block_4, Volume_Node_Set_14, truth_val))
    if abs(Volume_Block_4 - Volume_Node_Set_14) > tol:
        result = False
    if abs(Volume_Block_4 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_5 * element_vol
    print("Block 5 results:   {} should equal {} and should equal {}".format(Volume_Block_5, Volume_Node_Set_15, truth_val))
    if abs(Volume_Block_5 - Volume_Node_Set_15) > tol:
        result = False
    if abs(Volume_Block_5 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_6 * element_vol
    print("Block 6 results:   {} should equal {} and should equal {}".format(Volume_Block_6, Volume_Node_Set_16, truth_val))
    if abs(Volume_Block_6 - Volume_Node_Set_16) > tol:
        result = False
    if abs(Volume_Block_6 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_7 * element_vol
    print("Block 7 results:   {} should equal {} and should equal {}".format(Volume_Block_7, Volume_Node_Set_17, truth_val))
    if abs(Volume_Block_7 - Volume_Node_Set_17) > tol:
        result = False
    if abs(Volume_Block_7 - truth_val) > tol:
        result = False
    truth_val = num_elem_block_8 * element_vol
    print("Block 8 results:   {} should equal {} and should equal {}".format(Volume_Block_8, Volume_Node_Set_18, truth_val))
    if abs(Volume_Block_8 - Volume_Node_Set_18) > tol:
        result = False
    if abs(Volume_Block_8 - truth_val) > tol:
        result = False

    # Node Set 20 should equal Block_1 + Block_2 + Block_3
    block_val = Volume_Block_1 + Volume_Block_2 + Volume_Block_3
    print("\nNode Set 20 results:  {} should equal {}".format(block_val, Volume_Node_Set_20))
    if abs(block_val - Volume_Node_Set_20) > tol:
        result = False

    # Node Set 30 is the entire model
    # The volume should be 4.2e-6
    block_val =  Volume_Block_1 + Volume_Block_2 + Volume_Block_3 + Volume_Block_4 + Volume_Block_5 + Volume_Block_6 + Volume_Block_7 + Volume_Block_8
    truth_val = 4.2e-6
    print("Node Set 30 results:  {} should equal {} and should equal {}".format(block_val, Volume_Node_Set_30, truth_val))
    if abs(block_val - Volume_Node_Set_30) > tol:
        result = False
    if abs(truth_val - Volume_Node_Set_30) > tol:
        result = False

    # The minimim volume over the entire model should be 8.0e-9
    truth_val = 8.0e-9
    print("\nMin element volume: {} should equal {}".format(Min_Volume_Node_Set_30, truth_val))
    if abs(truth_val - Min_Volume_Node_Set_30) > tol:
        result = False

    # The maxmimum volume over the entire model should be 8.0e-9
    truth_val = 8.0e-9
    print("Max element volume: {} should equal {}".format(Max_Volume_Node_Set_30, truth_val))
    if abs(truth_val - Max_Volume_Node_Set_30) > tol:
        result = False

    if result == True:
        print("\nTest Passed.\n")
    else:
        print("\nTest FAILED.\n")
