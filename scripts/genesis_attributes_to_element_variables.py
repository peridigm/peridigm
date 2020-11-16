#!/usr/bin/env python

# This script requires exodus.py

import sys
import os
import math

# The following points to the location of exodus.py
path_to_exodus_py = 'trilinos_install_path/lib'

# The path above must also point to the netcdf libraries.  These libraries are not generally included in
# the Trilinos install.  One solution is to create symbolic links in this directory that point to the
# netcdf libraries (which you built as a TPL for Trilinos)

sys.path.append(path_to_exodus_py)

if sys.version_info >= (3, 0):
  import exodus3 as exodus
else:
  import exodus2 as exodus

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("\nUsage:  genesis_attributes_to_node_variables.py <input.g> <output.g>\n")
        sys.exit(1)

    input_file_name = sys.argv[1]
    output_file_name = sys.argv[2]

    print("\n-- genesis_attributes_to_node_variables.py --\n")
    print("Genesis input file:{}".format(input_file_name))
    print("Genesis output file:{}".format(output_file_name))

    if os.path.exists(output_file_name):
        os.remove(output_file_name)

    old_database = exodus.exodus(input_file_name, 'r')
    block_ids = old_database.get_elem_blk_ids()

    elem_var_names = {}
    elem_var_values = {}

    for block_id in block_ids:

        attr_names = old_database.get_element_attribute_names(block_id)
        for i in range(len(attr_names)):
            if attr_names[i] == '':
                attr_names[i] = "attribute_" + str(i+1)

        # For the case of sphere meshes, there are two unnamed attributes which are the sph_radius and the volume.
        # In this case, add radius as a variable.
        add_radius = False
        if len(attr_names) == 2:
            if attr_names[0] == 'attribute_1' and attr_names[1] == 'attribute_2':
                add_radius = True

        num_attr = len(attr_names)
        attr = old_database.get_elem_attr(block_id)

        elem_var_names[block_id] = attr_names
        elem_var_values[block_id] = []
        for i in range(num_attr):
            elem_var_values[block_id].append([])

        if add_radius == True:
            elem_var_names[block_id].append('radius')
            elem_var_values[block_id].append([])

        attr_index = 0
        for i_elem in range( old_database.num_elems_in_blk(block_id) ):
            for i_attr in range(num_attr):
                elem_var_values[block_id][i_attr].append(attr[attr_index])
                attr_index += 1
            if add_radius == True:
                volume = elem_var_values[block_id][1][-1]
                radius = math.pow((3.0*volume)/(4.0*math.pi), 1.0/3.0)
                elem_var_values[block_id][-1].append(radius)

    new_database = old_database.copy(output_file_name)
    old_database.close()

    elem_names_and_blocks = {}
    for block_id in elem_var_names.keys():
        for name in elem_var_names[block_id]:
            if name not in elem_names_and_blocks:
                elem_names_and_blocks[name] = []
            elem_names_and_blocks[name].append(block_id)

    g_var_names = []
    n_var_names = []
    e_var_names = []
    for var_name in elem_names_and_blocks.keys():
        e_var_names.append((var_name, elem_names_and_blocks[var_name]))

    new_database.put_time(1, 0.0)
    exodus.add_variables(new_database, g_var_names, n_var_names, e_var_names)

    for block_id in block_ids:
        for i_var in range(len(elem_var_names[block_id])):
            new_database.put_element_variable_values(block_id, elem_var_names[block_id][i_var], 1, elem_var_values[block_id][i_var])

    new_database.close()

