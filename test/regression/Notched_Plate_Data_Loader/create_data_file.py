#!/usr/bin/env python

# This script requires exodus.py

import sys
import os
import math

# The following points to the location of exodus.py
path_to_exodus_py = 'trilinos_install_path/lib'

sys.path.append(path_to_exodus_py)

if sys.version_info >= (3, 0):
  import exodus3 as exodus
else:
  import exodus2 as exodus

def FluidConcentration(x, y, z, t):

    val = 0.0
    if t > 0.0:
        val = math.erfc((2.0-x)/(4.0*math.sqrt(t/0.001)))

    return val

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("\nUsage:  create_data_file.py <input.g> <output.e>\n")
        sys.exit(1)

    input_file_name = sys.argv[1]
    output_file_name = sys.argv[2]

    print("\n-- create_data_file.py --\n")
    print("Genesis input file: {}".format(input_file_name))
    print("Exodus output file: {}".format(output_file_name))

    if os.path.exists(output_file_name):
        os.remove(output_file_name)

    old_database = exodus.exodus(input_file_name, 'r')
    x, y, z = old_database.get_coords()
    new_database = old_database.copy(output_file_name)
    old_database.close()

    field_name = "Fluid_Concentration"
    num_nodes = len(x)
    field_values = [0.0]*num_nodes

    g_var_names = []
    n_var_names = [field_name]
    e_var_names = []
    exodus.add_variables(new_database, g_var_names, n_var_names, e_var_names)

    num_steps = 11
    final_time = 0.001
    time_step = final_time/(num_steps-1)

    for step in range(num_steps):
        time = step*time_step
        new_database.put_time(step+1, time)
        for i in range(num_nodes):
            field_values[i] = FluidConcentration(x[i], y[i], z[i], time)
        new_database.put_node_variable_values(field_name, step+1, field_values)

    new_database.close()
