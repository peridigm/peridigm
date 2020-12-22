#! /usr/bin/env python

if __name__ == "__main__":

    # User settings
    bar_depth = 0.001
    aspect_ratio = 10
    num_points_through_depth = 10

    bar_length = aspect_ratio * bar_depth
    bar_volume = bar_length * bar_depth * bar_depth
    element_size = bar_depth/num_points_through_depth
    num_points_along_length = num_points_through_depth * aspect_ratio
    total_num_points = num_points_along_length * num_points_through_depth * num_points_through_depth
    element_volume = bar_volume / total_num_points
    horizon = 3 * element_size

    print("\nGenerating bar mesh:")
    print("  bar length {}".format(bar_length))
    print("  bar depth {}".format(bar_depth))
    print("  total number of cells {}".format(total_num_points))
    print("  number of cells through depth {}".format(num_points_through_depth))
    print("  number of cells along length {}".format(num_points_along_length))
    print("  element size {}".format(element_size))
    print("  horizon {}".format(horizon))

    half_elem_size = element_size / 2.0

    data = []
    node_set_1 = []
    node_set_2 = []
    id = 1
    for k in range(num_points_along_length):
        for j in range(num_points_through_depth):
            for i in range(num_points_through_depth):

                # Record the data for this element
                x = k * element_size + half_elem_size
                y = -bar_depth / 2.0 + i * element_size + half_elem_size
                z = -bar_depth / 2.0 + j * element_size + half_elem_size
                block_id = 1
                data.append([x, y, z, block_id, element_volume])

                if x < bar_length - 2.0 * horizon:
                    node_set_1.append(id)
                if x > bar_length - 2.0 * horizon:
                    node_set_2.append(id)
                id += 1

    with open("wave_in_bar.txt", 'w') as mesh_file:
        for datum in data:
            mesh_file.write(str(datum[0]) + " " + str(datum[1]) + " " + str(datum[2]) + " " + str(datum[3]) + " " + str(datum[4]) + "\n")

    with open("node_set_1.txt", 'w') as node_set_1_file:
        for id in node_set_1:
            node_set_1_file.write(str(id) + "\n")

    with open("node_set_2.txt", 'w') as node_set_2_file:
        for id in node_set_2:
            node_set_2_file.write(str(id) + "\n")

    print("\nDiscretization written to wave_in_bar.txt, node_set_1.txt, and node_set_2.txt\n")

