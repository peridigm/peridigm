#! /usr/bin/env python

import math
import random

if __name__ == "__main__":

    # Cylinder dimensions
    cylinder_height = 0.100
    cylinder_outer_radius = 0.025
    cylinder_inner_radius = 0.020

    # Number of cells through the wall thickness
    num_points_radius = 5

    element_size = (cylinder_outer_radius - cylinder_inner_radius)/num_points_radius
    num_points_circumference = int((2.0 * math.pi * (cylinder_outer_radius + cylinder_inner_radius)/2.0)/element_size + 0.5)
    num_points_height = int(cylinder_height/element_size + 0.5)
    arc_length = 2.0*math.pi/num_points_circumference

    print("\nGenerating cylinder mesh:")
    print("  number of points through wall thickness: {}".format(num_points_radius))
    print("  number of points around circumference: {}".format(num_points_circumference))
    print("  number of points along cylinder length: {}".format(num_points_height))
    print("  total number of points: {}".format(num_points_radius*num_points_circumference*num_points_height))

    data = []
    total_volume = 0.0
    for k in range(num_points_height):
        for j in range(num_points_circumference):
            for i in range(num_points_radius):

                # Geometric calculations make use of the fact that the element is a segment of a ring
                # The radial distance of the centroid of a circular sector from the origin is given by:
                #    (2 * r * sin(alpha)) / (3 * alpha)
                # where alpha is half the angle of the sector

                # Compute the centroid
                elem_inner_radius = cylinder_inner_radius + i*(cylinder_outer_radius - cylinder_inner_radius)/num_points_radius
                elem_outer_radius = cylinder_inner_radius + (i+1)*(cylinder_outer_radius - cylinder_inner_radius)/num_points_radius
                alpha = arc_length/2.0
                large_circular_sector_centroid_distance = 2.0*elem_outer_radius*math.sin(alpha)/(3.0*alpha)
                small_circular_sector_centroid_distance = 2.0*elem_inner_radius*math.sin(alpha)/(3.0*alpha)
                large_circular_sector_area = alpha*elem_outer_radius*elem_outer_radius
                small_circular_sector_area = alpha*elem_inner_radius*elem_inner_radius
                elem_area = large_circular_sector_area - small_circular_sector_area
                elem_centroid_distance = (large_circular_sector_centroid_distance*large_circular_sector_area - small_circular_sector_centroid_distance*small_circular_sector_area)/elem_area
                elem_centroid_angle = (j + 0.5)*arc_length
                elem_height = (k + 0.5)*cylinder_height/num_points_height
                elem_vol = elem_area*cylinder_height/num_points_height
                total_volume += elem_vol

                # Record the data for this element
                x = elem_centroid_distance * math.cos(elem_centroid_angle)
                y = elem_centroid_distance * math.sin(elem_centroid_angle)
                z = elem_height
                block_id = 1
                data.append([x, y, z, block_id, elem_vol])

    # Perturb the node locations randomly to break symmetry in the model
    magnitude = 0.001 * element_size
    random_number_generator_seed = 42
    random.seed(random_number_generator_seed)
    for i in range(len(data)):
        data[i][0] += (2.0*random.random() - 1.0) * magnitude
        data[i][1] += (2.0*random.random() - 1.0) * magnitude
        data[i][2] += (2.0*random.random() - 1.0) * magnitude

    mesh_file = open("fragmenting_cylinder.txt", 'w')
    for datum in data:
        mesh_file.write(str(datum[0]) + " " + str(datum[1]) + " " + str(datum[2]) + " " + str(datum[3]) + " " + str(datum[4]) + "\n")
    mesh_file.close()

    # Nodeset needed for initial conditions. Nodeset indexed from 1.
    nodeset_file = open("fragmenting_cylinder_nodeset.txt", 'w')
    for i in range(len(data)):
        nodeset_file.write(str(i+1) + "\n")
    nodeset_file.close()

#    correct_total_volume = (math.pi*cylinder_outer_radius*cylinder_outer_radius - math.pi*cylinder_inner_radius*cylinder_inner_radius)*cylinder_height
#    print "\nSanity check on volume:", total_volume, "=?", correct_total_volume

    print("\nDiscretization written to fragmenting_cylinder.txt\n")
    print("Node set written to fragmenting_cylinder_nodeset.txt\n")
