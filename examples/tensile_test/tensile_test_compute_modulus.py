#! /usr/bin/env python

import string

if __name__ == "__main__":

    truth_modulus = 180.0 # GPa
    cross_sectional_area = 0.1999996  # cm^2

    infile = open("tensile_test.csv")
    lines = infile.readlines()
    infile.close()

    engineering_stress = 0.0
    engineering_strain = 0.0

    for line in lines:
        vals = string.splitfields(line)
        entries = []
        if len(vals) == 7:
            for val in vals:
                entries.append( float( string.strip(val, "',") ) )
            engineering_stress = (entries[6] - entries[5])*0.5/cross_sectional_area
            # convert from dyne/cm^2 to GPa
            engineering_stress = engineering_stress * 1.0e-10
            engineering_strain = (entries[3] - entries[4])/(entries[1] - entries[2])
            if engineering_strain != 0.0:
                modulus = engineering_stress/engineering_strain
                error = modulus - truth_modulus
                percent_error = 100.0*abs(truth_modulus - modulus)/truth_modulus

    print "\nSimulation results:"
    my_string =  "  Input modulus    = %.2f GPa\n" % truth_modulus
    my_string += "  Computed modulus = %.2f GPa\n" % modulus
    my_string += "  Error            =  %.2f GPa\n" % error
    print my_string
