#!/usr/bin/env python

"""
image_to_txt.py:  Converts an image into a meshless grid in the Peridigm text file format."
"""

__author__ = "David Littlewood (djlittl@sandia.gov)"

# ************************************************************************
#
#
#                             Peridigm
#                 Copyright (2011) Sandia Corporation
#
# Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
# the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions?
# David J. Littlewood   djlittl@sandia.gov
# John A. Mitchell      jamitch@sandia.gov
# Michael L. Parks      mlparks@sandia.gov
# Stewart A. Silling    sasilli@sandia.gov
#
# ************************************************************************

import string
import sys
from PIL import Image

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print "\nUsage:  text_to_genesis.py <image_file> <pixel_edge_length> <extrusion_depth>\n"
        print "        The image file must be readable by the python PIL package."
        print "        pixel_edge_length is the physical length associated with a single pixel in the image."
        print "        extrusion_depth is the number of nodes to be created in the out-of-plane direction.\n"
        sys.exit(1)

    image_file_name = sys.argv[1]
    pixel_edge_length = float(sys.argv[2])
    depth = int(sys.argv[3])

    threshold_grayscale_value = 220

    img = Image.open('florida.tiff').convert('L')
    width, height = img.size

    txt_file_name = string.splitfields(image_file_name, '.')[0] + ".txt"

    print "\n--Image to Text--\n"
    print "  image file:", image_file_name
    print "  pixel edge length:", pixel_edge_length
    print "  image width:", width
    print "  image height:", height
    print "  extrusion depth:", depth

    data = list(img.getdata())
    data = [data[offset:offset+width] for offset in range(0, width*height, width)]

    # x y z block_number volume
    pts = []

    pixel_volume = pixel_edge_length*pixel_edge_length*pixel_edge_length

    for k in range(depth):
        for i in range(width):
            for j in range(height):
                pixel_value = data[height-j-1][i]
                if pixel_value < threshold_grayscale_value:
                    x = i*pixel_edge_length
                    y = j*pixel_edge_length
                    z = k*pixel_edge_length
                    block_id = 1
                    pts.append([x, y, z, block_id, pixel_volume])

    txt_file = open(txt_file_name, 'w')
    for pt in pts:
        txt_file.write(str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2]) + " " + str(pt[3]) + " " + str(pt[4]) + "\n")
    txt_file.close()

    print "\nText file written to", txt_file_name
    print "\nConvert to genesis with text_to_genesis.py\n"
