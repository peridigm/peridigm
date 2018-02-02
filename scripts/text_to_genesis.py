#!/usr/bin/env python

"""
text_to_genesis.py:  Converts a meshfree discretization from the Peridigm text file format to the genesis/exodus file format."
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

import sys
import string
import os

# The following points to the location of exodus.py
path_to_exodus_py = 'trilinos_install_path/lib'

# The path above must also point to the netcdf libraries.  These libraries are not generally included in
# the Trilinos install.  One solution is to create symbolic links in this directory that point to the
# netcdf libraries (which you built as a TPL for Trilinos)

sys.path.append("/Users/djlittl/Software/seacas/GCC_4.9.4_THREAD_SAFE/lib")
import exodus

def read_line(file):
    """Scans the input file and ignores lines starting with a '#' or '\n'."""

    buff = file.readline()
    if len(buff) == 0: return None
    while buff[0] == '#' or buff[0] == '\n':
        buff = file.readline()
        if len(buff) == 0: return None
    return buff

if __name__ == "__main__":

    print "\n---- Text to Genesis\n"

    if len(sys.argv) < 2:
        print "Usage:  text_to_genesis.py <discretization_file.txt> <nodeset_1.txt> <nodeset_2.txt> ... <nodeset_n.txt>\n"
        print "The discretization file lists the nodes as (x, y, z, block_id, volume)"
        print "The node set files list the node numbers in each node set (1-based indexing)\n"
        sys.exit(1)

    textFileName = sys.argv[1]
    nodeSetFileNames = []
    for i in range(len(sys.argv)-2):
        nodeSetFileNames.append(sys.argv[i+2])

    print "Discretization file:"
    print " ", textFileName
    print
    print "Node set files:"
    for i in range(len(nodeSetFileNames)):
        print " ", nodeSetFileNames[i]
    print

    textFile = open(textFileName)

    # The text file contains (x, y, z, block_id, vol)
    X = []
    Y = []
    Z = []
    vol = []
    blocks = {}

    nodeId = 0
    buff = read_line(textFile)
    while buff != None:
        nodeId += 1
        vals = string.splitfields(buff)
        X.append(float(vals[0]))
        Y.append(float(vals[1]))
        Z.append(float(vals[2]))
        block_id = int(vals[3])
        vol.append(float(vals[4]))
        if block_id not in blocks.keys():
            blocks[block_id] = []
        blocks[block_id].append(nodeId)
        buff = read_line(textFile)
    textFile.close()

    print "Read", len(X), "nodes and", len(blocks), "block from", textFileName

    nodeSets = []
    for fileName in nodeSetFileNames:
        nodeSets.append([])
        node_set_id = len(nodeSets)
        nodeSetFile = open(fileName)
        buff = read_line(nodeSetFile)
        while buff != None:
            vals = string.splitfields(buff)
            for val in vals:
                node_id = int(val)
                nodeSets[node_set_id-1].append(node_id)
            buff = read_line(nodeSetFile)
        nodeSetFile.close()
        print "Read", len(nodeSets[node_set_id-1]), "node ids from", fileName

    # Write the Exodus II file
    exodusFileName = textFileName[:-4] + ".g"
    exodusTitle = "Text to Genesis translation of " + textFileName[:-4]

    # Remove old version of the exodus file, if it exists
    if os.path.exists(exodusFileName):
        os.remove(exodusFileName)

    # Open the output Exodus file
    exodusNumberOfDimensions = 3
    exodusNumberOfNodes = len(X)
    exodusNumberOfElements = len(X)
    exodusNumberOfBlocks = len(blocks)
    exodusNumberOfNodeSets = len(nodeSets)
    exodusNumberOfSideSets = 0
    exodusFile = exodus.exodus( exodusFileName,
                                'w',
                                'ctype',
                                exodusTitle,
                                exodusNumberOfDimensions,
                                exodusNumberOfNodes,
                                exodusNumberOfElements,
                                exodusNumberOfBlocks,
                                exodusNumberOfNodeSets,
                                exodusNumberOfSideSets )

    # Write the nodal coordinates
    coordNames = ["X", "Y", "Z"]
    exodusFile.put_coord_names(coordNames)
    exodusFile.put_coords(X, Y, Z)

    exodusBlockIds = []
    for key in blocks.keys():
        exodusBlockIds.append(key)

    # Write the element block information
    exodusElementTypes = ["SPHERE"]*exodusNumberOfBlocks
    exodusNumberOfElementsPerBlock = []
    for blockId in exodusBlockIds:
        exodusNumberOfElementsPerBlock.append(len(blocks[blockId]))
    exodusNumberOfNodesPerElement = [1]*exodusNumberOfBlocks
    exodusNumberOfAttributes = [2]*exodusNumberOfBlocks
    exodusDefineMaps = 0

    exodusFile.put_concat_elem_blk(exodusBlockIds,
                                   exodusElementTypes,
                                   exodusNumberOfElementsPerBlock,
                                   exodusNumberOfNodesPerElement,
                                   exodusNumberOfAttributes,
                                   exodusDefineMaps)

    # Create attribute arrays for the volume and radius
    # We follow a long-standing convention for legacy codes and record two attributes
    # 1) The so-called SPH radius, which is NOT the actual radius of the sphere and is NOT used by Peridigm
    # 2) The element volume
    for blockId in exodusBlockIds:

        numberOfPoints = len(blocks[blockId])
        exodusElementConnectivity = blocks[blockId]

        numberOfAttributesPerElement = 2
        exodusElementAttributes = [0.0]*numberOfAttributesPerElement*numberOfPoints
        index = 0
        for nodeId in blocks[blockId]:
            volume = vol[nodeId-1]
            sphRadius = volume**(1./3.)
            exodusElementAttributes[2*index] = sphRadius
            exodusElementAttributes[2*index+1] = volume
            index += 1

        # Write connectivity array
        exodusFile.put_elem_connectivity(blockId, exodusElementConnectivity)

        # Write element attributes
        exodusFile.put_elem_attr(blockId, exodusElementAttributes)

    # Write the node sets
    #   node_set_params (id, numSetNodes, numSetDistFacts)
    #   node_set_nodes
    #   node_set_dist_fact
    #
    for i in range(len(nodeSets)):
        distributionFactors = [1.0]*len(nodeSets[i])
        exodusFile.put_node_set_params(i+1, len(nodeSets[i]), len(distributionFactors))
        exodusFile.put_node_set(i+1, nodeSets[i])
        exodusFile.put_node_set_dist_fact(i+1, distributionFactors)

    # Close the output Exodus file
    exodusFile.close()

    print "\nData written to Exodus II file", exodusFileName, "\n"
