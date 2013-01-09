#! /usr/bin/env python

def createEquallySpacedCubeMesh1000():
    
    xLow = 0.5
    xHigh = 9.5
    xNum = 10
    yLow = 0.5
    yHigh = 9.5
    yNum = 10
    zLow = 0.5
    zHigh = 9.5
    zNum = 10
    volume = 1.0
    material = 1

    mesh = []

    for i in range(xNum):
        for j in range(yNum):
            for k in range(zNum):
                x = xLow + i*(xHigh - xLow)/(xNum - 1)
                y = yLow + j*(yHigh - yLow)/(yNum - 1)
                z = zLow + k*(zHigh - zLow)/(zNum - 1)
                mesh.append( (x, y, z, material, volume) )

    return mesh

if __name__ == "__main__":

    message =  "\n--CreateTestMesh.py --\n\n"
    message += "This script generates a series of text files containing meshes\n"
    message += "for testing neighbor search algorithms.  The format of the mesh\n"
    message += "files is valid input for Peridigm, although for the search tests\n"
    message += "only the coordinates are used.\n"
    print message

    mesh = createEquallySpacedCubeMesh1000()
    fout = open("cube_1000.txt", 'w')
    fout.write("# x y z block_id volume\n");
    for pt in mesh:
        fout.write(str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2]) + " " + str(pt[3]) + " " + str(pt[4]) + "\n")
    fout.close()
    print "Wrote cube_1000.txt\n"

    print "Complete.\n"
