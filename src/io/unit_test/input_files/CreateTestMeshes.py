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

def createEquallySpacedCubeMesh8000():
    
    xLow = 0.25
    xHigh = 9.75
    xNum = 20
    yLow = 0.25
    yHigh = 9.75
    yNum = 20
    zLow = 0.25
    zHigh = 9.75
    zNum = 20
    volume = 0.5
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

def createEquallySpacedCubeMesh27000():
    
    xLow = 1.0/6.0
    xHigh = 10.0 - 1.0/6.0
    xNum = 30
    yLow = 1.0/6.0
    yHigh = 10.0 - 1.0/6.0
    yNum = 30
    zLow = 1.0/6.0
    zHigh = 10.0 - 1.0/6.0
    zNum = 30
    volume = 1.0/27.0
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

def writeTestMeshes():

    mesh = createEquallySpacedCubeMesh1000()
    fname = "cube_1000.txt"
    fout = open(fname, 'w')
    fout.write("# x y z block_id volume\n");
    for pt in mesh:
        fout.write(str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2]) + " " + str(pt[3]) + " " + str(pt[4]) + "\n")
    fout.close()
    print "Wrote", fname, "\n"

    mesh = createEquallySpacedCubeMesh8000()
    fname = "cube_8000.txt"
    fout = open(fname, 'w')
    fout.write("# x y z block_id volume\n");
    for pt in mesh:
        fout.write(str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2]) + " " + str(pt[3]) + " " + str(pt[4]) + "\n")
    fout.close()
    print "Wrote", fname, "\n"

    mesh = createEquallySpacedCubeMesh27000()
    fname = "cube_27000.txt"
    fout = open(fname, 'w')
    fout.write("# x y z block_id volume\n");
    for pt in mesh:
        fout.write(str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2]) + " " + str(pt[3]) + " " + str(pt[4]) + "\n")
    fout.close()
    print "Wrote", fname, "\n"

    return

if __name__ == "__main__":

    message =  "\n--CreateTestMesh.py --\n\n"
    message += "This script generates a series of text files containing meshes\n"
    message += "for testing neighbor search algorithms.  The format of the mesh\n"
    message += "files is valid input for Peridigm, although for the search tests\n"
    message += "only the coordinates are used.\n"
    print message

    writeTestMeshes()

    print "Complete.\n"
