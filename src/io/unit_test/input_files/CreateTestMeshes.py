#! /usr/bin/env python

from math import sqrt

def crossProduct(a, b):
    
    return ( (a[1]*b[2]-a[2]*b[1]), (a[2]*b[0]-a[0]*b[2]), (a[0]*b[1]-a[1]*b[0]) )

def norm(a):

    return sqrt( a[0]*a[0] + a[1]*a[1] + a[2]*a[2] )

class poorMansDeviceIndependentRandomNumber:

    def __init__(self):
        # Random number seeds
        self.m_w = 1234
        self.m_z = 1004321
 
    def rand(self):
        # 32-bit result
        self.m_z = 36969 * (self.m_z & 65535) + (self.m_z >> 16)
        self.m_w = 18000 * (self.m_w & 65535) + (self.m_w >> 16)
        val = (self.m_z << 16) + self.m_w
        return val

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

def createDumbbellMesh():
    
    xLow = 0.0
    xHigh = 10.0
    xNum = 60
    yLow = 0.0
    yHigh = 10.0
    yNum = 80
    zLow = 0.0
    zHigh = 10.0
    zNum = 100
    volume = 10.0*((float(xNum)+1.0)/float(xNum))/float(xNum) * \
        10.0*((float(yNum)+1.0)/float(yNum))/float(yNum) * \
        10.0*((float(zNum)+1.0)/float(zNum))/float(zNum) 
    material = 1

    diagonalThreshold = 0.2
    d1 = (xLow, yLow, zLow)
    d2 = (xHigh, yHigh, zHigh)
    
    radius1 = 1.5
    center1 = (xHigh - radius1, yHigh - radius1, zHigh - radius1)

    radius2 = 0.75
    center2 = (xLow + radius2, yLow + radius2, zLow + radius2)

    mesh = []

    temp = norm( (d2[0]-d1[0], d2[1]-d1[1], d2[2]-d1[2]) )

    for i in range(xNum):
        for j in range(yNum):
            for k in range(zNum):
                x = xLow + i*(xHigh - xLow)/(xNum - 1)
                y = yLow + j*(yHigh - yLow)/(yNum - 1)
                z = zLow + k*(zHigh - zLow)/(zNum - 1)

                distanceFromDiagonal = norm( crossProduct( (x-d1[0], y-d1[1], z-d1[2]), (x-d2[0], y-d2[1], z-d2[2]) ) ) / temp
                distanceFromCenter1 = norm( (x-center1[0], y-center1[1], z-center1[2]) )
                distanceFromCenter2 = norm( (x-center2[0], y-center2[1], z-center2[2]) )

                if distanceFromDiagonal < diagonalThreshold or distanceFromCenter1 < radius1 or distanceFromCenter2 < radius2:
                    mesh.append( (x, y, z, material, volume) )

    return mesh

def createRandomMesh():
    
    low = 0.0
    high = 10.0
    mesh = []
    numPts = 8000
    volume = 0.5
    material = 1
    r = poorMansDeviceIndependentRandomNumber()

    for i in range(numPts):
        x = r.rand()
        y = r.rand()
        z = r.rand()
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

    mesh = createDumbbellMesh()
    fname = "dumbbell.txt"
    fout = open(fname, 'w')
    fout.write("# x y z block_id volume\n");
    for pt in mesh:
        fout.write(str(pt[0]) + " " + str(pt[1]) + " " + str(pt[2]) + " " + str(pt[3]) + " " + str(pt[4]) + "\n")
    fout.close()
    print "Wrote", fname, "\n"

    mesh = createRandomMesh()
    fname = "random.txt"
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
