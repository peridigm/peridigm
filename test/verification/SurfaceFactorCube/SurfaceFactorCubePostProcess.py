#! /usr/bin/env python

import string

def read_line(file):
    """Scans the input file and ignores lines starting with a '#' or '\n'."""
    
    buff = file.readline()
    if len(buff) == 0: return None
    while buff[0] == '#' or buff[0] == '\n':
        buff = file.readline()
        if len(buff) == 0: return None
    return buff

if __name__=='__main__':

    SCF_Unrotated_Block_Center = 0.0
    SCF_Unrotated_Block_Face = 0.0
    SCF_Unrotated_Block_Edge = 0.0
    SCF_Unrotated_Block_Corner = 0.0
    SCF_Rotated_Block_Center = 0.0
    SCF_Rotated_Block_Face = 0.0
    SCF_Rotated_Block_Edge = 0.0
    SCF_Rotated_Block_Corner = 0.0

    infileName = "SurfaceFactorCube.csv"
    dataFile = open(infileName)
    buff = read_line(dataFile)
    while buff != None:
        vals = string.splitfields(buff)
        for i in range(len(vals)):
            vals[i] = float(string.rstrip(vals[i], ','))
        if len(vals) > 1:
            SCF_Unrotated_Block_Center = float(vals[1])
            SCF_Unrotated_Block_Face   = float(vals[2])
            SCF_Unrotated_Block_Edge   = float(vals[3])
            SCF_Unrotated_Block_Corner = float(vals[4])
            SCF_Rotated_Block_Center   = float(vals[5])
            SCF_Rotated_Block_Face     = float(vals[6])
            SCF_Rotated_Block_Edge     = float(vals[7])
            SCF_Rotated_Block_Corner   = float(vals[8])
        buff = read_line(dataFile)
    dataFile.close()

    print
    print "Surface Correction Factors"
    print
    print "Rotation    Location    Analytic_Value   Computed_Value"
    print "-------------------------------------------------------"
    print "unrotated   center      1.0              %.6f" % SCF_Unrotated_Block_Center
    print "unrotated   face        TBD              %.6f" % SCF_Unrotated_Block_Face
    print "unrotated   edge        TBD              %.6f" % SCF_Unrotated_Block_Edge
    print "unrotated   corner      TBD              %.6f" % SCF_Unrotated_Block_Corner
    print "rotated     center      1.0              %.6f" % SCF_Rotated_Block_Center
    print "rotated     face        TBD              %.6f" % SCF_Rotated_Block_Face
    print "rotated     edge        TBD              %.6f" % SCF_Rotated_Block_Edge
    print "rotated     corner      TBD              %.6f" % SCF_Rotated_Block_Corner
    print
