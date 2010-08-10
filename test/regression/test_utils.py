#! /usr/bin/env python

def read_line(file):
    """Scans the input file and ignores lines starting with a '#' or '\n'."""
    
    buff = file.readline()
    if len(buff) == 0: return None
    while buff[0] == '#' or buff[0] == '\n':
        buff = file.readline()
        if len(buff) == 0: return None
    return buff

def fequal_tol(a, b, tol):
    if abs(a-b) < tol:
        return True
    else:
        return False
