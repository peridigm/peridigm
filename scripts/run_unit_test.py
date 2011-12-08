#! /usr/bin/env python

import sys
import os
import string
from subprocess import Popen

if __name__ == "__main__":

    executable = sys.argv[-1]
    base_name = string.splitfields(executable, '/')[-1]
    logfile = open(base_name + ".log", 'w')

    if not os.path.exists(executable):
        logfile.write("\nError:  " + executable + " not found!\n\n")
        sys.exit(1)

    command = []
    for i in range(len(sys.argv)-1):
        command.append(sys.argv[1+i])

    p = Popen(command, stdout=logfile, stderr=logfile)
    return_code = p.wait()

    sys.exit(return_code)
