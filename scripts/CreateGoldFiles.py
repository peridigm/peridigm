#! /usr/bin/env python

import os
import shutil

if __name__ == "__main__":

    message =  "\nThis script creates a set of gold files from a set\n"
    message += "of standard output files.  The string '_gold' is inserted\n"
    message += "into the file names and all references to files within the\n"
    message += "*.pvd and *.pvtu files are altered to include '_gold'\n"
    print message

    nameChanges = {}
    pvdFiles = []
    pvtuFiles = []
    vtuFiles = []
    for file in os.listdir(os.getcwd()):
        if (file[-4:] == '.pvd' or file[-5:] == '.pvtu' or file[-4:] == '.vtu') and file.find('_gold') == -1:
            oldName = file
            if file[-4:] == '.pvd':
                newName = file[:-4] + '_gold' + file[-4:]
                pvdFiles.append(oldName)
            elif file[-5:] == '.pvtu':
                insertLocation = file.rfind('_t')
                newName = file[:insertLocation] + '_gold' + file[insertLocation:]
                pvtuFiles.append(oldName)
            elif file[-4:] == '.vtu':
                insertLocation = file.rfind('_t')
                newName = file[:insertLocation] + '_gold' + file[insertLocation:]
                vtuFiles.append(oldName)
            nameChanges[oldName] = newName

    for oldName in pvdFiles:
        oldFile = open(oldName)
        newFile = open(nameChanges[oldName], 'w')
        oldLine = oldFile.readline()
        while len(oldLine) != 0:
            insertLocation = oldLine.rfind('_t')
            newLine = oldLine[:insertLocation] + '_gold' + oldLine[insertLocation:]
            newFile.write(newLine)
            oldLine = oldFile.readline()
        newFile.close()

    for oldName in pvtuFiles:
        oldFile = open(oldName)
        newFile = open(nameChanges[oldName], 'w')
        oldLine = oldFile.readline()
        while len(oldLine) != 0:
            insertLocation = oldLine.rfind('_t')
            newLine = oldLine[:insertLocation] + '_gold' + oldLine[insertLocation:]
            newFile.write(newLine)
            oldLine = oldFile.readline()
        newFile.close()

    for oldName in vtuFiles:
        shutil.copyfile(oldName, nameChanges[oldName])

    print len(nameChanges), "gold files created.\n"
