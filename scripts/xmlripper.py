#!/usr/bin/env python

''' strip information from an xml file and import in to a txt file name, value, type seperated by comma's
'''
import xml.etree.ElementTree as et, sys

def main(flname):
	tree = et.parse(flname + '.xml')
	root = tree.getroot()
	fl = open(flname+'.txt','w+')
	for child in root.iter():
		if(child.get('name')):
			if(child.get('value')):
				fl.write('  ' + child.attrib['name'] + ' "' + child.attrib['value'] + '"\n')
			else:
				fl.write(child.attrib['name'] + '\n')
	fl.close()
def __init__():
	flname = sys.argv[-1]
	print(flname)
	main(flname)
if __name__ == '__main__':
	__init__()
