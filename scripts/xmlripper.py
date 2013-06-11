#!/usr/bin/env python

''' parse Teuchos compatible xml based ParameterList to txt based ParameterList: standard xml files accepted also No Error Handling
'''
import xml.etree.ElementTree as et, sys


def sTn(val):
	if val.isdigit():
		return float(val)
	else:
		return val

def main(flname):
	tree = et.parse(flname)
	root = tree.getroot()
	fl = open(flname.replace('.xml', '.txt'),'w+')
	dataVec = []
	i = 1
	for param in root.findall("./"):
#		print param.attrib
		if param.get('value'):
			if param.attrib['type'] == 'string' or param.attrib['type'] == 'bool':
				fl.write(param.attrib['name'] + ' "' + param.attrib['value'] + '"\n')
			elif param.attrib['type'] == 'int' or param.attrib['type'] == 'double':
				fl.write(param.attrib['name'] + ' ' + param.attrib['value'] + '\n')
		if not param.get('value'):
			fl.write(param.attrib['name'] + '\n')
			j = 1
			for subparam in root.findall("./ParameterList[%s]/"%i):
#				print subparam.attrib
				if subparam.get('value'):
					if subparam.attrib['type'] == 'string' or subparam.attrib['type'] == 'bool':
						fl.write('  ' + subparam.attrib['name'] + ' "' + subparam.attrib['value'] + '"\n')
					elif subparam.attrib['type'] == 'int' or subparam.attrib['type'] == 'double':
						fl.write('  ' + subparam.attrib['name'] + ' ' + subparam.attrib['value'] + '\n')
				if not subparam.get('value'):
					fl.write('  ' + subparam.attrib['name'] + '\n')
					for sub_subparam in root.findall("./ParameterList[%s]/ParameterList[%s]/" %(i, j)):
#						print sub_subparam.attrib
						if sub_subparam.get('value'):
							if sub_subparam.attrib['type'] == 'string' or sub_subparam.attrib['type'] == 'bool':
								fl.write('  ' + '  ' + sub_subparam.attrib['name'] + ' "' + sub_subparam.attrib['value'] + '"\n')
							elif sub_subparam.attrib['type'] == 'int' or sub_subparam.attrib['type'] == 'double':
								fl.write('  ' + '  ' + sub_subparam.attrib['name'] + ' ' + sub_subparam.attrib['value'] + '\n')
					j = j+1
			i = i+1
def __init__():
	flname = sys.argv[-1]
	print(flname)
	main(flname)
if __name__ == '__main__':
	__init__()
