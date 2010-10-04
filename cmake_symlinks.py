#! /usr/bin/env python

import os
import sys
import getopt

def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "e:s:b:d:h")
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)
	
	# check for options
	src=None
	bin=None
	dirs=[]
	exclude=[".project", ".cproject","cmake_symlinks.py"]
	# note default exclusions above
	for o, a in opts:
		if o=="-s" and src:
			print "symlinks ERROR"
			sys.exit("option [-s source directory] can only occur once on command line")
		elif o=="-s":
			src=a
		elif o=="-b" and bin:
			print "symlinks ERROR"
			sys.exit("option [-b binary directory] can only occur once on command line")
		elif o=="-b":
			bin=a
		elif o=="-d":
			dirs.append(a)
		elif o=="-e":
			exclude.append(a)
		elif o=="-h":
			usage()
			sys.exit(0)
		else:
			print "ERROR: invalid option(s)"
			usage()
			sys.exit(2)

	print 'symlinks creating symbolic links from \'dirs\' in \'src\' to \'dirs\' in \'bin\''
	print 'src = ',src
	print 'bin = ',bin
	print 'dirs = ',dirs
	print 'exclude files = ',exclude,'\n'
	src_dirs=[]
	bin_dirs=[]
	for d in dirs:
		target_dir = os.path.join(src,d)
		link_dir   = os.path.join(bin,d)
		os.chdir(target_dir)
		src_files=os.listdir('.')
		print 'Creating symbolic links for out-of-source build: dir = '+d+'\n'
		for f in src_files:
			if '~'==f[-1] or not os.path.isfile(f):
				continue
			elif f in exclude:
				continue
			target=os.path.join(target_dir,f)
			link  =os.path.join(link_dir,f)
#			print "target = ",target
#			print "link = ",link
			os.symlink(target,link)

def usage():
	print """
	
	Name: cmake_symlinks.py
	
	Description: creates symbolic links of FILES (ONLY) from 'target' dir to 'link' dir
	
	Usage: symlink -b target_directory -s link_directory -d directory_relative_to_link_directory [Options] 
	
	Options:
		-h            Print this message.
		-e filename   Exclude this file in any target directories.
		              Note default exclusion: ['.project', '.cproject', 'cmake_symlinks.py']
		              
	Example:
		Create links in PdITI binary directory for out-of-source build.
		Only create links in directories 'operator', 'intrepid' and 'operator/unitTest'.
		Exclude build.home.sh -- ie do not create a link to this file.
	
		$cmake_symlinks.py  -s /home/awesome/c++/eclipseProjects/pimp \\
		          -b /home/awesome/c++/eclipseProjects/pimp.build \\ 
		          -d operator -d intrepid -d operator/unitTests -e build.home.sh
	"""
	
if __name__=='__main__':
	main()


