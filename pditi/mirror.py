#! /usr/bin/env python

import os
import sys
import getopt
import string

class OutOfSourceSymlinks():
	def __init__(self, opts, args):
		
		# check for options
		src=None
		bin=None
		dirs=[]
		exclude_files=[".project", ".cproject","cmake.pditi.sh", "mirror.py"]
		exclude_dirs=['.svn', '.metadata','.plugins']
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
				exclude_files.append(a)
			elif o=="-h":
				usage()
				sys.exit(0)
			else:
				print "ERROR: invalid option(s)"
				usage()
				sys.exit(2)
		self._src = src
		self._bin = bin
		self._dirs = dirs
		self._exclude_files=exclude_files
		self._exclude_dirs =exclude_dirs

	def mirror(self):
		print 'Mirroring directory structure (with FILE links ONLY) from \'dirs\' in \'src\' to \'dirs\' in \'bin\''
		print 'src = ',self._src
		print 'bin = ',self._bin
		print 'dirs = ',self._dirs
		print 'exclude directories = ',self._exclude_dirs,'\n'
		print 'exclude files = ',self._exclude_files,'\n'

		# loop over input dirs and create structure
		for d in self._dirs:
			src_dir=os.path.join(self._src,d)
			bin_dir=os.path.join(self._bin,d)
			mirror(src_dir,bin_dir,self._exclude_dirs,self._exclude_files)

def mirror(src,bin,exclude_dirs,exclude_files):

	# check that top level 'bin' exists, if not then create
	doesNotExist=not os.path.exists(bin)
	if doesNotExist:
		os.mkdir(bin)

	# create links to files in top level 'bin'
	create_file_links(src,bin,'.',exclude_files)
	
	# walk directory structure and create folders/links 
	for root, dirs, files in os.walk(src):
		(drive,tail)=os.path.split(root)
		if tail in exclude_dirs:
			continue
#		print 'dir = ', root
		
		relpath=os.path.relpath(root,src)
		binpath=os.path.join(bin,relpath)
		
		# remove excluded dirs from directories
		for e in exclude_dirs:
			if e in dirs:
				dirs.remove(e)
				
		for d in dirs:
			new_dir=os.path.join(binpath,d)
			doesNotExist=not os.path.exists(new_dir)
			if doesNotExist:
	#			print '\tCreate nested dir = ',new_dir,'\n'
				os.mkdir(new_dir)

		# create links to all files in dirs
		create_file_links(root,binpath,dirs,exclude_files)


def create_file_links(src,bin,dirs,exclude):
	src_dirs=[]
	bin_dirs=[]
	for d in dirs:
		target_dir = os.path.join(src,d)
		link_dir   = os.path.join(bin,d)
		os.chdir(target_dir)
		src_files=os.listdir('.')
#		print 'Creating symbolic links for out-of-source build: dir = '+d+'\n'
		for f in src_files:
			if '~'==f[-1] or not os.path.isfile(f):
				continue
			elif f in exclude:
				continue
			target=os.path.join(target_dir,f)
			link  =os.path.join(link_dir,f)
#			print '\t\tCreate link to file = ',f,'\n'
#			print '\t\t',link,'-->',target,'\n'
			if os.path.lexists(link):
				print '\t\t\tremoving existing link: ',link
				os.remove(link)
			os.symlink(target,link)

		
def main():
	try:
		opts, args = getopt.getopt(sys.argv[1:], "e:s:b:d:h")
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)

	# Check options and get 'target' and 'link' directories
	oss = OutOfSourceSymlinks(opts,args)

	# mirror src dirs in bin; create file links in bin dirs to src dirs
	oss.mirror()
	


def usage():
	print """
	
	Name: cmake_symlinks.py
	
	Description: Mirrors directory structure of <target_directory> in <link_directory>.
	             Then it creates symbolic links of FILES (ONLY) from 'target' dir(s) to 'link' dir(s)
	
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


