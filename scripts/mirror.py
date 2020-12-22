#! /usr/bin/env python

import os
import sys
import getopt

class OutOfSourceSymlinks():
  def __init__(self, opts, args):
    """
      Collect all command line arguments and set the related options
    """
    # check for options
    src=None
    bin=None
    dirs=[]
    exclude_files=[".project", ".cproject","cmake_symlinks.py", "mirror.py"]
    exclude_dirs=['.svn', '.metadata','.plugins']
    # note default exclusions above
    for o, a in opts:
      if o=="-s" and src:
        print("symlinks ERROR")
        sys.exit("option [-s source directory] can only occur once on command line")
      elif o=="-s":
        src=a
      elif o=="-b" and bin:
        print("symlinks ERROR")
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
        print("ERROR: invalid option(s)")
        usage()
        sys.exit(2)
    self._src = src
    self._bin = bin
    self._dirs = dirs
    self._exclude_files=exclude_files
    self._exclude_dirs =exclude_dirs

  def mirror(self):
    print('Mirroring directory structure (with FILE links ONLY) from \'dirs\' in \'src\' to \'dirs\' in \'bin\'')
    print('src = {}'.format(self._src))
    print('bin = {}'.format(self._bin))
    print('dirs = {}'.format(self._dirs))
    print('exclude directories = {}'.format(self._exclude_dirs))
    print('exclude files = {}'.format(self._exclude_files))

    # loop over input dirs and create structure
    num_removed_links = 0
    for d in self._dirs:
      src_dir=os.path.join(self._src,d)
      bin_dir=os.path.join(self._bin,d)
      num_removed_links += mirror(src_dir,bin_dir,self._exclude_dirs,self._exclude_files)

    print('total number of symbolic links removed/replaced = {}\n'.format(num_removed_links))

def mirror(src,bin,exclude_dirs,exclude_files):

  # check that top level 'bin' exists, if not then create
  doesNotExist=not os.path.exists(bin)
  if doesNotExist:
    os.mkdir(bin)

  # create links to files in top level 'bin'
  num_removed_links = create_file_links(src,bin,'.',exclude_files)
  
  # walk directory structure and create folders/links 
  for root, dirs, files in os.walk(src):
    (drive,tail)=os.path.split(root)
    if tail in exclude_dirs:
      continue
#    print('dir = {}'.format(root))
    
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
  #      print('\tCreate nested dir = {}\n'.format(new_dir))
        os.mkdir(new_dir)

    # create links to all files in dirs
    num_removed_links += create_file_links(root,binpath,dirs,exclude_files)

  return num_removed_links

def create_file_links(src,bin,dirs,exclude):
  num_removed_links = 0
  for d in dirs:
    target_dir = os.path.join(src,d)
    link_dir   = os.path.join(bin,d)
    os.chdir(target_dir)
    src_files=os.listdir('.')
#    print('Creating symbolic links for out-of-source build: dir = {}\n'.format(d))
    for f in src_files:
      if '~'==f[-1] or not os.path.isfile(f):
        continue
      elif f in exclude:
        continue
      target=os.path.join(target_dir,f)
      link  =os.path.join(link_dir,f)
#      print('\t\tCreate link to file = {}\n'.format(f))
#      print('\t\t {} --> {}\n'.format(link,target)))
      if os.path.lexists(link):
        os.remove(link)
        num_removed_links += 1
      os.symlink(target,link)

  return num_removed_links

def main():
  """
    Main routine
  """
  try:
    opts, args = getopt.getopt(sys.argv[1:], "e:s:b:d:h")
  except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(2)

  # Check options and get 'target' and 'link' directories
  oss = OutOfSourceSymlinks(opts,args)

  # mirror src dirs in bin; create file links in bin dirs to src dirs
  oss.mirror()


def usage():
  print("""
  
  Name: cmake_symlinks.py
  
  Description: Mirrors directory structure of <target_directory> in <link_directory>.
               Then it creates symbolic links of FILES (ONLY) from 'target' dir(s) to 'link' dir(s)
  
  Usage: symlink -b target_directory -s link_directory -d directory_relative_to_link_directory [Options] 
  
  Options:
    -h            Print this message.
    -e filename   Exclude this file in any target directories.
                  Note default exclusion: ['.project', '.cproject', 'cmake_symlinks.py']

  Example:
    Create links in binary directory for out-of-source build.
    Only create links in directories 'operator', 'intrepid' and 'operator/unitTest'.
    Exclude build.home.sh -- ie do not create a link to this file.
  
    $cmake_symlinks.py  -s /home/awesome/c++/eclipseProjects/binary \\
              -b /home/awesome/c++/eclipseProjects/pimp.build \\ 
              -d operator -d intrepid -d operator/unitTests -e build.home.sh
  """)
  
if __name__=='__main__':
  main()


