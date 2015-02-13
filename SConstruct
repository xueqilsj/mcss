import os
import sys

import glob

version=1.3 #oct. 27,2009
version=1.4 #oct. 31,2009, unsigned int
version=1.5 #Nov. 3,2009, Temp OpenMP SEO
version=1.6 #Nov. 13,2009, Signle par,Heat bath ising
version=1.7 #Nov. 17,2009, GRCE extend Wollf cluster,fix Delt():Q-1
version=1.8 #Nov. 20,2009, GCE extend Wollf cluster,format input
version=1.9 #Nov. 20,2009, wxPython GUI,LJ-Fluid

env = Environment()   
env.Append(ENV=os.environ)
   	
SRC=[glob.glob('./core-cpp/*.h'),glob.glob('./core-cpp/*.cpp'),glob.glob('./pure-cpp/*.h'),glob.glob('./pure-cpp/*.cpp'),'./pure-cpp/SConsMPI','./pure-cpp/SConsOpenMP','./pure-cpp/SConstruct',\
glob.glob('./matlab-mex/*.cpp'),glob.glob('./matlab-mex/*.m'),glob.glob('./python-swig/*.py'),glob.glob('./python-swig/*.i'),'./python-swig/SConstruct','README','INSTALL','SConstruct']


SRC_pure_cpp=[glob.glob('./core-cpp/*.h'),glob.glob('./core-cpp/*.cpp'),glob.glob('./pure-cpp/*.h'),glob.glob('./pure-cpp/*.cpp'),'./pure-cpp/SConsMPI','./pure-cpp/SConsOpenMP','./pure-cpp/SConstruct']

PACKAGE='mcss%g-src'%version
PACKAGE_pure_cpp='mcss%gpure-cpp-src'%version
	
if sys.platform[:3]=='win':
	f=env.Zip('dist'+os.sep+PACKAGE+'.zip',SRC)
	Clean(f, 'dist')
else:
	env=Environment(TARFLAGS='-c -z')
	f=env.Tar('dist'+os.sep+PACKAGE+'.tar.gz',SRC)
	#f=env.Tar('dist'+os.sep+PACKAGE_pure_cpp+'.tar.gz',SRC_pure_cpp)
	Clean(f, 'dist')
	