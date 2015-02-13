# Monte Carlo Statistical Simulation in Lattice model

python+swig c++ extension module
&&
matlab+mex c++ extension module  

## Prerequisite
=====================
* For ubuntu/debian
g++ 4.3.2
	sudo apt-get install g++
	http://gcc.gnu.org/
python 2.5
	sudo apt-get install python-dev
	http://docs.python.org/index.html
numpy 1.1.1
	sudo apt-get install python-numpy
	http://docs.scipy.org/doc/
swig  1.3.39
	sudo apt-get install swig
	http://www.swig.org/doc.html
matplotlib 0.98.3
	sudo apt-get install python-matplotlib
	http://matplotlib.sourceforge.net/
TkInter,Revision: 50704
	sudo apt-get install python-tk
	http://wiki.python.org/moin/TkInter
scons v0.98.5
	sudo apt-get install scons
	http://www.scons.org/doc/production/HTML/scons-user.html

* For Windows XP
Microsoft Visual Studio 2008
	http://msdn.microsoft.com/en-us/vstudio/default.aspx
python 2.5
swigwin-1.3.39
c:\dev\swigwin-1.3.39  for swig.exe
c:\Python26\Scripts  for scons.bat

## INSTALLATION
==========================
* for python-swig
cd <MCSS_SRC>/python-swig
scons

* for cpps
open Microsoft Visual Studio 2008
click 'Tools'->'Visual Studio 2008 Command Prompt'

* for matlab 
open Matlab2009
set Current Directory: <MCSS_SRC>/matlab-mex
setup

## RUN 
==========================
* for python
$ python wxMCSS.py

* for matlab 
>> check NVTIsing.m NVTPotts.m NVTLinear.m

## UNINSTALL
==========================
* for python
uninstall
scons --clean

* for matlab
uninstall
setup clean
mex -setup

	
	




