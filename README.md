# MESA2HYDRO


## UPDATE
MESA2HYDRO is now written in Python3. The most up-to-date documentation is here on GitHub.

The cython package is no longer a dependency!  

## Documentation
Read the paper

	Refereed: https://iopscience.iop.org/article/10.3847/1538-4357/ab3405/meta

	Free: https://arxiv.org/abs/1907.09062

Read the User's Guide (WARNING: updated less frequently):

	https://github.com/mjoyceGR/MESA2HYDRO/blob/master/DOCUMENTATION/MESA2HYDRO_users_guide.pdf

## Introduction
From discrete stellar structure data in the style of MESA (http://mesa.sourceforge.net/) density profiles, this package creates initial conditions (IC) files readable by smoothed-particle hydrodynamics (SPH) codes such as Phantom (https://phantomsph.bitbucket.io/).  

The primary purpose of this software is to generate "NR" files, which contain coordinates parameterizing the mass distribution of an arbitrary model star to arbitrary depth.

An NR file is a basic text file containing 4 columns of numerical data: 
(1) integer N (dimensionless),
(2) radii (cm),
(3) masses M (grams), and
(4) internal energy (ergs).

The first column contains the integer used by HEALPix to determine the number of particles (np) distributed over a spherical shell--see the HEALPix papers for more information. 

The radii represent midpoint values rmid = (ru + rl) / 2 , in physical units (cm), corresponding to the radius at which a given shell must be placed in order to recast the mass contained in a region ru - rl as a shellular distribution of particles. The third column contains the mass contained in the computed interval. The fourth contains the internal energy (cgs units) at the midpoint radius.
The length of the file corresponds to the number of shells needed to reconstruct the sampled 1D density profile as a 3D particle distribution built of concentric shells. 


## Recommended Installation: Git/tarball
Cloning this git repository (or downloading the tarball) is the safest way to collect all software and dependencies: 

	git clone https://github.com/mjoyceGR/MESA2HYDRO.git

Download and upack the .zip file via

	unzip MESA2HYDRO-master.zip

Download and upack a .tar file via

	tar -xvf MESA2HYDRO-master.tar

After cloning or unpacking, move to the head directory via

	cd MESA2HYDRO/

Issue the following command in the top level directory

	python3 local_setup.py

to install the package locally.


To compile the Fortran components separately before installation, run


	./setup.py build_ext

followed by

	./setup.py install --user



## Installation via pip

MESA2HYDRO is also theoretically available via pip, but this will not be updated as frequently as the git repo. 
WARNING: some issues may occur if you have multiple versions of pip installed.

Using pip3, run

	pip3 install MESA2HYDRO


In some cases, "pip" and "python" can/should replace "pip3" and "python3", respectively (i.e. if your default version is Python3)

To upgrade, run

	pip install --upgrade MESA2HYDRO==1.1.XXXX

where 0.1.XXXX is the lastest version number on https://pypi.org/manage/project/MESA2HYDRO/releases/

Sometimes the error 

	ERROR: Cannot uninstall 'MESA2HYDRO'. It is a distutils installed project and thus we cannot accurately determine which files belong to it which would lead to only a partial uninstall.

is triggered by an attempt to use pip install --upgrade. As a workaround, remove any current installation of MESA2HYDRO, followed by 

	pip install MESA2HYDRO==0.1.XXXX

or 

	pip install MESA2HYDRO

for the most recent version.


If you've failed to specify pip3 and have other versions of Python installed, the package may end up somewhere unexpected.
You can locate the package by typing

	pip show MESA2HYDRO

You may want to copy or move the entire MESA2HYDRO root directory somewhere where it is easier to work out of directly 

	mv MESA2HYDRO/ /home/your_name/

From there, you should set the $MESA2HYDRO_ROOT environment variable in your .bashrc to point to the root directory:

	export MESA2HYDRO_ROOT=/home/your_name/MESA2HYDRO

Then

	cd MESA2HYDRO/work/

and follow the instructions for Running a Test Case. No additional setup or compilation should be required if installed via pip. 

## Prerequisites
For a complete list of Python and external software dependencies, see the Prerequisites section of the user's guide.
This project was developed with MESA version 10398 and uses example inlists and input models (e.g. wd.mod) included in mesa-r10398. Neither Phantom nor MESA is required for operation, but we recommend these two packages. 


## Running a Test Case
To run a basic MESA2HYDRO calculation, move to the work directory via

	cd MESA2HYDRO/work/

Basic operation proceeds via

	python3 run_conversion.py mainsequence.cfg

where mainsequence.cfg is a configuration file. These are the files that the user should modify to suit their problem. 



## Troubleshooting
Failure may be caused by missing packages. The user must have, at minimum, the following Python3 and external packages installed on their machine:

	healpy
	matplotlib
	numpy
  	python-tk
  	scipy.interpolate
	scipy.optimize


If they were not automatically installed via pip or the setup procedure, these can be installed from the command line via, e.g., 

	sudo apt-get install python3-numpy

or similarly via, e.g., 

	pip3 install numpy 

The authors have had better results installing these additional libraries via command line rather than through pip. 

In particular, one must have numerical and other libraries required by healpy/HEALPix (https://healpy.readthedocs.io/en/latest/) installed on their machine. 
More information on dealing with healpy can be found at https://healpy.readthedocs.io/en/latest/install.html


Some failures to install can be rectified by installing python3-dev via 
	
	sudo apt-get install python3-dev

In your Python 3 environment is out of date, one can install the required dependencies via 

	sudo apt-get install python3-matplotlib
	sudo apt-get install python3-scipy
	sudo apt-get install python3-h5py
	sudo apt-get install python3-healpy

or equivalent. See "requirements.txt" in the top MESA2HYDRO directory for all dependencies