# MESA2HYDRO

## Documentation
Read the paper

	Refereed: https://iopscience.iop.org/article/10.3847/1538-4357/ab3405/meta

	Free: https://arxiv.org/abs/1907.09062

Read the User's Guide:

	https://github.com/mjoyceGR/MESA2HYDRO/blob/master/DOCUMENTATION/MESA2HYDRO_users_guide.pdf

## Introduction
From discrete stellar structure data in the style of MESA (http://mesa.sourceforge.net/) density profiles, this package creates initial conditions (IC) files readable by smoothed-particle hydrodynamics (SPH) codes such as Phantom (https://phantomsph.bitbucket.io/).  

The primary purpose of this software is to generate "NR" files that parameterize the mass distribution in an arbitrary model star.

An NR file is a basic text file containing 4 columns of numerical data: 
(1) integer N (dimensionless),
(2) radii (cm),
(3) masses M (grams), and
(4) internal energy (ergs).

The first column contains the integer used by HEALPix to determine the number of particles (np) distributed over a spherical shell--see the HEALPix papers for more information. 

The radii represent midpoint values rmid = (ru + rl) / 2 , in physical units (cm), corresponding to the radius at which a given shell must be placed in order to recast the mass contained in a region ru - rl as a shellular distribution of particles. The third column contains the mass contained in the computed interval. The fourth contains the internal energy (cgs units) at the midpoint radius.
The length of the file corresponds to the number of shells needed to reconstruct the sampled 1D density profile as a 3D particle distribution built of concentric shells. 


## Installation
Via pip:

	pip install MESA2HYDRO

Alternatively, clone the git repository, unpack, and run the following from the top level directory: 

	python setup.py install

In some cases, "pip2" and "python2" should replace "pip" and "python", respectively (i.e. if your default version is Python3--we will be releasing this package for Python3 shortly)

To upgrade, run

	pip install --upgrade MESA2HYDRO==0.1.24

where 0.1.24 is the lastest version number on https://pypi.org/manage/project/MESA2HYDRO/releases/


## Running a test case
To successfully run a basic MESA2HYDRO instance

	./run.py mainsequence_test.cfg

the user must have the following Python and external packages installed on their machine:

 	argparse
	cython
	healpy
	matplotlib.pyplot
	numpy
  	python-tk
	random
  	scipy.interpolate
	scipy.optimize
	time

These can be installed from the command line via, e.g., 

	sudo apt-get install python-tk

or similarly via, e.g., 

	pip install numpy 

The authors have had better results installing these additional libraries via command line rather than through pip. 

In particular, one must have numerical and other libraries required by healpy/HEALPix (https://healpy.readthedocs.io/en/latest/) installed on their machine. The installation command for healpy is:

	pip install healpy

More information on dealing with healpy can be found at https://healpy.readthedocs.io/en/latest/install.html


## Prerequisites
For a complete list of Python and external software dependencies, see the Prerequisites section of the user's guide.

This project was designed to interface with MESA version 10398, using example inlists and input models (e.g. wd.mod) sourced from data provided within mesa-r10398. See 1D-MESA2HYDRO-3D user's guide.


## Operation
Basic operation proceeds from the "work" subdirectory via

	python2 run.py <filename>.cfg

See 1D-MESA2HYDRO-3D user's guide for more detailed information.

## Common Installation Issues
Some failures to install can be rectified by installing python-dev via 
	
	sudo apt-get install python-dev
