# MESA2SPH

## Introduction
From discreet MESA stellar density profile data, this package creates an Initial Conditions (IC) file readable by smoothed-particle hydrodynamics (SPH) codes such as GADGET-2.  
Generating MESA data requires basic operational understanding of MESA. Running GADGET-2, AREPO, or other SPH codes are the user's responsibility. 

The primary purpose of this software is to generate "NR" files which characterize the mass distribution in the outer regions of an arbitrary model star.

An NR file is a basic text file containing 4 columns of numerical data: 
(1) integer N (dimensionless),
(2) radii (cm),
(3) masses M (grams), and
(4) internal energy values as computed by MESA (cgs units).

The first column contains the integer used by HEALPix to determine the number of particles (np) distributed over a spherical shell--see the HEALPix papers for more information. 

The radii represent midpoint values rmid = (ru + rl) / 2 , in physical units (cm), corresponding to the radius at which a given shell must be placed in order to recast the mass contained in a region ru - rl as a shellular distribution of particles. The third column contains the mass contained in the computed interval. The fourth contains the internal energy (cgs units) at the midpoint radius.
The length of the file corresponds to the number of shells needed to reconstruct the sampled 1D density profile as a 3D particle distribution built of concentric shells. 


## Installation
From the directory in which this is unpacked, run 

	./setup.py local

Further information on installation is described in the 1D-MESA2SPH-3D user's guide.

## Prerequisites
This project was interfaced with mesa version 10398 with example inlists and input models (e.g. wd.mod) sourced from data provided within mesa-r10398. See 1D-MESA2SPH-3D user's guide.


## System requirements:
	mesa-r10398
		-requires latest version of the mesasdk
	HEALPix numerical libraries	
	hdf5 libraries 


	See 1D-MESA2SPH-3D user's guide for a list of Python dependencies


## Operation
Basic operation proceeds from the "work" subdirectory via

	./run test.cfg

See 1D-MESA2SPH-3D user's guide for more detailed information.