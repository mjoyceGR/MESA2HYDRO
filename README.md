# MESA2GADGET

## Introduction
From a MESA stellar density profile, creates an Initial Conditions (IC) file readable by sph codes such as GADGET2  
Generating MESA data requires basic operational understanding of MESA

The primary purpose of this software is to generate "NR" files corresponding to the outer regions of an arbitrary model star.
An NR file is a basic text file containing 3 columns of numerical data. The first contains an integer N (always repeating), the second contains radii (in cm), and the third contains a mass M (in grams). The first column corresponds to the integer used by HEALPix to determine the number of particles placed on a spherical shell, where n_particles=12 N^2. This is a constraint imposed by HEALPix, see ref{} for more information.
The radii represent midpoint values rmind = (r_upper + r_lower)/ 2, in physical units (cm), corresponding to the radius at which a given shell must be placed in order to recast the mass contained in a region (r_upper - r_lower) as a shellular distribution of particles. The third column contains the mass contained in the computed interval. The length of the file corresponds to the number of shells needed to reconstruct the sampled 1D density profile as a 3D particle distribution built of concentric shells. 

## Installation
???

## Prerequisites
This project was interfaced with mesa version 10398 with example inlists and input models (e.g. wd.mod) sourced from data provided within mesa-r10398
You must set the environment variable 'MESA2GADGET_ROOT' in your bashrc


## System requirements:
	mesa-r10398
		-requires latest version of the mesasdk
	healpy
		-contains necessarily HEALPix libs?
	hdf5 if you hate yourself

	The following python libraries are required:



## Operation
You must set the following environment variables in your bashrc:

'MESASDK_ROOT'

'MESA2GADGET_ROOT'

'MESA_DIR'

e.g. 

	export MESASDK_ROOT=/home/usr/mesasdk

	export MESA_DIR=/home/usr/mesa-r10398

	export OP_NUM_THREADS=4 
	
	export MESA2GADGET_ROOT=/home/usr/MESA2GADGET

and add this source command to your .bashrc to initialize the MESA software development kit
	source $MESASDK_ROOT/bin/mesasdk_init.sh