# MESA2GADGET

## Introduction
Creates Gadget IC file from a MESA stellar density profile 

## Installation
Install every numerical python library known to man and then have them battle to the death.

## Example
Put in numbers and you will get numbers you should look at

## Prerequisites
This project was interfaced with mesa version 10398 with example inlists and input models (e.g. wd.mod) sourced from those provided within mesa-r10398

System requirements:
	mesa-r1098
		requires latest version of the mesasdk


## Operation
You must set the following environment variables in your bashrc:

'MESASDK_ROOT'

'MESA2GADGET_ROOT'

'MESA_DIR'

e.g. 

	export MESASDK_ROOT=/home/meridith/mesasdk

	export MESA_DIR=/home/meridith/MESA/mesa-r10398

	export OP_NUM_THREADS=4 
	
	export MESA2GADGET_ROOT=/home/meridith/UCT_SAAO/detached_shells/MESA2GADGET

and add this source command to your .bashrc to initialize the MESA software development kit
	source $MESASDK_ROOT/bin/mesasdk_init.sh