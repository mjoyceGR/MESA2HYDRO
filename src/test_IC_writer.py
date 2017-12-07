#!/usr/bin/env python
import numpy as np
import os
import sys
import math
import matplotlib.pyplot as plt
import hdf5lib ## 9/5/17

import read_write_HDF5 as rw
import h5py
import MESAlibjoyce as MJ

######################################
#
# (4) write the Gadget/GIZMO initial conditions file:
# 	see write_GIZMO_ICs.py
#
######################################
out_fname='box_3d_r32.hdf5'
var=MJ.make_IC_box_hdf5(out_fname)
print var
print type(var)