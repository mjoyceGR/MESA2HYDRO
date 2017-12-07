#!/usr/bin/env python
import numpy as np
import os
import sys
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
out_fname='test.hdf5'

#function which takes N, gets healpix data, appends to mega file

# 10 18.0230974563
# 11 23.971146512
# 13 31.3128187525
# 15 40.35578738
# 16 51.6825885374
# 17 66.0363797552
# 18 84.291397913
# 18 108.998058194
# 14 149.5410294
# 0 176.018543271


NSIDE=16
r_mid=

x,y,z=MJ.get_coords(NSIDE,r_mid)

for i in range(len(x)):
	print x[i]**2.0 + y[i]**2.0  + z[i]**2.0



healpix_file='healpix_to_gadget_shell_1.dat'

coord_file='combined.dat'

mp=1e-8

var=MJ.make_IC_box_hdf5(out_fname,mp, coord_file)
print var
print type(var)



