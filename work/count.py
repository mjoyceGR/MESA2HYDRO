#!/usr/bin/env python
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import sys 
import converge_funcs as cf

#
import io_lib as rw
import mainlib as mainlib

IC_file="IC_ms_cgs_bad_hsml" #"IC_ms_90"
format_type='gadget_binary'

r, m = mainlib.reload_IC( IC_file, format_type, which_dtype='f')

print "total num particles: ", len(m)




f=open(IC_file+'.bin','r')
ptype=0
header=rw.load_gadget_binary_header(f)
attribute_dictionary=rw.load_gadget_binary_particledat(f, header, ptype, skip_bh=0, which_dtype='f')

positions=attribute_dictionary['Coordinates']
masses=attribute_dictionary['Masses'] ##all of these are the same value, mp
x=positions[:,0]
y=positions[:,1]
z=positions[:,2]


######## FIGURE OUT WHAT THIS KWARG IS
hsml=attribute_dictionary['HSML']
print hsml