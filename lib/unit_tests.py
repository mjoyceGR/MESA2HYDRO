#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#import pygadgetreader as pgr # works- credit this person
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/lib/')
import MESA2GADGET.mesalib.MESAhandling as MJ
import datetime as dt 
import random as rand
import healpy as hp
import MESA2GADGET.mesalib.converge_funcs as cf  



###################################################
#
# Runge Kutta
#
##################################################
def f(x,placeholder):
	y=1.0/x**2.0
	return y
step=0.001
x=1
x_upper=3
y=0
while x<x_upper:
	x,y=cf.RK1(x,y,9999,f,step)
print "obtained value:", y, "\ncorrect value: ", 2.0/3.0


###################################################
#
# Rotate Arrays
#
##################################################


