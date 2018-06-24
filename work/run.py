#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
MESA_PKG_DIR = os.path.abspath(
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '..'))
sys.path.insert(0, os.path.join(MESA_PKG_DIR, 'lib'))
import MESAhandling as MJ
import converge_funcs as cf
import io_lib as rw
import mainlib as mn
from cfg_parser import *

import time
start_time = time.time()


#######################################################
#
# LOAD PARAMETER VALUES
#
########################################################

## To add new parameters add name and default value
SCRIPT_CONFIGS = {
    'check_MESA_profile': True,
    'MESA_file': 'data/sample_MESA_output/profile_mainsequence_logE.data',
    'make_NR_file': False,
    'loaded_NR_filename': 'work/NR_files/saveNR_ms.dat',
    'new_NR_filename': 'latest_NR.dat',

    'make_IC_file': False,
    'loaded_IC_filename': 'ms_logE',
    'new_IC_filename': 'latest_IC',
    'IC_format_type': 'binary',

    'masscut': 0.95,
    'N': 8,
    'mp': 1e-7,
    'mp_cgs': 1.988e26,
    'stepsize': 2.45e6,

    'try_reload': False,
    'png_tag': 'latest',

    'reload_bin_number': 70.0,
    'use_bins_from_NR_file': False,
    'which_dtype':'f'}
 

# Auto generates command line arguments and configuration file reading from SCRIPT_CONFIG definition
options = OptionInputs(SCRIPT_CONFIGS, description='Program for converting MESA output to Gadget simulation files')

# Only arguments in SCRIPT_CONFIGS can be found in user_configs
user_configs = options.get_configs()


check_MESA_profile = user_configs['check_MESA_profile']
make_NR_file = user_configs['make_NR_file']
make_IC_file = user_configs['make_IC_file']
try_reload = user_configs['try_reload']
IC_format_type = user_configs['IC_format_type']
MESA_file = relative_to_root(user_configs['MESA_file'])
masscut = user_configs['masscut']
N = user_configs['N']
mp = user_configs['mp']
mp_cgs = user_configs['mp_cgs']
# doesn't seem to be used anymore
#startype = user_configs['startype']
stepsize = user_configs['stepsize']
png_tag = user_configs['png_tag']
new_NR_filename=user_configs['new_NR_filename']
loaded_NR_filename = user_configs['loaded_NR_filename']
new_IC_filename = user_configs['new_IC_filename']
loaded_IC_filename = user_configs['loaded_IC_filename']

## meridith additional parameters
reload_bin_number = user_configs['reload_bin_number']
use_bins_from_NR_file=user_configs['use_bins_from_NR_file']
which_dtype=user_configs['which_dtype']



if make_NR_file:
    nrfile=new_NR_filename #loaded_NR_filename
else:
    nrfile=loaded_NR_filename
if make_IC_file:
    icfile=new_IC_filename
else:        
    icfile=loaded_IC_filename

# makes the path relative to MESA2GADGET_ROOT directory
nrfile = relative_to_root(nrfile)
icfile = relative_to_root(icfile)

############################################################
###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################
mp=mp*M_to_solar
#try:
mp_cgs=mp_cgs    

if mp != mp_cgs:
    print '\n\nWARNING! Inconsistent values of mp and mp_cgs!'
    print 'using mp'


##############################################################
#
# make sure MESA density profile being loaded isn't nonsense
#
#############################################################
if check_MESA_profile:
    print 'Confirm sensible profile data...'
    mn.check_MESA(MESA_file, masscut,uselog=True, save=True)


##############################################################
#
# Generate shell placement radii
#
#############################################################
### need to include internal energy here
if make_NR_file:
    t1=time.time()
    print '\n\nGenerating NR file...'
    if os.path.exists(nrfile):
        print("NR file already exists. Cannot overwrite {}.".format(nrfile))
        sys.exit(1)
    outf=open(nrfile,"w")
    mn.make_NR_file(MESA_file,masscut,N,mp, stepsize,outf)
    outf.close()
    print 'NR file generation complete!'
    print("--- %s seconds ---" % (t1 - start_time))


##############################################################
#
# Convert placement radii into HEALPix shells and send to IC
#
#############################################################
#which_dtype='d'
if make_IC_file:

    t2=time.time()
    print '\n\nGenerating IC file...'
    in_file=nrfile
    out_file=icfile
    mn.get_IC(in_file,out_file,mp,format_type=IC_format_type,which_dtype=which_dtype)
    print 'IC file generation complete!'
    print("--- %s seconds ---" % (time.time() - t2))


###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################
if try_reload:
    r_recovered, masses_recovered=mn.reload_IC(icfile,IC_format_type, which_dtype=which_dtype)
    print "min and max radius from recovered data: ", r_recovered.min(), r_recovered.max()

    if use_bins_from_NR_file:
        r_reloaded,rho_reloaded=mn.bins_from_NR(nrfile,r_recovered,masses_recovered[0])

    else:
        nbin=reload_bin_number
        print "plotting data using ", reload_bin_number, " bins"
        r_reloaded,rho_reloaded=mn.binned_r_rho(r_recovered, masses_recovered[0], reload_bin_number)


    mn.quick_plot(MESA_file,masscut, r_reloaded,rho_reloaded,IC_format_type,png_tag=png_tag)   


print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
