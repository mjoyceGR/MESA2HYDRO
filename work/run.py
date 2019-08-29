#!/usr/bin/env python

#************************************************************************
#
#   Copyright (C) 2019  M. Joyce, L. Lairmore, D. J. Price
#
#   See MESA2HYDRO/LICENSE
#
#************************************************************************

'''

program execution script

'''


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
    'IC_format_type': 'phantom_binary',

    'masscut': 0.95,
    'r_depth':0.5,

    'N': 8,
    'mp': 1e-7,
    'mp_cgs': 1.988e26,
    'stepsize': 2.45e6,
    
    'which_dtype':'f',
    'TOL':0.01}
 
# Auto generates command line arguments and configuration file reading from SCRIPT_CONFIG definition
options = OptionInputs(SCRIPT_CONFIGS, description='Program for converting MESA output to SPH-compatible IC files')

# Only arguments in SCRIPT_CONFIGS can be found in user_configs
user_configs = options.get_configs()

check_MESA_profile = user_configs['check_MESA_profile']
make_NR_file = user_configs['make_NR_file']
make_IC_file = user_configs['make_IC_file']
IC_format_type = user_configs['IC_format_type']
MESA_file = relative_to_root(user_configs['MESA_file'])
masscut = user_configs['masscut']
r_depth = user_configs['r_depth']

N = user_configs['N']
mp = user_configs['mp']
mp_cgs = user_configs['mp_cgs']
stepsize = user_configs['stepsize']
new_NR_filename=user_configs['new_NR_filename']
loaded_NR_filename = user_configs['loaded_NR_filename']
new_IC_filename = user_configs['new_IC_filename']
loaded_IC_filename = user_configs['loaded_IC_filename']

which_dtype=user_configs['which_dtype']
TOL = user_configs['TOL']

## add 'lognorm' boolean here 7/30/19



if make_NR_file:
    nrfile=new_NR_filename 
else:
    nrfile=loaded_NR_filename

if make_IC_file:
    icfile=new_IC_filename
else:        
    icfile=loaded_IC_filename

# makes the path relative to MESA2HYDRO_ROOT directory
nrfile = relative_to_root(nrfile)
icfile = relative_to_root(icfile)



###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################



##################################################################
#
# convert radial depth to masscut
#
# if a radial preference has been passed, override mass depth and
# define mass proportion relative to r_depth
#
# set r_depth to -1.0 in config file to rely on mass depth
#
###################################################################
if r_depth != -1.0:
    fit_region_R = 10.0**(MJ.get_quantity(MESA_file, 'logR').astype(np.float)) #mn.MESA_r(MESA_file, 0) #whole star
    fit_region_M = MJ.get_quantity(MESA_file,'mass').astype(np.float)*M_to_solar  #mn.MESA_m(MESA_file, 0)
    
    
    print "len(R), len(M): ", len(fit_region_R), len(fit_region_M)


    r_bound = r_depth*fit_region_R.max()

    assoc_region=np.where( fit_region_R >= r_bound)[0]
    print "len(assoc_region)", len(assoc_region)
    mass_assoc = fit_region_M[assoc_region]

    masscut = mass_assoc.min()/mass_assoc.max()
    print "estimated masscut from radial input = ", masscut


else:
    masscut = masscut
print "masscut being used is: ", masscut


mp=mp*M_to_solar
mp_cgs=mp_cgs    

if mp != mp_cgs:
    print '\nWARNING! Inconsistent values of mp and mp_cgs! Choosing mp'


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
Romberg=False#True

if make_NR_file:
    t1=time.time()
    print '\n\nGenerating NR file...'
    if os.path.exists(nrfile):
        print("NR file already exists. Cannot overwrite {}.".format(nrfile))
        sys.exit(1)
    outf=open(nrfile,"w")
    mn.make_NR_file(MESA_file,masscut,N,mp, stepsize, TOL, outf, Romberg=Romberg)
    outf.close()
    print 'NR file generation complete!'
    print("--- %s seconds ---" % (t1 - start_time))



##############################################################
#
# Convert placement radii into HEALPix shells and send to IC
#
#############################################################
if make_IC_file:
    t2=time.time()
    print '\n\nGenerating IC file...'
    in_file=nrfile
    out_file=icfile
    if IC_format_type=="phantom_binary":
        print "setting lognorm..."
        #sys.exit()
        lognorm= True ## make this a passable parameter
        mn.get_IC(MESA_file, masscut, in_file, out_file, mp,\
                format_type=IC_format_type,which_dtype=which_dtype, lognorm=lognorm)
        time.sleep(2)

    else:
        mn.get_IC(MESA_file, masscut, in_file, out_file, mp,\
                format_type=IC_format_type,which_dtype=which_dtype)
    print 'IC file generation complete!'
    print("--- %s seconds ---" % (time.time() - t2))


###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################

# ########## TEMPORARY ###############
# IC_format_type=
# ########## TEMPORARY ###############

# if try_reload:
#     r_recovered, masses_recovered=mn.reload_IC(icfile,IC_format_type, which_dtype=which_dtype)
#     print "min and max radius from recovered data: ", r_recovered.min(), r_recovered.max()

#     if use_bins_from_NR_file:
#         r_reloaded,rho_reloaded=mn.bins_from_NR(nrfile,r_recovered,masses_recovered[0])

#     else:
#         nbin=reload_bin_number
#         print "plotting data using ", reload_bin_number, " bins"
#         r_reloaded,rho_reloaded=mn.binned_r_rho(r_recovered, masses_recovered[0], reload_bin_number)


#     mn.quick_plot(MESA_file,masscut, r_reloaded,rho_reloaded,IC_format_type,png_tag=png_tag)   


print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))

# end run script