#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
sys.path.append(os.path.abspath('../lib'))
import MESAlibjoyce as MJ
import converge_funcs as cf
import io_lib as rw
import mainlib as mn
from cfg_parser import *

MESA_PKG_DIR = os.path.abspath(os.path.dirname(__file__))

m2g_save_path=os.path.abspath(os.environ.get('MESA2GADGET_ROOT'))

if m2g_save_path is None:
    print("Environmental variable 'MESA2GADGET_ROOT' isn't set")
    print("Storing output data in default directory root {}"
          .format(MESA_PKG_DIR))
    m2g_save_path = MESA_PKG_DIR
else:
    print("Files will be saved in {}\n\n".format(m2g_save_path))
    if not os.path.exists(m2g_save_path):
        os.makedirs(m2g_save_path)

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
    'MESA_file': os.path.join(MESA_PKG_DIR, 'out/sample_MESA_output/profile_mainsequence_logE.data'),
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

    'reload_bin_size': 70.0,
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
MESA_file = user_configs['MESA_file']
masscut = user_configs['masscut']
N = user_configs['N']
mp = user_configs['mp']
mp_cgs = user_configs['mp_cgs']
# doesn't seem to be used anymore
#startype = user_configs['startype']
stepsize = user_configs['stepsize']
png_tag = user_configs['png_tag']
loaded_NR_filename = user_configs['loaded_NR_filename']
new_IC_filename = user_configs['new_IC_filename']
loaded_IC_filename = user_configs['loaded_IC_filename']

## meridith additional parameters
reload_bin_size = user_configs['reload_bin_size']
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

# If the given path exists use that otherwise prepend the package path
if not os.path.exists(nrfile):
    if make_NR_file:
        nrfile = os.path.join(m2g_save_path, nrfile)
    else:
        nrfile = path_from_package(nrfile)
if not os.path.exists(icfile):
    if make_IC_file:
        icfile = os.path.join(m2g_save_path, icfile)
    else:
        icfile = path_from_package(icfile)



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

    print "min and max radius from r recovered data: ", r_recovered.min(), r_recovered.max()

    #print "\n\ndatatype found before loading binned data:   which_dtype=" , which_dtype, '\n\n'


    if use_bins_from_NR_file:
        print 'using radial binning from NR file'
        r_reloaded,rho_reloaded=mn.bins_from_NR(nrfile,r_recovered,masses_recovered[0])

        r_norm=r_reloaded/r_reloaded.min()
        rho_norm=rho_reloaded/rho_reloaded.min()


        A_guess=35.0#50
        B_guess=4.0##44
        C_guess=8#0
        D_guess=4.8#45

        try:
            A,B,C,D=cf.get_curve(r_norm,rho_norm,A_guess,B_guess,C_guess,D_guess, cf.one_over_r)
        except RuntimeError:
            # A=1.0
            # B=1.0
            # C=1.0
            # D=1.0
            A=A_guess
            B=B_guess
            C=C_guess
            D=D_guess

        print "found A,B,C,D:", A,B,C,D
        rho_fit_found=cf.one_over_r(r_norm,A,B,C,D)

        A=A_guess
        B=B_guess
        C=C_guess
        D=D_guess
        #D=8.2e8
        rho_fit_forced=cf.one_over_r(r_norm,A,B,C,D)#,D)
        # #print rho_fit

        plt.plot(r_norm,rho_norm,'yo', color='orange', label='loaded data')
        plt.plot(r_norm,rho_fit_found, 'r-', linewidth=3, label=r'fit to rho w/ optimize')#r_reloaded
        plt.plot(r_norm,rho_fit_forced, 'g-', linewidth=3, label=r'fit to rho by hand')#r_reloaded
        plt.ylim(-1,55)
        plt.legend(loc=1)
        plt.savefig('fit.png')
        #plt.show()
        plt.close()


    else:
        nbin=reload_bin_size
        print "plotting data using binsize=", reload_bin_size
        r_reloaded,rho_reloaded=mn.binned_r_rho(r_recovered, masses_recovered[0], reload_bin_size)


    #mn.quick_plot(MESA_file,masscut, r_reloaded,rho_reloaded,IC_format_type,png_tag=png_tag)   



print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
