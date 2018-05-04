#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/lib/')
try:
    import MESA2GADGET.lib.MESAlibjoyce as MJ
    import MESA2GADGET.lib.converge_funcs as cf
    import MESA2GADGET.lib.io_lib as rw
    import MESA2GADGET.lib.mainlib as mn
    from MESA2GADGET.lib.cfg_parser import *

    from MESA2GADGET import MESA_PKG_DIR
except ImportError as err:
    print("Err: {}".format(err))
    print("Problem with MESA2GADGET installation")
    print("To use this package please run sudo python setup.py install")
    print("or set your PYTHONPATH environment variable to the directory")
    print("MESA2GADGET is in (pointing it directly to MESA2GADGET still causes problems)")
    exit(1)

import time
start_time = time.time()


#######################################################
#
# LOAD PARAMETER VALUES
#
########################################################
VALID_CONFIGS = {
    'check_MESA_profile': True,
    'MESA_file': path_from_package('out/sample_MESA_output/profile_mainsequence_logE.data'),
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
    

parser=command_line_parser()
args = parser.parse_args()
if args.defaults:
    for config, value in VALID_CONFIGS.items():
        print("{} = {}".format(config, value))
    exit(0)


###################### THIS IS THE MAIN THING I WANT TO DEAL WITH #######################
elif args.config_file:
    if not os.path.exists(args.config_file):
        args.config_file = path_from_package(args.config_file)

    if not os.path.exists(args.config_file):
        print("Config file does not exist")

    with open(args.config_file, 'r') as f:
        config_file = f.read()

    lines = config_file.splitlines()
    lines = [line for line in lines
             if not re.match('\s*#', line) and '=' in line]
    config_file = "\n".join(lines)

    check_MESA_profile = get_bool_option(config_file,'check_MESA_profile', VALID_CONFIGS['check_MESA_profile'])
    make_NR_file = get_bool_option(config_file,'make_NR_file', VALID_CONFIGS['make_NR_file'])
    make_IC_file = get_bool_option(config_file,'make_IC_file', VALID_CONFIGS['make_IC_file'])
    try_reload = get_bool_option(config_file,'try_reload', VALID_CONFIGS['try_reload'])
    IC_format_type = get_str_option(config_file,'IC_format_type', VALID_CONFIGS['IC_format_type'])
    MESA_file = get_path_option(config_file,'MESA_file', VALID_CONFIGS['MESA_file'])
    masscut = get_float_option(config_file,'masscut', VALID_CONFIGS['masscut'])
    N = get_int_option(config_file,'N', VALID_CONFIGS['N'])
    mp = get_float_option(config_file,'mp', VALID_CONFIGS['mp']) ##IN UNITES OF Msolar!!!
    mp_cgs = get_float_option(config_file,'mp_cgs', VALID_CONFIGS['mp_cgs'])
    stepsize = get_float_option(config_file,'stepsize', VALID_CONFIGS['stepsize']) ##IN UNITES OF cm
    png_tag = get_str_option(config_file,'png_tag',VALID_CONFIGS['png_tag'])

    if make_NR_file:
        new_NR_filename = get_path_option(config_file,'new_NR_filename',VALID_CONFIGS['new_NR_filename'])
    else:
        loaded_NR_filename = get_path_option(config_file,'loaded_NR_filename', VALID_CONFIGS['loaded_NR_filename'])
    if make_IC_file:
        new_IC_filename = get_path_option(config_file,'new_IC_filename',VALID_CONFIGS['new_IC_filename'])
    else:
        loaded_IC_filename = get_path_option(config_file,'loaded_IC_filename', VALID_CONFIGS['loaded_IC_filename'])

    ## Meridith adding more parameters
    reload_bin_size = get_float_option(config_file,'reload_bin_size',VALID_CONFIGS['reload_bin_size'])
    use_bins_from_NR_file= get_bool_option(config_file,'use_bins_from_NR_file', VALID_CONFIGS['use_bins_from_NR_file'])
    which_dtype= get_str_option(config_file,'which_dtype', VALID_CONFIGS['which_dtype'])


else:
    check_MESA_profile = args.check_MESA_profile
    make_NR_file = args.make_NR_file
    make_IC_file = args.make_IC_file
    try_reload = args.try_reload
    IC_format_type = args.IC_format_type
    MESA_file = args.MESA_file
    masscut = args.masscut
    N = args.N
    mp = args.mp
    mp_cgs = args.mp_cgs
    startype = args.star_type
    stepsize = args.step_size
    png_tag = args.png_tag
    if args.make_NR_file:
        nrfile = args.new_NR_filename
    else: 
        nrfile = args.loaded_NR_filename

    ## meridith additional parameters
    reload_bin_size = args.reload_bin_size
    use_bins_from_NR_file=args.use_bins_from_NR_file
    which_dtype=args.which_dtype



if make_NR_file:
    nrfile=new_NR_filename#loaded_NR_filename
else:
    nrfile=loaded_NR_filename
if make_IC_file:
    icfile=new_IC_filename
else:        
    icfile=loaded_IC_filename
# If the given path exists use that otherwise prepend the package path
if not os.path.exists(nrfile):
    nrfile = path_from_package(nrfile)
if not os.path.exists(icfile):
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

    else:
        nbin=reload_bin_size
        print "plotting data using binsize=", reload_bin_size
        r_reloaded,rho_reloaded=mn.binned_r_rho(r_recovered, masses_recovered[0], reload_bin_size)


    mn.quick_plot(MESA_file,masscut, r_reloaded,rho_reloaded,IC_format_type,png_tag=png_tag)   



print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
