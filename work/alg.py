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

# ###################################################
# M_to_solar=1.988*10.0**33.0 ## g/Msolar
# R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
# ###################################################

def path_from_package(path):
    return os.path.join(MESA_PKG_DIR, path)

#######################################################
#
# user controls
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
    'png_tag': 'latest'}
    
##############################################################
#
# Argument input
#
##############################################################

parser = argparse.ArgumentParser(description='Program for converting MESA output to Gadget simulation files')
parser.add_argument('--config-file',
                    help='Path to configuration file')
parser.add_argument('--defaults', action='store_true',
                    help='Lists the scripts default parameters and exits')
config_args = parser.add_argument_group("Configuration")

# Boolean option flags set opposite of file default
config_args.add_argument('--check-MESA',
                         action='store_true' if VALID_CONFIGS['check_MESA_profile'] else 'store_false',
                         help='Sets check-MESA value to {}'
                              .format(not VALID_CONFIGS['check_MESA_profile']))
config_args.add_argument('--make-NR-file',
                         action='store_true' if VALID_CONFIGS['make_NR_file'] else 'store_false',
                         help='Sets make-NR-file to {}'
                              .format(not VALID_CONFIGS['make_NR_file']))
config_args.add_argument('--make-IC-file',
                         action='store_true' if VALID_CONFIGS['make_IC_file'] else 'store_false',
                         help='Sets make-IC-file value to {}'
                              .format(not VALID_CONFIGS['make_IC_file']))
config_args.add_argument('--try-reload',
                         action='store_true' if VALID_CONFIGS['try_reload'] else 'store_false',
                         help='Sets try-reload to {}'
                              .format(not VALID_CONFIGS['try_reload']))
config_args.add_argument('--format-type', default=VALID_CONFIGS['IC_format_type'],
                         help='Type of the format (binary...)')
config_args.add_argument('--MESA-file',
                         default=VALID_CONFIGS['MESA_file'],

                         help='Path to input MESA output')
config_args.add_argument('--masscut', default=VALID_CONFIGS['masscut'],
                         type=float,
                         help='Sets masscut')
config_args.add_argument('--N', default=VALID_CONFIGS['N'],
                         type=int,
                         help='Sets N')
config_args.add_argument('--mp', default=VALID_CONFIGS['mp'],
                         type=float,
                         help='Set the mp value in Msolar units')
config_args.add_argument('--mp_cgs', default=VALID_CONFIGS['mp_cgs'],
                         type=float,
                         help='Set the mp value in cgs units')
config_args.add_argument('--step-size', default=VALID_CONFIGS['stepsize'],
                         type=float,
                         help='Sets the step size to the given value. Units in cm')
# config_args.add_argument('--star-type', default=VALID_CONFIGS['startype'],
#                          help='Set start type')
config_args.add_argument('--loaded_NR_filename', default=None,
                         help='Set the NR file')
config_args.add_argument('--loaded_IC_filename', default=None,
                         help='Set the IC file')

args = parser.parse_args()
if args.defaults:
    for config, value in VALID_CONFIGS.items():
        print("{} = {}".format(config, value))
    exit(0)

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

    def get_path_option(name, default):
        m = re.search('^\s*{}\s*=\s*(?P<value>[\w\./]+)'.format(name),
                      config_file, re.MULTILINE)
        if m:
            value = m.group('value')
            print("Using config file value for {} of {}".format(name, value))
            return value
        else:
            if re.search(name, config_file):
                print("{} was malformed in config file".format(name))
            return default

    def get_str_option(name, default):
        m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                      config_file, re.MULTILINE)
        if m:
            value = m.group('value')
            print("Using config file value for {} of {}".format(name, value))
            return value
        else:
            if re.search(name, config_file):
                print("{} was malformed in config file".format(name))
            return default

    def get_float_option(name, default):
        m = re.search('^\s*{}\s*=\s*(?P<value>[\w\.-]+)'.format(name),
                      config_file, re.MULTILINE)
        if m:
            try:
                value = m.group('value')
                value = float(value)
                print("Using config file value for {} of {}".format(name, value))
                return value
            except ValueError:
                print("{} has a misformed value. {} is not a recognized float"
                      .format(name, value))
        else:
            if re.search(name, config_file):
                print("{} was malformed in config file".format(name))
            return default

    def get_int_option(name, default):
        m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                      config_file, re.MULTILINE)
        if m:
            try:
                value = m.group('value')
                value = int(value)
                print("Using config file value for {} of {}".format(name, value))
                return value
            except ValueError:
                print("{} has a misformed value. {} is not a recognized float"
                      .format(name, value))
        else:
            if re.search(name, config_file):
                print("{} was malformed in config file".format(name))
            return default

    def get_bool_option(name, default):
        m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                      config_file, re.MULTILINE)
        if m:
            # case insensitive true
            value = m.group('value')
            if value.upper() == "TRUE":
                print("Using config file value for {} of {}".format(name, True))
                return True
            else:
                if value.upper() == "FALSE":
                    print("Using config file value for {} of {}".format(name, False))
                    return False
            print("{} needs to have true or false and is {}".format(name, value))
            return default
        else:
            if re.search(name, config_file):
                print("{} was malformed in config file".format(name))
            return default

    check_MESA_profile = get_bool_option('check_MESA_profile', VALID_CONFIGS['check_MESA_profile'])
    make_NR_file = get_bool_option('make_NR_file', VALID_CONFIGS['make_NR_file'])
    make_IC_file = get_bool_option('make_IC_file', VALID_CONFIGS['make_IC_file'])
    try_reload = get_bool_option('try_reload', VALID_CONFIGS['try_reload'])
    IC_format_type = get_str_option('IC_format_type', VALID_CONFIGS['IC_format_type'])


    MESA_file = get_path_option('MESA_file', VALID_CONFIGS['MESA_file'])
    masscut = get_float_option('masscut', VALID_CONFIGS['masscut'])
    N = get_int_option('N', VALID_CONFIGS['N'])
    mp = get_float_option('mp', VALID_CONFIGS['mp']) ##IN UNITES OF Msolar!!!
    mp_cgs = get_float_option('mp_cgs', VALID_CONFIGS['mp_cgs'])
    stepsize = get_float_option('stepsize', VALID_CONFIGS['stepsize']) ##IN UNITES OF cm

    #startype = get_str_option('startype', VALID_CONFIGS['startype'])
    png_tag = get_str_option('png_tag',VALID_CONFIGS['png_tag'])
    


    if make_NR_file:
        new_NR_filename = get_path_option('new_NR_filename',VALID_CONFIGS['new_NR_filename'])
    else:
        new_NR_filename = get_path_option('loaded_NR_filename', VALID_CONFIGS['loaded_NR_filename'])

    if make_IC_file:
        new_IC_filename = get_path_option('new_IC_filename',VALID_CONFIGS['new_IC_filename'])
    else:
        new_IC_filename = get_path_option('loaded_IC_filename', VALID_CONFIGS['loaded_IC_filename'])

else:
    check_MESA_profile = args.check_MESA_profile
    make_NR_file = args.make_NR_file
    make_IC_file = args.make_IC_file
    #new_NR_filename = args.new_NR_filename
    #new_IC_filename = args.new_IC_filename

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


    ## this is probably wrong
    if args.make_NR_file:
        nrfile = args.new_NR_filename
    else: 
        nrfile = args.loaded_NR_filename#"work/NR_files/saveNR_ms_logE.dat"



# 
#
### Meridith screwing things up
if make_NR_file:
    nrfile=new_NR_filename#loaded_NR_filename
else:
    nrfile=loaded_NR_filename

if make_IC_file:
    icfile=new_IC_filename
else:        
    icfile=loaded_IC_filename
#
#
#


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
    print '\n\nInconsistent values of mp specified!'
    sys.exit()

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
if make_IC_file:
    t2=time.time()
    print '\n\nGenerating IC file...'
    in_file=nrfile
    out_file=icfile
    mn.get_IC(in_file,out_file,mp,format_type=IC_format_type)
    print 'IC file generation complete!'
    print("--- %s seconds ---" % (time.time() - t2))


###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################
if try_reload:
    r_temp, rho_temp=mn.reload_IC(icfile,IC_format_type)

    nbin=70.
    r_reloaded,rho_reloaded=mn.binned_r_rho(r_temp, nbin, mp)

    fit_region_R=mn.MESA_r(MESA_file, masscut)
    fit_region_rho=mn.MESA_rho(MESA_file, masscut)

    plt.plot(r_reloaded, rho_reloaded,'r.', markersize=6, label='GADGET data')
    plt.plot(fit_region_R, fit_region_rho, "b.", markersize=4, label='MESA data') #cf.to_log()
    plt.xlabel("R")
    plt.ylabel("test density")
    plt.legend(loc=1)
    if IC_format_type=='hdf5':
    	plt.savefig('lin_'+png_tag+'_hdf5.png')
    else:
    	plt.savefig('lin_'+png_tag+'_bin.png')
    plt.close()

    plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=4, label='MESA data') #cf.to_log()
    plt.plot(r_reloaded, cf.to_log(rho_reloaded),'r.', markersize=6, label='GADGET data')
    #plt.ylim(-2.5,0.3)
    plt.xlabel("R")
    plt.ylabel("log(test density)")
    plt.legend(loc=1)
    if IC_format_type=='hdf5':
    	plt.savefig('log_'+png_tag+'_hdf5.png')
    else:
    	plt.savefig('log_'+png_tag+'_bin.png')
        plt.close()


print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
