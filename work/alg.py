#!/usr/bin/env python
import ConfigParser
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import re
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/src/')
try:
    import MESA2GADGET.mesalib.MESAlibjoyce as MJ
    import MESA2GADGET.mesalib.converge_funcs as cf
    import MESA2GADGET.mesalib.io_lib as rw
    import MESA2GADGET.mesalib.mainlib as mn
    from MESA2GADGET import MESA_PKG_DIR
except ImportError:
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
    'check_MESA': False,
    'make_NR_file': True,
    'make_IC_file': True,
    'try_reload': False,
    'format_type': 'binary',
    'MESA_file': path_from_package('out/sample_MESA_output/profile_OB.data'),
    'masscut': 0.95,
    'N': 16,
    'mp': 1e-7,
    'startype': 'OB_latest',
    'saveNR': 'work/NR_files/saveNR_ms.dat',
    'tag': 'main_sequence'}
    
#check_MESA=False
#make_NR_file=True
#make_IC_file=True
#try_reload=True
#format_type='binary' 
#MESA_file='../out/sample_MESA_output/profile_OB.data'
#masscut=0.95
#N=32
#mp=1e-5 ##IN UNITES OF Msolar!!!

#startype='OB_latest'#'wd_from_mod'
#tag=startype+'_m'+str(masscut)+'_N'+str(N)+'_'+'mp'+str(mp)
#outname=tag

#saveNR="./NR_files/saveNR_"+startype+".dat"

##############################################################
#
# Argument input
#
##############################################################

### edited defaults!!!!

parser = argparse.ArgumentParser(description='Program for converting MESA output to Gadget simulation files')
parser.add_argument('--config-file',
                    help='Path to configuration file which should be in \'INI\' format')
config_args = parser.add_argument_group("Configuration")

# Boolean option flags set opposite of file default
config_args.add_argument('--check-MESA',
                         action='store_true' if VALID_CONFIGS['check_MESA'] else 'store_false',
                         help='Sets check-MESA value to {}'
                              .format(not VALID_CONFIGS['check_MESA']))
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
config_args.add_argument('--format-type', default=VALID_CONFIGS['format_type'],
                         help='Type of the format (binary...)')
config_args.add_argument('--MESA-file',
                         default=VALID_CONFIGS['MESA_file'],
                         help='Path to input MESA output')
config_args.add_argument('--masscut', default=VALID_CONFIGS['masscut'],
                         help='Sets masscut')
config_args.add_argument('--N', default=VALID_CONFIGS['N'],
                         help='Sets N')
config_args.add_argument('--mp', default=VALID_CONFIGS['mp'],
                         help='Set the mp value in Msolar units')
config_args.add_argument('--star-type', default=VALID_CONFIGS['startype'],
                         help='Set start type')
config_args.add_argument('--saveNR', default=None,
                         help='Set the saveNR file')

args = parser.parse_args()
if args.config_file:
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

    def get_str_option(name, default):
        m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                      config_file, re.MULTILINE)
        if m:
            value = m.group('value')
            print("Using config file value for {} of {}".format(name, value))
            return value
        else:
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
            return default

    check_MESA = get_bool_option('check_MESA', VALID_CONFIGS['check_MESA'])
    make_NR_file = get_bool_option('make_NR_file', VALID_CONFIGS['make_NR_file'])
    make_IC_file = get_bool_option('make_IC_file', VALID_CONFIGS['make_IC_file'])
    try_reload = get_bool_option('try_reload', VALID_CONFIGS['try_reload'])
    format_type = get_str_option('format_type', VALID_CONFIGS['format_type'])


    MESA_file = get_str_option('MESA_file', VALID_CONFIGS['MESA_file'])
    masscut = get_float_option('masscut', VALID_CONFIGS['masscut'])
    N = get_int_option('N', VALID_CONFIGS['N'])
    mp = get_float_option('mp', VALID_CONFIGS['mp']) ##IN UNITES OF Msolar!!!

    startype = get_str_option('startype', VALID_CONFIGS['startype'])
    tag = VALID_CONFIGS['tag']
    outname = tag

    saveNR=get_str_option('saveNR', VALID_CONFIGS['saveNR'])
    
else:
    check_MESA = args.check_MESA
    make_NR_file = args.make_NR_file
    make_IC_file = args.make_IC_file
    try_reload = args.try_reload
    format_type = args.format_type
    MESA_file = args.MESA_file
    masscut = args.masscut
    N = args.N
    mp = args.mp
    startype = args.star_type
    #tag = startype + '_m' + str(masscut) + '_N' + str(N) + '_' + 'mp' + str(mp)
    tag = 'main_sequence'
    outname = tag
    # If a saveNR path is given use that
    # otherwise construct one from other variables
    if args.saveNR:
       saveNR = args.saveNR
    else: 
        saveNR = "work/NR_files/saveNR_ms.dat"

# If the given path exists use that otherwise
# prepend the package path
if not os.path.exists(saveNR):
    saveNR = path_from_package(saveNR)
        
        
    

##############################################################
#
# Generate shell placement radii
#
#############################################################
rough_Nshells=1000.
#outer_step=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
stepsize=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
print "estimated stepsize: ", stepsize

stepsize=2.45e6 #87580000 #(cm)


if make_NR_file:
	outf=open(saveNR,"w")
	mn.make_NR_file(MESA_file,masscut,N,mp, stepsize,outf,check_MESA=check_MESA)
	outf.close()

##############################################################
#
# Convert placement radii into HEALPix shells and send to IC
#
#############################################################
fit_region_R=mn.MESA_r(MESA_file, masscut)
fit_region_rho=mn.MESA_rho(MESA_file, masscut)
rmax=fit_region_R.max()

print "rmax: ", rmax

if make_IC_file:
	mn.get_IC(saveNR,outname,mp,format_type=format_type)


###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################
if try_reload:
	r_temp, rho_temp=mn.reload_IC(outname,format_type)


	plt.plot(r_temp, rho_temp,'r.', markersize=6, label='GADGET data')
	plt.plot(fit_region_R, fit_region_rho, "b.", markersize=4, label='MESA data') #cf.to_log()
	plt.xlabel("R")
	plt.ylabel("test density")
	plt.legend(loc=1)
	if format_type=='hdf5':
		plt.savefig('lin_'+outname+'_hdf5.png')
	else:
		plt.savefig('lin_'+outname+'_bin.png')
	plt.close()


	plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=4, label='MESA data') #cf.to_log()
	plt.plot(r_temp, cf.to_log(rho_temp),'r.', markersize=6, label='GADGET data')
	#plt.ylim(-2.5,0.3)
	plt.xlabel("R")
	plt.ylabel("log(test density)")
	plt.legend(loc=1)
	if format_type=='hdf5':
		plt.savefig('log_'+outname+'_hdf5.png')
	else:
		plt.savefig('log_'+outname+'_bin.png')
	plt.close()


print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
