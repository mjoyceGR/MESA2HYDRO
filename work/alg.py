#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/src/')
import MESA2GADGET.mesalib.MESAlibjoyce as MJ
import MESA2GADGET.mesalib.converge_funcs as cf
import MESA2GADGET.mesalib.io_lib as rw
import MESA2GADGET.mesalib.mainlib as mn
from MESA2GADGET import MESA_PKG_DIR

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
def defaults():
    check_MESA=False
    make_NR_file=False
    make_IC_file=True
    try_reload=True
    format_type='binary' 


    MESA_file=os.path.join(MESA_PKG_DIR, 'out/sample_MESA_output/profile140.data')
    masscut=0.95
    N=32
    mp=1e-7 ##IN UNITES OF Msolar!!!

    startype='p140_test'
    tag=startype+'_m'+str(masscut)+'_N'+str(N)+'_'+'mp'+str(mp)
    outname=tag

    saveNR=os.path.join(MESA_PKG_DIR, "work/NR_files/saveNR_ms.dat")


##############################################################
#
# Argument input
#
##############################################################
parser = argparse.ArgumentParser(description='Program for converting MESA output to Gadget simulation files')
parser.add_argument('--config-file', type=file,
                         help='Path to configuration file')
config_args = parser.add_argument_group("Configuration")
config_args.add_argument('--check-MESA', action='store_true',
                         help='Sets check-MESA value')
config_args.add_argument('--make-NR-file', action='store_true',
                         help='Sets make-NR-file')
config_args.add_argument('--no-IC-file', action='store_true',
                         help='Sets make-IC-file value')
config_args.add_argument('--no-reload', action='store_true',
                         help='Sets try-reload')
config_args.add_argument('--format-type', default='binary',
                         help='Type of the format (binary...)')
config_args.add_argument('--MESA-file',
                         default=path_from_package('out/sample_MESA_output/profile140.data'),
                         help='Path to input MESA output')
config_args.add_argument('--masscut', default=0.95,
                         help='Sets masscut')
config_args.add_argument('--N', default=32,
                         help='Sets N')
config_args.add_argument('--mp', default=1e-7,
                         help='Set the mp value in Msolar units')
config_args.add_argument('--start-type', default='p140_test',
                         help='Set start type')
config_args.add_argument('--saveNR', default=None,
                         help='Set the saveNR file')

args = parser.parse_args()

if args.config_file:
    # TODO: Read configurations from a file
    pass
else:
    check_MESA = args.check_MESA
    make_NR_file = args.make_NR_file
    make_IC_file = not args.no_IC_file
    try_reload = not args.no_reload
    format_type = args.format_type
    MESA_file = args.MESA_file
    masscut = args.masscut
    N = args.N
    mp = args.mp
    startype = args.start_type
    tag = startype + '_m' + str(masscut) + '_N' + str(N) + '_' + 'mp' + str(mp)
    outname = tag
    # If a saveNR path is given use that
    # otherwise construct one from other variables
    if args.saveNR:
       saveNR = args.saveNR
    else: 
        saveNR = "saveNR_"+startype+".dat"
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
stepsize=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
print "estimated stepsize: ", '%1.5e'%stepsize

if make_NR_file:
	mn.make_NR_file(MESA_file,masscut,N,mp,stepsize,saveNR,check_MESA=check_MESA)


##############################################################
#
# Convert placement radii into HEALPix shells and send to IC
#
#############################################################
fit_region_R=mn.MESA_r(MESA_file, masscut)
fit_region_rho=mn.MESA_rho(MESA_file, masscut)
rmax=fit_region_R.max()

if make_IC_file:
	mn.get_IC(saveNR,outname,rmax,mp,format_type=format_type)


###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################
if try_reload:
	r_temp, rho_temp=mn.reload_IC(rmax, saveNR,outname,format_type)


	plt.plot(r_temp, rho_temp,'r.', markersize=4, label='GADGET data')
	plt.plot(fit_region_R, fit_region_rho, "b.", markersize=6, label='MESA data') #cf.to_log()
	plt.xlabel("R")
	plt.ylabel("test density")
	plt.legend(loc=1)
	if format_type=='hdf5':
		plt.savefig('lin_'+outname+'_hdf5.png')
	else:
		plt.savefig('lin_'+outname+'_bin.png')
	plt.close()


	plt.plot(r_temp, cf.to_log(rho_temp),'r.', markersize=4, label='GADGET data')
	plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=6, label='MESA data') #cf.to_log()
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
