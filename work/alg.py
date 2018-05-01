#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/src/')
import MESAlibjoyce as MJ
import converge_funcs as cf
import io_lib as rw
import mainlib as mn

import time
start_time = time.time()

# ###################################################
# M_to_solar=1.988*10.0**33.0 ## g/Msolar
# R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
# ###################################################


#######################################################
#
# user controls
#
########################################################
check_MESA=False
make_NR_file=True
make_IC_file=True
try_reload=True
format_type='binary' 


MESA_file='../out/sample_MESA_output/profile_OB.data'
masscut=0.95
N=32
mp=1e-5 ##IN UNITES OF Msolar!!!

startype='OB_latest'#'wd_from_mod'
tag=startype+'_m'+str(masscut)+'_N'+str(N)+'_'+'mp'+str(mp)
outname=tag

saveNR="./NR_files/saveNR_"+startype+".dat"


##############################################################
#
# Generate shell placement radii
#
#############################################################
rough_Nshells=1000.
#outer_step=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
stepsize=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
print "estimated stepsize: ", stepsize

stepsize=87580000 #(cm)


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
	mn.get_IC(saveNR,outname,rmax,mp,format_type=format_type)


###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################
if try_reload:
	r_temp, rho_temp=mn.reload_IC(rmax, saveNR,outname,format_type)


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