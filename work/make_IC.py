#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#m2g_path=os.environ['MESA2GADGET_ROOT']
MESA_PKG_DIR = os.path.abspath(
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '..'))
sys.path.insert(0, os.path.join(MESA_PKG_DIR, 'lib'))
import MESAhandling as MJ
import converge_funcs as cf
import io_lib as rw
import time
start_time = time.time()

###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################


make_IC_file=True
check_MESA=False
filetype='binary'#'hdf5' #or 'hdf5'

MESA_file='../out/sample_MESA_output/profile_whitedwarf_from_mod.data'
saveNR='./NR_files/saveNR_wd_mod.dat'

masscut=0.95
mp=1e-07*M_to_solar
outname='wd_mod'


fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_R=cf.unlog(fit_region_R)*R_to_solar
fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)*M_to_solar
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
fit_region_rho=cf.unlog(fit_region_rho)

if check_MESA:
	plt.plot(fit_region_R, fit_region_rho,'c.')
	plt.show()
	plt.close()

rmax=fit_region_R.max()
rmin=fit_region_R.min()

if make_IC_file:
	rw.get_IC(saveNR,outname,rmax,mp,filetype=filetype)
	print "IC execution length: "
	print("--- %s seconds ---" % (time.time() - start_time))

N,r_set,m_cont=np.loadtxt(saveNR, usecols=(0,1,2), unpack=True)
if filetype=='hdf5':
	hdf5_file=outname+'.hdf5'
	PartType=0
	masses=rw.read_block_single_file(hdf5_file,'Masses',PartType)[0][:]
	x=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,0]
	y=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,1]
	z=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,2]
else:	
	f=open(outname+'.bin','r')
	ptype=0
	header=rw.load_gadget_binary_header(f)
	attribute_dictionary=rw.load_gadget_binary_particledat(f, header, ptype, skip_bh=0)

	positions=attribute_dictionary['Coordinates']
	masses=attribute_dictionary['Masses'] ##all of these are the same value, mp
	x=positions[:,0]
	y=positions[:,1]
	z=positions[:,2]

r_recovered= np.sqrt(x**2.0 + y**2.0 + z**2.0)*rmax
p_mass=masses[5]

r1=r_recovered.min() 
r_temp=[]
rho_temp=[]
for i in range(len(r_set)-1):
	r1=r_set[i]
	r2=r_set[i+1]
	region=np.where( (r1<=r_recovered) &(r2>r_recovered))

	if len(r_recovered[region])==0:
		break
	#print "length of r_recovered*mp", len(r_recovered[region])*p_mass
	#print "volume of shell r2-r1: ", cf.volume(r2)-cf.volume(r1)
	#print "r_rec[ r1 to r2 ]*mp / vol(r2-r1): ", len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))
	r_temp.append(r2)
	rho_temp.append( len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))  )



plt.plot(r_temp, rho_temp,'r.', markersize=4, label='GADGET data')
plt.plot(fit_region_R, fit_region_rho, "b.", markersize=6, label='MESA data') #cf.to_log()
plt.xlabel("R")
plt.ylabel("test density")
plt.legend(loc=1)
if filetype=='hdf5':
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
if filetype=='hdf5':
	plt.savefig('log_'+outname+'_hdf5.png')
else:
	plt.savefig('log_'+outname+'_bin.png')
plt.close()

print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
