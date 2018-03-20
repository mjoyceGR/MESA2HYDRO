#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf
#import read_write_HDF5 as rw
import io_lib as rw


import time
start_time = time.time()


#########################
#
# copied from validate_numeric
#
##########################


#-------------------
Rsolar = 6.955e10 #[cgs]
Msolar= 2.0e30 #[kg]
Msolar_cgs= 2.0e33 #[cgs]
#----------------


run=True
#run=True
#make_hdf5=False
make_hdf5=True

fname='../data/profile140.data' #140, 32
MESA_file="{}".format(fname)
guess_a=1e-7#1000e5
guess_b=0.68
guess_c=1

masscut=0.95

fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
fit_region_R=cf.unlog(fit_region_R)
fit_region_rho=cf.unlog(fit_region_rho)


A,B,C,y_func=cf.get_curve(fit_region_R,fit_region_rho,guess_a,guess_b,guess_c)

fit= r"g(y) = %s [1/(x + %s)] + %s"%(A, B, C)
#fit=r"analytical fit $g(y)$ = %s [1/($r$ + %s)] + %s"%(4.9,60.9,-11)

#print "type of fit_region_R: ", type(fit_region_R)
analytic_rho=y_func(fit_region_R)


rl=fit_region_R.min()#2.9 #from fit range to tail of MESA profile; CAN'T BE THE SAME AS ru[0] or log will freak out
rmax=fit_region_R.max()

saveNR="../out/saveNR_03.dat"


########################
#
# working: stepsize=10e-6
# masscut=0.91
# mp=10e-6
# force_N=64
#
########################

stepsize=10e-3
mp=10e-6
force_N=16

test_ru=rl + stepsize


hdf5file='mp'+str(mp)+ '_ss' +str(stepsize)+ '.hdf5'#'lowerN.0.hdf5'

#---------------------------------------------------------------
if run:
	outf=open(saveNR,"w")
	N_r_dict={}
	shell_ru=rl
	while shell_ru < rmax:
		rl, shell_ru, rmid, N= cf.get_N_continuous(shell_ru, rmax, A, B, C, mp, stepsize)
		N_r_dict[rmid]=int(N)
		print >> outf, N, rmid, rl, shell_ru

	outf.close()


N,rmid, rl_list, ru_list=np.loadtxt(saveNR, usecols=(0,1,2,3), unpack=True)

if make_hdf5:
	N=force_N + 0.*N
	test=[0]
	super_x=[]
	super_y=[]
	super_z=[]

	super_rho=[] #pull out rough density associated with radii that were saved in 
	#NEED TO GET THIS DIRECTLY FROM MESA

	for i in range(len(N)):
		print N[i], rmid[i]#N_r_dict[key], key

		NSIDE= N[i]#+1#N_r_dict[key]
		r_mid= rmid[i]#key

		######## THIS IS THE WHOLE THING THAT DETERMINES WHETHER IT'S ANALYTIC OR NUMERICAL
		found_rho_idx=np.where( (r_mid>fit_region_R-2) & (r_mid<fit_region_R+2) )[0]
		found_rho=np.average( fit_region_rho[found_rho_idx])


		x,y,z=MJ.get_coords(NSIDE,r_mid, rmax) #now rmax included in coord function
		rho =found_rho+0.0*x

		print "\nshould be: ", r_mid, "\tr_mid/rmax: ",r_mid/rmax,\
		"\tis: ", np.sqrt(x[test]**2.0 + y[test]**2.0  + z[test]**2.0),\
		"\t with rho: ", rho[test]
		#"\n\nlen(super_x): ", len(super_x)
		#print x

		super_x=super_x+list(x)
		super_y=super_y+list(y)
		super_z=super_z+list(z)
		super_rho = super_rho+list(rho)

	print len(super_x), type(super_x)


	super_x=MJ.to_array(super_x)
	super_y=MJ.to_array(super_y)
	super_z=MJ.to_array(super_z)
	super_rho=MJ.to_array(super_rho)

	out_fname=hdf5file
	var=MJ.make_IC_hdf5(out_fname, mp, super_x, super_y, super_z, super_rho, userho=True)
	print var, type(var)

hdf5_file=hdf5file 
PartType=0

density=rw.read_block_single_file(hdf5_file,"Density",PartType)[0][:]
#temp = rw.read_block_single_file(hdf5_file, "Temperature",PartType)[0][:]

x=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,0]
y=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,1]
z=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,2]
r= np.sqrt(x**2.0 + y**2.0 + z**2.0)*rmax



#numeric_healpix = 
#print "numeric_healpix: ", numeric_healpix


plt.plot(fit_region_R, fit_region_rho, "b.", markersize=2, label='MESA data')
plt.plot(fit_region_R, analytic_rho,"y-", linewidth=2, label="Analytic fit")

#plt.plot(r, analytic_healpix,'k--',linewidth=0.5,label='Profile returned from HEALPix')
# DON'T NEED A HEALPIX VALIDATOR IN THIS CASE and also I have no idea how to make one; 
# the devations should be captured in the HDF5 file

plt.plot(r, density, 'ms', markersize=3,label='Profile returned from GADGET2')




plt.legend(loc=1, fontsize='small')
plt.xlabel(r'$R_{\odot}$')
plt.ylabel(r'$\rho$')
plt.savefig('NUMERIC_rho_v_r_Gadget_hdf5.png')
plt.close()

print "execution length: "
print("--- %s seconds ---" % (time.time() - start_time))