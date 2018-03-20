#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf
import read_write_HDF5 as rw

import time
start_time = time.time()


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

#fit= "g(y) = %s [1/(x + %s)] + %s"%(A, B, C)
fit=r"analytical fit $g(y)$ = %s [1/($r$ + %s)] + %s"%(4.9,60.9,-11)

print "type of fit_region_R: ", type(fit_region_R)
analytic_rho=y_func(fit_region_R)

# fig,ax=MJ.plotter(5,1,0.5,0.1,xf='%1.0f')
# plt.plot(fit_region_R, analytic_rho*10.e6,"r-", markersize=1, label=fit)
# plt.plot(fit_region_R, fit_region_rho*10.e6, "bo", markersize=2, alpha=0.5, label='MESA data')
# plt.legend(loc=1, fontsize='small')
# plt.xlabel(r"Radius $R_{\odot}$")# (unlogged)")
# plt.ylabel(r"$\rho$ [10$^{-6}$ g/cm$^2$]")
# plt.savefig('fit_region.png')
# plt.close()	



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

# stepsize=10e-6
# mp=10e-9
# force_N=16


stepsize=10e-6# 10e-6 #then to -7 and -8
test_ru=rl + stepsize
mp=10e-7#8
force_N=512

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
	# dkeys=N_r_dict.keys()
	# dkeys.sort()
	# for key in dkeys:
	# 	print >> outf, N_r_dict[key], key
	outf.close()
#N,r=np.loadtxt('test_N_r.dat', usecols=(0,1), unpack=True)


N,rmid, rl_list, ru_list=np.loadtxt(saveNR, usecols=(0,1,2,3), unpack=True)

if make_hdf5:
	#N,r=np.loadtxt(saveNR, usecols=(0,1), unpack=True)
	N=force_N + 0.*N
	test=[0]#np.arange(0,200,5)
	super_x=[]
	super_y=[]
	super_z=[]

	super_rho=[] #pull out rough density associated with radii that were saved in 

	for i in range(len(N)):
		print N[i], rmid[i]#N_r_dict[key], key

		NSIDE= N[i]#+1#N_r_dict[key]
		r_mid= rmid[i]#key

		#region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		found_rho_idx=np.where( (r_mid>fit_region_R-2) & (r_mid<fit_region_R+2) )[0]
		#found_rho=np.average( fit_region_rho[found_rho_idx])
		found_rho=np.average( analytic_rho[found_rho_idx])


		x,y,z=MJ.get_coords(NSIDE,r_mid, rmax) #now rmax included in coord function
		rho =found_rho+0.0*x

		print "\nshould be: ", r_mid, "\tr_mid/rmax: ",r_mid/rmax,\
		"\tis: ", np.sqrt(x[test]**2.0 + y[test]**2.0  + z[test]**2.0),\
		"\t with rho: ", rho[test]
		#"\n\nlen(super_x): ", len(super_x)
		#print x

		super_x=super_x+list(x)#np.array(x).astype(float)
		super_y=super_y+list(y)#np.array(y).astype(float)
		super_z=super_z+list(z)#np.array(z).astype(float)
		#print "len(super_x): ", len(super_x)
		super_rho = super_rho+list(rho)

	print len(super_x), type(super_x)
	##super_x, 

	super_x=np.array(super_x).astype(float)
	super_y=np.array(super_y).astype(float)
	super_z=np.array(super_z).astype(float)
	super_rho=np.array(super_rho).astype(float)

	out_fname=hdf5file#'multi_shell_test_small_atm.hdf5'
	#'multi_shell_test_NSIDE_512.hdf5'
	var=MJ.make_IC_hdf5(out_fname, mp, super_x, super_y, super_z, super_rho, userho=True)
	print var, type(var)


#------------------------- plots -----------------------------------------------
# shell_radii=ru_list#[1:50:10]#[np.arange(0,50,5).astype(int)] #from loaded saveNR.dat

# #stepsize=10

# ru_full=np.arange(fit_region_R.min(), rmax, stepsize)

# sample=shell_radii[0:35:5]
# print "sample: ", sample 

# for r in range(len(sample)):
# 	rl=sample#shell_radii[r]
# 	ru_list=np.arange(rl,rmax,stepsize)
# 	#plt.plot(ru_list, cf.calc_n1(cf.Mshell_from_integral(rl,ru_list,A,B,C),mp), 'g-', label='n1(ru)')
# 	plt.plot(ru_list, cf.calc_n2(rl,ru_list), 'm-', label='n2($r_u$) for $r_u=$'+str(rl))

# plt.plot(ru_full, cf.calc_n1(cf.Mshell_from_integral(fit_region_R.min(),ru_full,A,B,C),mp), 'g-', label='n1($r_u$)')
# plt.legend(loc=4, fontsize='x-small')
# plt.xlabel("Bounding Radius $r_u$")
# plt.ylabel("$N$")
# plt.ylim(20,70)
# plt.savefig('find_intersection.png')
# plt.close()
 

# plt.plot(ru_list, cf.Mshell_from_integral(rl,ru_list,A,B,C), 'g-', label='mass of shell')
# #plt.plot(ru_list, cf.calc_n2(rl,ru_list), 'm-', label='n2(ru)')
# plt.legend(loc=1)#, fontsize='x-small')
# plt.xlabel("Bounding Radius $r_u$")
# plt.ylabel("Mass of shell")
# #plt.ylim(0,20)
# plt.savefig('Mshell.png')
# plt.close()



hdf5_file=hdf5file #'multi_shell_test_small_atm.hdf5'
PartType=0

density=rw.read_block_single_file(hdf5_file,"Density",PartType)[0][:]
#temp = rw.read_block_single_file(hdf5_file, "Temperature",PartType)[0][:]

x=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,0]
y=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,1]
z=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,2]
r= np.sqrt(x**2.0 + y**2.0 + z**2.0)*rmax



analytic_healpix = y_func((MJ.to_array(r)))
print "analytic_healpix: ", analytic_healpix

plt.plot(fit_region_R, fit_region_rho, "b.", markersize=1, label='MESA data')
plt.plot(fit_region_R, analytic_rho,"y-", linewidth=2, label="Analytic fit")

#plt.plot(r, density,'r-',linewidth=0.5,label='Gadget profile')
plt.plot(r, analytic_healpix,'k--',linewidth=0.5,label='Profile returned from HEALPix')

plt.plot(r, density, 'ms', markersize=3,label='Profile returned from GADGET2')
#plt.ylim(-2,10)
plt.legend(loc=1, fontsize='small')
#plt.plot()
#plt.xlim(0,1)
plt.xlabel(r'$R_{\odot}$')
plt.ylabel(r'$\rho$')
plt.savefig('rho_v_r_Gadget_hdf5.png')
plt.close()

print "execution length: "
print("--- %s seconds ---" % (time.time() - start_time))