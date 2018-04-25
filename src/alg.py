#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import MESAlibjoyce as MJ
import converge_funcs as cf
import io_lib as rw
import time
start_time = time.time()

###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################


#######################################################
#
# switchboard
#
########################################################
check_MESA=True
run=True
make_IC_file=True
try_reload=True
filetype='binary' 


MESA_file='../out/raw_MESA_output/profile_agb_timmes.data'#profile_mainsequence.data'
masscut=0.95
startype='agb'
N=8
mp=1e-7
tag=startype+'_m'+str(masscut)+'_N'+str(N)+'_'+'mp'+str(mp)
outname=tag

#######################################################
#
# load density, mass, radius from MESA
# and convert to cgs units
#
########################################################
fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_R=cf.unlog(fit_region_R)*R_to_solar

fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)*M_to_solar

fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
fit_region_rho=cf.unlog(fit_region_rho)


if check_MESA:
	plt.plot(fit_region_R, fit_region_rho,'c.')
	plt.show()
	plt.close()



##############################################################
#
# initial conditions
# see N_mp_combinations.dat for suggested step sizes
#
#############################################################
mp=mp*M_to_solar
rl=fit_region_R.min()
rmax=fit_region_R.max()
RKstep=(rmax-rl)/10000. ## divisor here will be roughly the number of shells in the model
outer_step=RKstep

print "target Mshell", cf.target_Mshell(N,mp),\
 		"\ttarget Mshell/totalM, solar conv",cf.target_Mshell(N,mp)/(fit_region_M.max()),\
		"\ttarget Mshell/totalM, not solar",cf.target_Mshell(N,mp/M_to_solar)/(fit_region_M.max()/M_to_solar)


##############################################################
#
# Generate shell placement radii
#
#############################################################
ru=rl
saveNR="saveNR_"+tag+".dat"
if run:
	outf=open(saveNR,"w")
	print >> outf, '## fname',MESA_file ,' masscut',masscut,'   N', N, '  mp', mp/M_to_solar,\
	   'mp_solar', mp,'  RKstep',('%1.3e'% RKstep), '   outer_step', ('%1.3e'%outer_step)
	ru=rl
	while ru <= rmax:
		try:
			ru, Mshell=cf.get_placement_radii(rl, ru, RKstep, N, mp, MESA_file,masscut, suppress_output=False, load_unlogged=False)
			print >> outf, N, ru, Mshell
			rl=ru
			ru=ru + outer_step
		except TypeError:
			print 'reached', ('%1.5f'% (ru/rmax)), r'% of outer radius' 
			break

	print 'runtime: ', time.time()-start_time
	print >> outf, '#\n#\n# runtime: ', time.time()-start_time
	outf.close()




##############################################################
#
# Generate shell placement radii
#
#############################################################
if make_IC_file:
	rw.get_IC(saveNR,outname,rmax,mp,filetype=filetype)

# tag=saveNR.split('.dat')[0]

# hdf5file=tag+ '.hdf5'
# binaryfile=tag+ '.bin'
# N,rmid=np.loadtxt(saveNR, usecols=(0,1), unpack=True) 

# print "number of shells in NR file: ", len(N)

# if make_IC_file:
# 	super_x=[]
# 	super_y=[]
# 	super_z=[]

# 	for i in range(len(N)):
# 		NSIDE= N[i]
# 		r_mid= rmid[i]
# 		radius=float(rmid[i])/float(rmax)

# 		###############################################################
# 		#
# 		# Arbitrarily rotate shells
# 		#
# 		###############################################################
# 		x,y,z=cf.get_coords(NSIDE)

# 		x=x*radius
# 		y=y*radius
# 		z=z*radius

# 		super_x=np.concatenate((super_x,x),axis=0)
# 		super_y=np.concatenate((super_y,y),axis=0)
# 		super_z=np.concatenate((super_z,z),axis=0)

# 		r_super=np.sqrt(super_x**2.0+ super_y**2.0 + super_z**2.0)
		

# 	super_x=cf.to_array(super_x)
# 	super_y=cf.to_array(super_y)
# 	super_z=cf.to_array(super_z)
# 	r_super=np.sqrt(super_x**2.0+ super_y**2.0 + super_z**2.0)

# 	if filetype=='hdf5':
# 		outf_name=hdf5file
# 	else:
# 		out_fname=binaryfile
# 		var=rw.make_IC_binary(out_fname, mp, super_x, super_y, super_z, central_mass=1) 



# 	if filetype=='hdf5':
# 		var=rw.make_IC_hdf5(hdf5file, mp, super_x, super_y, super_z, userho=False) #super_rho
# 	else:
# 		var=rw.make_IC_binary(binaryfile, mp, super_x, super_y, super_z, central_mass=1) 
# 	print var, type(var)



###############################################
#
# test raw profile recovery
# ---NOT a substitution for checking with an SPH viewer!
#
###############################################
if try_reload:
	print "\n\n\n<----testing reload---->"
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
	plt.ylim(-2.5,0.3)
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