#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import MESAlibjoyce as MJ
import converge_funcs as cf
import io_lib as rw
import time
start_time = time.time()

M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar

run=False
make_IC_file=True
tryreload=True

#fname='../data/profile140.data' #140, 32
#tag='profile140'

#MESA_file='../data/AGB/profile31.data'
MESA_file='../data/profile140.data'
#'profile_mainsequence.data'
masscut=0.95#999
tag='p140_m'+str(masscut)


fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
#fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='rho', masscut=masscut ,strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
fit_region_rho=cf.unlog(fit_region_rho)

############# WARNING! NOT CONVERTED!
fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)



fit_region_R=cf.unlog(fit_region_R)*R_to_solar


# #plt.plot(fit_region_R, fit_region_M*M_to_solar,'c.')
# plt.plot(fit_region_R, fit_region_rho,'c.')
# plt.show()
# plt.close()

##############################################################
#
# initial conditions, see N_mp_combinations.dat for suggested step sizes
#
#############################################################

test_r=0.75e13#0.76e12
test_rho=cf.rho_r(test_r, MESA_file,0.95,load_unlogged=False)
test_mass= fit_region_M[cf.find_nearest(fit_region_R,test_r)[1]]*M_to_solar
print "rho at r=",test_r ,test_rho
print "mass at r=",test_r, test_mass

print "force density: ", test_mass/( (4.0/3.0)*np.pi * (test_r)**3.0)

N=8#16#8
mp=1e-7#2e-   # in terms of solar masses OR grams, but need rest of units to match
mp=mp*M_to_solar
print "mass per particle in solar: ", mp
#print 'mp', mp


#print cf.target_Mshell(N,mp)
print "target Mshell", cf.target_Mshell(N,mp),\
 		"\ttarget Mshell/totalM, solar conv",cf.target_Mshell(N,mp)/(fit_region_M.max()*M_to_solar),\
		"\ttarget Mshell/totalM, not solar",cf.target_Mshell(N,mp/M_to_solar)/fit_region_M.max()

target_mass_per_shell=cf.target_Mshell(N,mp)#*M_to_solar
#sys.exit()


rl=fit_region_R.min()
rmax=fit_region_R.max()
print "rl, rmax: ", rl,rmax

ru=rl
RKstep=(rmax-rl)/10000.#target_mass_per_shell/10.#10e-4
outer_step=RKstep

print "RKstep, outer step, percentage step of Rtot ", RKstep,outer_step, RKstep/rmax
#sys.exit()

#---------------------------------------------------------------
saveNR="saveNR_"+tag+".dat" #_long
if run:
	outf=open(saveNR,"w")
	print >> outf, '## fname',MESA_file ,' masscut',masscut,'   N', N, '  mp', mp, '  RKstep', RKstep, '   outer_step', outer_step
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

print time.time()-start_time


hdf5file=tag+ '.hdf5'
binaryfile=tag+ '.bin'
N,rmid=np.loadtxt(saveNR, usecols=(0,1), unpack=True) #rl_list, ru_list

if make_IC_file:
	#test=0
	super_x=[]
	super_y=[]
	super_z=[]

	#temp=open('temp.dat',"w")
	for i in range(len(N)):
		NSIDE= N[i]
		r_mid= rmid[i]
		radius=float(rmid[i])/float(rmax)

		###############################################################
		#
		# Rotation takes place here
		#
		###############################################################
		x,y,z=cf.get_coords(NSIDE)#*(r_mid/rmax)
		#,r_mid, rmax) #now rmax included in coord function
		#print "x,y,z: ", x,y,z

		x=x*radius
		y=y*radius
		z=z*radius

		# print "actual x values: ", x[0:5:1]
		# print "actual y values: ", y[0:5:1]
		# print "actual z values: ", z[0:5:1]


		r_gc=np.sqrt(x**2.0 + y**2.0 + z**2.0)#*rmax#[0:10:1]
		print "\n\nr before get_coords:", r_mid, "\nr after get_coords:", r_gc[-1] , "loop index", i
		#######################################################
		#
		# CHECK IF SCREWUP IS HERE
		#
		######################################################
		#print "super_x", super_x[0:10:1]

		super_x=np.concatenate((super_x,x),axis=0)
		super_y=np.concatenate((super_y,y),axis=0)
		super_z=np.concatenate((super_z,z),axis=0)
		# print "\nsuper_x after concat: ", super_x[0:5:1]
		# print "\super_y after concat: ", super_y[0:5:1]
		# print "\super_z after concat: ", super_z[0:5:1]


		r_super=np.sqrt(super_x**2.0+ super_y**2.0 + super_z**2.0)#*rmax
		print "\nr after append to supers:", r_super[-1]
		print '\n'
		#break

	super_x=cf.to_array(super_x)
	super_y=cf.to_array(super_y)
	super_z=cf.to_array(super_z)
	r_super=np.sqrt(super_x**2.0+ super_y**2.0 + super_z**2.0)
	print "\nr after append to ALL supers:", r_super

	#print "yxz[6] sent to binary", super_x[6], super_y[6], super_z[6]
	out_fname=binaryfile
	var=rw.make_IC_binary(out_fname, mp, super_x, super_y, super_z, central_mass=1,userho=True) #super_rho
	print var, type(var)




if tryreload:
	N,r_set,m_cont=np.loadtxt(saveNR, usecols=(0,1,2), unpack=True)

	print "\n\n\n<----testing reload---->"

	f=open(tag+'.bin','r')#open('ms.bin',"r")
	ptype=0

	header=rw.load_gadget_binary_header(f)
	attribute_dictionary=rw.load_gadget_binary_particledat(f, header, ptype, skip_bh=0)
	positions=attribute_dictionary['Coordinates']
	p_masses=attribute_dictionary['Masses'] ##all of these are the same value, mp
	x_recovered=positions[:,0]
	y_recovered=positions[:,1]
	z_recovered=positions[:,2]

	r_recovered=np.sqrt(x_recovered**2.0 + y_recovered**2.0 + z_recovered**2.0)*rmax
	p_mass=p_masses[5]


	###############################################
	#
	# sort the data???
	#
	###############################################
	#r=r_recovered.sort()

	r1=r_recovered.min() #0.70e13#r_recovered[i+1]
	#rbin=0.0096e13
	r_temp=[]
	rho_temp=[]
	for i in range(len(r_set)-1):
	#range(int(np.floor(r_recovered.max()/rbin  ))):
		#r2=r1+rbin
		r1=r_set[i]
		r2=r_set[i+1]
		print "r1, r2", r1, r2
		region=np.where( (r1<=r_recovered) &(r2>r_recovered))
		print "length of r_recovered", len(r_recovered[region])

		if len(r_recovered[region])==0:
			break
		print "length of r_recovered*mp", len(r_recovered[region])*p_mass
		print "volume of shell r2-r1: ", cf.volume(r2)-cf.volume(r1)
		print "r_rec[ r1 to r2 ]*mp / vol(r2-r1): ", len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))

		r_temp.append(r2)#( (r2+r1)/2.0 )
		rho_temp.append( len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))  )

		#r1=r1+rbin

	print "\n\n"


	print 'length of x_recovered', len(x_recovered)
	print 'length of p_masses', len(p_masses)
	print 'length of recovered radii list', len(r_recovered)
	print 'rmin recovered',r_recovered.min(), '   rmax recovered', r_recovered.max()



	#plt.plot(test_r, test_cm,'r.', markersize=4, label='GADGET data')
	#plt.plot(fit_region_R, fit_region_M*M_to_solar, "b.", markersize=6, label='MESA data')
	plt.plot(r_temp, rho_temp,'r.', markersize=4, label='GADGET data')
	plt.plot(fit_region_R, fit_region_rho, "b.", markersize=6, label='MESA data')
	#plt.ylim(0,1e-7)
	#plt.ylim(0,1e-5)
	plt.xlabel("R")
	plt.ylabel("test density")
	plt.legend(loc=1)
	plt.savefig(tag+'recovered.png')
	#plt.show()
	plt.close()

	# print "execution length: "
	# print("--- %s seconds ---" % (time.time() - start_time))