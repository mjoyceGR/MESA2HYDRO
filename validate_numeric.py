#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import converge_funcs as cf
import MESAlibjoyce as MJ

#############################################################################
#
# numerical 
#
##############################################################################

fname='profile175.data'#'profile175.data' #175, 140, 32
MESA_file="{}".format(fname)

masscut=0.95

#mass cut is what does it.
r_array = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut, strip=False)
M_array = cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut, strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=0.95 ,strip=False)

r_array=cf.unlog(r_array)
fit_region_rho=cf.unlog(fit_region_rho)

rl=r_array.min()
rmax = r_array.max() 


#weirdn=8.0
#n_p_initial=12.0*weirdn**2.0 #3072 #72#3072#100,000 


#first_mp=cf.get_first_mp(MESA_file, n_p_initial)
#first_mp=1e-8

# no. Np isn't physical. Mp is.

mp=1e-8#cf.get_mp_given_N(r_array, M_array, n_p_initial)


stepsize=10e-6#0.001#

#to get N to go down by a factor of 10, probably decrease step size by factor of 10

#mp=1e-5#first_mp#1e-05#



print "\nmp: ", mp#, " first np: ", n_p_initial
ru = rl + stepsize


N_r_dict={}
while ru < rmax:
	new_vals= cf.get_N(r_array, M_array, ru, ru, mp, stepsize) 
	try:
		rl=new_vals[0]
	except TypeError:
		print '\n\nRoutine terminated\n\n'
		break	
	ru=new_vals[1]
	rmid=new_vals[2]
	N=new_vals[3]
	n_p=new_vals[4]

	#shell_radii_vals.append()
	N_r_dict[rmid]=int(N)

	#print rl,'\t',ru, "\t", N, "\t", n_p, "\t", mp, "\t", stepsize
#outf.close()
#print 'profile data generated in file', filename

# cf.do_converge(MESA_file, r_array, M_array, n_p_initial,stepsize,\
# outputfile='test_vals_after_overhaul.dat', mp=1e-5)#0.000102319572081 )#)0.1)
# ## making it pick "mp" with my estimator is definitely breaking this


dkeys=N_r_dict.keys()
dkeys.sort()

super_x=[]
super_y=[]
super_z=[]
out_fname='test.hdf5'
for key in dkeys:
	print N_r_dict[key], key

	NSIDE= N_r_dict[key]
	r_mid= key

	x,y,z=MJ.get_coords(NSIDE,r_mid)

	#for i in range(len(x)):
	print "should be: ", r_mid, " is: ", x[3]**2.0 + y[3]**2.0  + z[3]**2.0
	super_x=super_x+x
	super_y=super_y+y
	super_z=super_z+z 
	print "len(super_x): ", len(super_x)



#healpix_file='healpix_to_gadget_shell_1.dat'


# coord_file='combined.dat'

# mp=1e-8

# var=MJ.make_IC_box_hdf5(out_fname,mp, coord_file)
# print var
# print type(var)



