#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import converge_funcs as cf

Rsolar = 6.955e10 #[cgs]
Msolar= 2.0e30 #[kg]
Msolar_cgs= 2.0e33 #[cgs]

#############################################################################
#
# numerical 
#
##############################################################################

fname='profile175.data'#'profile175.data' #175, 140, 32
MESA_file="{}".format(fname)


r_array = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=0.65, strip=False)
M_array = cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=0.65, strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=0.95 ,strip=False)

r_array=cf.unlog(r_array)
fit_region_rho=cf.unlog(fit_region_rho)

rl=r_array.min()
rmax = r_array.max() 

#first_mp=get_first_mp(MESA_file, n_p_initial)

n_p_initial=3000 #72#3072#100,000 
stepsize=0.001#

mp=1e-05#

print "\nmp: ", mp, " first np: ", n_p_initial
ru = rl + stepsize

# outf=open('vals.dat',"a")
# print >> outf, 'MESA_file: ', MESA_file, ' on ', dt.datetime.now()
# print >> outf, 'rl				ru				N		n_p				mp			stepsize'
# print >> outf, rl,'\t',ru, "\tN", "\t", n_p_initial, "\t", mp, "\t", stepsize
shell_radii_num=[]

while ru < rmax:
	new_vals= cf.get_N(r_array, M_array, ru, ru, mp, stepsize, n_p_initial) 
	try:
		rl=new_vals[0]
	except TypeError:
		print '\n\nRoutine terminated\n\n'
		break	
	ru=new_vals[1]
	rmid=new_vals[2]
	N=new_vals[3]
	n_p=new_vals[4]

	shell_radii_num.append(ru)
	#print >> outf, rl,'\t',ru, "\t", N, "\t", n_p, "\t", mp, "\t", stepsize
#outf.close()
#print 'profile data generated in file', filename



# cf.do_converge(MESA_file, r_array, M_array, n_p_initial,stepsize,\
# outputfile='test_vals_after_overhaul.dat', mp=1e-5)#0.000102319572081 )#)0.1)
# ## making it pick "mp" with my estimator is definitely breaking this


