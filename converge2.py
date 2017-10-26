#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import converge_funcs as cf

Rsolar = 6.955e10 #[cgs]
Msolar= 2.0e30 #[kg]
Msolar_cgs= 2.0e33 #[cgs]

#############################################################################
#
# things that apply to both numerical and analytical
#
##############################################################################

fname='profile175.data'#'profile175.data' #175, 140, 32
MESA_file="{}".format(fname)


r_array = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=0.95, strip=False)
M_array = cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=0.95, strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=0.95 ,strip=False)

r_array=cf.unlog(r_array)
fit_region_rho=cf.unlog(fit_region_rho)

rl=r_array.min()
rmax = r_array.max() 

################################################################################
################################################################################
################################################################################
#first_mp=get_first_mp(MESA_file, n_p_initial)

n_p_initial=3000 #72#3072#100,000 
stepsize=0.001#

mp=1e-8#
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


# #------------------------------- ANALYTIC -------------------------------------


# guess_a=1e-7#1000e5
# guess_b=0.68
# guess_c=1


# #fit_region_R = r_array#cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=0.95,strip=False)
# #fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=0.95 ,strip=False)

# fit_region_R=r_array
# #fit_region_rho=cf.unlog(fit_region_rho)

# A,B,C,y=cf.get_curve(fit_region_R,fit_region_rho,guess_a,guess_b,guess_c)

# fit= "g(y) = %s [1/(x + %s)] + %s" % (A, B, C)

# plt.plot(fit_region_R, y(fit_region_R),"r-", markersize=1, label=fit)
# plt.plot(fit_region_R, fit_region_rho, "b.", markersize=1, label='MESA data')
# #plt.ylim(-2,10)
# plt.legend(loc=1, fontsize='x-small')
# plt.xlabel("radius R (unlogged)")
# plt.ylabel("density (unlogged)")
# plt.savefig('fit_region.png')
# plt.close()	


# # rl=fit_region_R.min()#2.9 #from fit range to tail of MESA profile; CAN'T BE THE SAME AS ru[0] or log will freak out
# # rmax=fit_region_R.max()

# test_num_part=3072#16#10.0e7
# stepsize=10e-3# 10e-6
# test_ru=rl + stepsize

# mp = cf.approximate_mp(rl, test_ru, test_num_part, A, B, C)
# print "value of mp: ", mp

# tolerance = 0.1
# shell_radii=[]

# shell_ru=rl
# while shell_ru < rmax:
# 	shell_ru= cf.get_N_continuous(shell_ru, rmax, A, B, C, mp, stepsize, tolerance)
# 	print 'ru found: ', shell_ru
# 	shell_radii.append(shell_ru)

# ru_full=np.arange(fit_region_R.min(),rmax,stepsize)


# for r in range(len(shell_radii)):
# 	rl=shell_radii[r]
# 	ru_list=np.arange(rl,rmax,stepsize)
# 	#plt.plot(ru_list, cf.calc_n1(cf.Mshell_from_integral(rl,ru_list,A,B,C),mp), 'g-', label='n1(ru)')
# 	plt.plot(ru_list, cf.calc_n2(rl,ru_list), 'm-', label='n2($r_u$) for $r_u=$'+str(rl))

# plt.plot(ru_full, cf.calc_n1(cf.Mshell_from_integral(fit_region_R.min(),ru_full,A,B,C),mp), 'g-', label='n1($r_u$)')
# plt.legend(loc=3, fontsize='x-small')
# plt.xlabel("Bounding Radius $r_u$")
# plt.ylabel("$N$")
# plt.ylim(0,30)
# plt.savefig('find_intersection.png')
# plt.close()
 


# print "shell radii through analytical fit: ", shell_radii
# print "shell radii through numerical fit: ", shell_radii_num

