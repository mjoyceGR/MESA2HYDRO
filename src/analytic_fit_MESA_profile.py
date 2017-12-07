#!/usr/bin/env python
import numpy as np
# import os
# import sys
# import math
import matplotlib.pyplot as plt
# import hdf5lib ## 9/5/17
# import read_write_HDF5 as rw
# import h5py
import MESAlibjoyce as MJ
import converge_funcs as cf

fname='profile140.data' #140, 32
MESA_file="{}".format(fname)
guess_a=1e-7#1000e5
guess_b=0.68
guess_c=1


fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=0.95,strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=0.95 ,strip=False)

fit_region_R=cf.unlog(fit_region_R)
fit_region_rho=cf.unlog(fit_region_rho)

A,B,C,y=cf.get_curve(fit_region_R,fit_region_rho,guess_a,guess_b,guess_c)

fit= "g(y) = %s [1/(x + %s)] + %s" % (A, B, C)

plt.plot(fit_region_R, y(fit_region_R),"r-", markersize=1, label=fit)
plt.plot(fit_region_R, fit_region_rho, "b.", markersize=1, label='MESA data')
#plt.ylim(-2,10)
plt.legend(loc=1, fontsize='x-small')
plt.xlabel("radius R (unlogged)")
plt.ylabel("density (unlogged)")
plt.savefig('fit_region.png')
plt.close()	


rl=fit_region_R.min()#2.9 #from fit range to tail of MESA profile; CAN'T BE THE SAME AS ru[0] or log will freak out
rmax=fit_region_R.max()

test_num_part=3072#16#10.0e7
stepsize=10e-3# 10e-6
test_ru=rl + stepsize

mp = cf.approximate_mp(rl, test_ru, test_num_part, A, B, C)
print "value of mp: ", mp

tolerance = 0.1
shell_radii=[]

shell_ru=rl
while shell_ru < rmax:
	shell_ru= cf.get_N_continuous(shell_ru, rmax, A, B, C, mp, stepsize, tolerance)
	print 'ru found: ', shell_ru
	shell_radii.append(shell_ru)

ru_full=np.arange(fit_region_R.min(),rmax,stepsize)


for r in range(len(shell_radii)):
	rl=shell_radii[r]
	ru_list=np.arange(rl,rmax,stepsize)
	#plt.plot(ru_list, cf.calc_n1(cf.Mshell_from_integral(rl,ru_list,A,B,C),mp), 'g-', label='n1(ru)')
	plt.plot(ru_list, cf.calc_n2(rl,ru_list), 'm-', label='n2($r_u$) for $r_u=$'+str(rl))

plt.plot(ru_full, cf.calc_n1(cf.Mshell_from_integral(fit_region_R.min(),ru_full,A,B,C),mp), 'g-', label='n1($r_u$)')
plt.legend(loc=3, fontsize='x-small')
plt.xlabel("Bounding Radius $r_u$")
plt.ylabel("$N$")
plt.ylim(0,30)
plt.savefig('find_intersection.png')
plt.close()
 

# plt.plot(ru_list, cf.Mshell_from_integral(rl,ru_list,A,B,C), 'g-', label='mass of shell')
# #plt.plot(ru_list, cf.calc_n2(rl,ru_list), 'm-', label='n2(ru)')
# plt.legend(loc=1)#, fontsize='x-small')
# plt.xlabel("Bounding Radius $r_u$")
# plt.ylabel("Mass of shell")
# #plt.ylim(0,20)
# plt.savefig('Mshell.png')
# plt.close()
