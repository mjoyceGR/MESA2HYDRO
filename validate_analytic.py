#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf


run=False



fname='profile140.data' #140, 32
MESA_file="{}".format(fname)
guess_a=1e-7#1000e5
guess_b=0.68
guess_c=1

masscut=0.95

fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)

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


stepsize=10e-6#10e-3# 10e-6
test_ru=rl + stepsize
mp=10e-8
if run:
	N_r_dict={}

	shell_ru=rl
	while shell_ru < rmax:
		rl, shell_ru, rmid, N= cf.get_N_continuous(shell_ru, rmax, A, B, C, mp, stepsize)
		N_r_dict[rmid]=int(N)

	dkeys=N_r_dict.keys()
	dkeys.sort()
	for key in dkeys:
		print N_r_dict[key], key
		NSIDE= N_r_dict[key]

#N = 7+0*N
N,r=np.loadtxt('test_N_r.dat', usecols=(0,1), unpack=True)

test=[0]#np.arange(0,200,5)

super_x=[]
super_y=[]
super_z=[]
out_fname='test.hdf5'
for i in range(len(N)):
	print N[i], r[i]#N_r_dict[key], key

	NSIDE= N[i]+1#N_r_dict[key]
	r_mid= r[i]#key

	x,y,z=MJ.get_coords(NSIDE,r_mid, rmax) #now rmax included in coord function

	print "\nshould be: ", r_mid, "\tr_mid/rmax: ",r_mid/rmax,\
	"\nis: ", np.sqrt(x[test]**2.0 + y[test]**2.0  + z[test]**2.0)#,\
	#"\n\nlen(super_x): ", len(super_x)
	#print x

	super_x=super_x+list(x)#np.array(x).astype(float)
	super_y=super_y+list(y)#np.array(y).astype(float)
	super_z=super_z+list(z)#np.array(z).astype(float)
	#print "len(super_x): ", len(super_x)

print len(super_x), type(super_x)
##super_x, 

super_x=np.array(super_x).astype(float)
super_y=np.array(super_y).astype(float)
super_z=np.array(super_z).astype(float)

out_fname='multi_shell_test.hdf5'
var=MJ.make_IC_hdf5(out_fname, mp, super_x, super_y, super_z)
print var, type(var)



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
 

# plt.plot(ru_list, cf.Mshell_from_integral(rl,ru_list,A,B,C), 'g-', label='mass of shell')
# #plt.plot(ru_list, cf.calc_n2(rl,ru_list), 'm-', label='n2(ru)')
# plt.legend(loc=1)#, fontsize='x-small')
# plt.xlabel("Bounding Radius $r_u$")
# plt.ylabel("Mass of shell")
# #plt.ylim(0,20)
# plt.savefig('Mshell.png')
# plt.close()
