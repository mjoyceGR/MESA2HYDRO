#!/usr/bin/env python
import numpy as np
import os
import sys
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import hdf5lib ## 9/5/17

import read_write_HDF5 as rw
import h5py
import MESAlibjoyce as MJ

def fit_MESA_analytic(MESA_file, **kwargs):
	strip=bool(kwargs.get('strip',False))
	unlog=bool(kwargs.get('unlog',False))

	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]	
	# try:
	# 	quantity=MJ.get_quantity(MESA_file,keyword)
	# except:
	# 	print "Quantity keyword not found. Allowed keywords are:"
	# 	print MJ.show_allowed_MESA_keywords(MESA_file)
	##################################################################3	

	mass = MJ.get_quantity(MESA_file,'mass').astype(np.float)
	logR = MJ.get_quantity(MESA_file,'logR').astype(np.float)
	logrho = MJ.get_quantity(MESA_file,'logRho').astype(np.float)

	Mtot=mass[np.argmax(mass)]
	print "Mtot: ", Mtot
	bound = 0.95*Mtot

	region=np.where( (mass>=bound) )[0]
	print "total length: ",len(mass), "  selected length: ",len(region)
	fit_region_R=logR[region]
	fit_region_rho=logrho[region]

	# print "last mass entry (total M??): ", Mtot, "\tlargest mass value: ", max(mass),\
	# "\tmass value at most inner R: ", mass[np.argmin(logR)]
	if unlog:
		print '\n\nusing unlogged variables!\n\n'
		fr_rho = [10.0**p for p in fit_region_rho]
		fr_R = [10.0**p for p in fit_region_R]
		return np.array(fr_R).astype(float), np.array(fr_rho).astype(float) 
	else:
		return np.array(fit_region_R).astype(float), np.array(fit_region_rho).astype(float)




def get_curve(r_array,rho_array,guess_a,guess_b,guess_c):
	def model_func(x,a,b,c):
			## this is a 1/r density profile
		return a*(1.0/(x+b)) + c

	p=[guess_a,guess_b,guess_c]
	params, cov = curve_fit(model_func, r_array, rho_array, p0=p)
	a, b, c = params

	def y(x):
		a,b,c=params
		return a*(1.0/(x+b)) + c

	return a,b,c,y 


#------------------------------------------------------------------------------------
def density_integral(rl,ru, A,B,C):
	rdiff=(ru - rl)
	f1=0.5*A*rdiff**2.0  + A*B*rdiff  + A*B**2.0 * np.log10( abs(rdiff - B) ) + (1.0/3.0)*C*rdiff**3.0
	return f1

### is it necessary to do this density integral or should I just fit the mass profile and pass that?

# def RHS(rl, ru, mp):
# 	f2=(mp/12.0) *( (rl + ru)/(ru - rl) )**2.0
# 	return f2

def calc_n1(rl, ru, mp, A,B,C):
	#this must take ru in unlogged form because that's the space in which the fit function was defined
	Mshell=4.0*np.pi*density_integral(rl, ru, A, B, C)
	return np.sqrt(Mshell/(12.0*mp))


def calc_n2(rl, ru):
	return np.sqrt(np.pi/12.0)*(ru + rl)/(ru - rl)



def approximate_mp(rl, ru, n_p, A, B, C):
	Mshell=4.0*np.pi*density_integral(ru,rl,A,B,C)
	mp=Mshell/n_p
	return mp


fname='profile140.data' #140, 32
MESA_file="{}".format(fname)
guess_a=1e-7#1000e5
guess_b=0.68
guess_c=1

fit_region_R, fit_region_rho=fit_MESA_analytic(MESA_file, strip=False, unlog=True)
A,B,C,y=get_curve(fit_region_R,fit_region_rho,guess_a,guess_b,guess_c)

fit= "g(y) = %s [1/(x + %s)] + %s" % (A, B, C)

plt.plot(fit_region_R, y(fit_region_R),"r-", markersize=1, label=fit)
plt.plot(fit_region_R, fit_region_rho, "b.", markersize=1)
#plt.ylim(-2,10)
plt.legend(loc=1, fontsize='x-small')
plt.xlabel("radius R (unlogged)")
plt.ylabel("density (unlogged)")
plt.savefig('fit_region.png')
plt.close()	



#### once we have the analytic fit, move back into log space?
log_analytic_rho=np.log10(y(fit_region_R))
log_analytic_R=np.log10(fit_region_R)

plt.plot(log_analytic_R, y(fit_region_R),"m-", markersize=1)#, label=fit)
#plt.plot(fit_region_R, fit_region_rho, "b.", markersize=1)
#plt.ylim(-2,10)
plt.legend(loc=1, fontsize='x-small')
plt.xlabel("radius R (logged)")
plt.ylabel("density FIT (unlogged)")
plt.savefig('fit_region_semilog.png')
plt.close()	


#rl=log_analytic_R.min() #
#rmax=log_analytic_R.max() 
rl=fit_region_R.min()#2.9 #from fit range to tail of MESA profile; CAN'T BE THE SAME AS ru[0] or log will freak out
rmax=fit_region_R.max()

test_num_part=3072#16#10.0e7

stepsize=10e-3# 10e-6
test_ru=rl + stepsize

mp = approximate_mp(rl, test_ru, test_num_part, A, B, C)
print "value of mp: ", mp

ru_list=np.arange(rl,rmax,stepsize)
#print ru

shell_ru=0
for i in range(len(ru_list)):
	#while (n2 - n1) >0:
	ru = ru_list[i]
	n1=calc_n1(rl,ru, mp, A, B, C)
	n2=calc_n2(rl, ru)
	#while 
	#if (n1 -n2) >= 0:
	if abs(n1 -n2) <= 0.1:#0.0001:
		print "tolerance met\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru
		#break
		shell_ru=ru


print 'ru found: ', shell_ru

plt.plot(ru_list, calc_n1(rl,ru_list,mp,A,B,C), 'g-', label='n1(ru)')
plt.plot(ru_list, calc_n2(rl,ru_list), 'm-', label='n2(ru)')
plt.legend(loc=1)#, fontsize='x-small')
plt.xlabel("Bounding Radius $r_u$")
plt.ylabel("$N$")
plt.ylim(0,20)
plt.savefig('find_intersection.png')
plt.close()
 