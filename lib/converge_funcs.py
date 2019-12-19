#!/usr/bin/env python3

#************************************************************************
#
#   Copyright (C) 2019  M. Joyce, L. Lairmore, D. J. Price
#
#   See MESA2HYDRO/LICENSE
#
#************************************************************************


'''

Contains: numerical methods

'''



import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
from MESA2HYDRO.lib import MESAhandling as MJ
from MESA2HYDRO.lib import mainlib as mn
import datetime as dt 
import random as rand
import healpy as hp
from MESA2HYDRO.lib import constants as const


M_to_solar=const.Msun
R_to_solar=const.Rsun


#######################################################################################
#
# data scanning functions
#
#######################################################################################
def rho_r(r, MESA_file, masscut, *args, **kwargs):
	############################################################
	#
	# Fixed in cgs units!!!
	#
	############################################################
	fit_region_R=mn.MESA_r(MESA_file,masscut)


	fit_region_rho=mn.MESA_rho(MESA_file,masscut)
	r0,idx=find_nearest(fit_region_R,r)
	# WARNING! THIS RELIES ON LOADED DATA BEING SORTED! DO NOT TAMPER!

	if r0 <= r:
		rho0=fit_region_rho[idx]
		r1=fit_region_R[idx-1]
		rho1=fit_region_rho[idx-1]
	else:
		r0=fit_region_R[idx+1]
		rho0=fit_region_rho[idx+1]
		r1=fit_region_R[idx]
		rho1=fit_region_rho[idx]

	rrho_r= (  (r1-r)*rho0 + (r-r0)*rho1 ) /(r1-r0)
	if (r0 <= r <= r1):
		return float(rrho_r)
	else:
		print("converge_funcs WARNING: r not found between r0 and r1....end of MESA density values\n")
		return  


		
def m_r(r, MESA_file, *args, **kwargs):

	all_MESA_m  = MJ.get_quantity(MESA_file,'mass').astype(np.float)
	all_MESA_r  = R_to_solar_f(unlog(MJ.get_quantity(MESA_file,'logR').astype(np.float)))
	fit_region_R = all_MESA_r
	fit_region_m = all_MESA_m

	r0,idx=find_nearest(all_MESA_r,r)
	# WARNING! THIS RELIES ON LOADED DATA BEING SORTED! DO NOT TAMPER!

	if r0 <= r:
		m0=fit_region_m[idx]
		r1=fit_region_R[idx-1]
		m1=fit_region_m[idx-1]
	else:
		r0=fit_region_R[idx+1]
		m0=fit_region_m[idx+1]
		r1=fit_region_R[idx]
		m1=fit_region_m[idx]

	rm_r= (  (r1-r)*m0 + (r-r0)*m1 ) /(r1-r0)
	if (r0 <= r <= r1):
		return float(rm_r)
	else:
		print("converge_funcs WARNING: could not compute m(r) for point")
		return 




#######################################################################################
#
# Numerical integration/integral interpolation
#
#######################################################################################
def target_Mshell(N,mp):
	Mshell=N**2.0*12.0*mp
	return Mshell

def density_integral_numeric(r, rho_r): 
	dmdr=4.0*np.pi*(r**2.0)*float(rho_r)
	return float(dmdr)


def Mshell_from_RK(rl, rmax, step, MESA_file,masscut, *args,**kwargs):
	#use_unlog=bool(kwargs.get('load_unlogged',False))
	### rmax is NOT MODIFIED BY THIS ROUTINE, it is simply a STOPPING CRITERION
	r = rl      
	m=0
	while r<rmax:
	    r, m= RK1(r,m, density_integral_numeric, step, MESA_file,masscut)#, load_unlogged=use_unlog)
	Mshell=m
	return Mshell



def Mshell_from_Romberg(rl, rmax, MESA_file, masscut,  *args,**kwargs): 
	from scipy import integrate
	from scipy.integrate import romberg, quad 
	
	steps=float(kwargs.get("steps", 1)) ## global parameter representing order?
	def fn(r):
		return density_integral_numeric(r, rho_r(r,MESA_file,masscut))
	Mshell = integrate.romberg(fn, rl, rmax, show=False)
	return  Mshell



def Mshell_from_quad(rl, rmax, MESA_file, masscut,  *args,**kwargs): #step, MESA_file,masscut,
	#steps=float(kwargs.get("steps", 1)) ## global parameter representing order?
	def fn(r):
		return density_integral_numeric(r, rho_r(r,MESA_file,masscut))

	Mshell = integrate.quad(fn, rl, rmax)[0]
	return  Mshell


def Newton_Raphson():
	#############################################
	#
	# TBD
	#
	############################################
	return


def RK1(r, m, fx, h, MESA_file,masscut, *args, **kwargs):
	#need to pass the value of rho roughly at r but don't update it, just need it for calculation
	#use_unlog=bool(kwargs.get('load_unlogged',False))
	#print("in RK routine")

	h=float(h)

	rho_k0=rho_r(r,MESA_file,masscut)#,load_unlogged=use_unlog)
	rho_k12=rho_r(r+0.5*h,MESA_file,masscut)#,load_unlogged=use_unlog)
	rho_k3=rho_r(r+h,MESA_file,masscut)#,load_unlogged=use_unlog)

	k0=float(fx(r, rho_k0 )*h)
	k1=float(fx(r + 0.5*h, rho_k12 )*h)  #+ 0.5*k0
	k2=float(fx(r + 0.5*h, rho_k12 )*h)  #+ 0.5*k1
	k3=float(fx(r + h, rho_k3 )*h)  #+ k2
	r     = r     + h
	m = m + (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0
	return float(r), float(m)



def Romberg_integrate(r, m, fn,\
                      stepsize,  steps, MESA_file, masscut, debug=False, exact=None): #a, b,
	
	### change these input parameters, should not be passing m or step size; use a, b to represent limits of integral
	def local_rho(nn):
			#print("local rho activated")
		return rho_r(nn, MESA_file, masscut)

	steps=int(steps)
	table = np.zeros((steps, steps))
	pow_4 = 4 ** np.arange(steps, dtype=np.float64) - 1

	# trapezoidal rule
	a = r
	b = r + stepsize

	h = (b - a)
	table[0, 0] = h * (fn(a,  local_rho(a)) + fn(b, local_rho(b))) / 2

	for j in range(1, steps):
	    h /= 2

	    # extended trapezoidal rule
	    table[j, 0] = table[j - 1, 0] / 2
	    table[j, 0] += h * np.sum(
	        fn(a + i * h, local_rho(a + i * h)) for i in range(1, 2 ** j + 1, 2)
	    )
	    # richardson extrapolation
	    for k in range(1, j + 1):
	        table[j, k] = table[j, k - 1] + \
	            (table[j, k - 1] - table[j - 1, k - 1]) / pow_4[k]
	mass_val=table[-1, -1]
	return b, mass_val ## WARNING returning b not r


##################################################################################
def get_placement_radii(rmin, rmax, RKstep, TOL, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):
	lower_convergence_limit = 1.0 - float(TOL)
	upper_convergence_limit = 1.0 + float(TOL)	

	Romberg=kwargs.get("Romberg", False)
	input_RKstep=RKstep
	fit_region_R   =mn.MESA_r(MESA_file, masscut)
	fit_region_E   =mn.MESA_E(MESA_file, masscut) 
	
	Mshell_target=target_Mshell(force_N,mp)

	rl = rmin
	ru_mass_loop = rl # RKstep 
	Mshell_integral = 0.0	

	while ru_mass_loop <= rmax:	
		while ( (Mshell_integral/Mshell_target) <= lower_convergence_limit)\
		   or ( (Mshell_integral/Mshell_target) >= upper_convergence_limit): 

			try:		
				Mshell_integral = Mshell_from_RK(rl, ru_mass_loop, RKstep, MESA_file, masscut)
					#Mshell_integral= Mshell_from_quad(rl, ru_mass_loop, MESA_file, masscut)
				#print("Mshell_integral, ru_mass_loop, stepsize:             ",\
			 	#	  ('%.3f'%(100.0*Mshell_integral/Mshell_target)), "    ", ru_mass_loop, "   ", RKstep)
			
			except TypeError:
				print("converge_funcs ERROR: TypeError in cf.get_placement_radii()")
				break  	 		 

			### adapative step size
			if (Mshell_integral/Mshell_target) >= upper_convergence_limit:
				RKstep=RKstep-0.5*RKstep
				#print("too high, new step=", RKstep)
				ru_mass_loop = ru_mass_loop - RKstep
				Mshell_integral =0.0
			

			elif (Mshell_integral/Mshell_target) <= lower_convergence_limit:
				# experimental
				#print("WARNING! no RK doubling")
				RKstep = RKstep + RKstep   #my version 
				
				ru_mass_loop = ru_mass_loop + RKstep
				Mshell_temp=Mshell_integral
				Mshell_integral =0.0
			else:
				pass
				
		###############################################################
		#
		# reset parameters for journey to next shell
		#
		###############################################################
		## find central point between upper and lower integral limits for converged value
		r_print=(ru_mass_loop + rl)/2.0
		r_nearest,rdex=find_nearest(fit_region_R,r_print)
		u_local=fit_region_E[rdex]
		outf.write("{} {} {} {}\n".format(force_N, r_print, Mshell_integral, u_local))

		RKtemp = RKstep	
		Mshell_temp=Mshell_integral

		m_percent = '%.3f'%(100.0*Mshell_temp/Mshell_target)

		print( 	m_percent,\
				"agreement  for Mshell =",Mshell_temp,\
		       " at ", ru_mass_loop, "    ",\
		        "%.5f"%(ru_mass_loop/rmax),\
		        r"%radius  ...current RKstep=",\
		        "%1.5e"%RKtemp,\
		        " from" , "%1.5e"%input_RKstep)#, "  from ", "%1.5e"%RKtemp

		## reset
		RKstep=input_RKstep
		rl = ru_mass_loop
		Mshell_integral = 0
	return #
##################################################################################





###########################################################################
#
# load MESA data in correct format
#
###########################################################################
def get_MESA_profile_edge(MESA_file,**kwargs):

	strip=bool(kwargs.get('strip',False))
	keyword=str(kwargs.get('quantity','zone'))
	masscut=float(kwargs.get('masscut',0.65))

	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]
	
	try:
		quantity=MJ.get_quantity(MESA_file,keyword)
	except KeyError:
		print("Quantity keyword not found. Allowed keywords are:")
		print(MJ.show_allowed_MESA_keywords(MESA_file))
		sys.exit(1)

	masses = MJ.get_quantity(MESA_file,'mass').astype(np.float)
	Mtot=masses[0]

	bound = masscut*Mtot
	fit_region_indices=np.where( (masses>=bound) )[0]	
	fit_region  = quantity[fit_region_indices]

	if keyword =='mass':
		fit_region = outer_mass(Mtot, fit_region)#[Mtot-p for p in mf]

	return np.array(fit_region).astype(float)





def outer_mass(Mtot,fit_region):
	mf=fit_region
	fit_region = [Mtot-p for p in mf]
	#print("fit_region in outer_mass: ", fit_region)
	fit_region=np.array(fit_region)

	select=np.where( (fit_region>0.0) )[0]
	nfit_region=fit_region[select]

	return np.array(nfit_region).astype(float)



###########################################################################
#
# HEALPix functions
#
###########################################################################
def healpixify(N):
	NSIDE=int(N)
	#this is just a straight up array of particle IDs
	ipix_array=np.arange(hp.nside2npix(NSIDE)) 
	x=[]
	y=[]
	z=[]
	for i in range(len(ipix_array)):
	    ipix=ipix_array[i]
	    coord=hp.pixelfunc.pix2vec(NSIDE, ipix, nest=True)

	    x.append(coord[0])
	    y.append(coord[1])
	    z.append(coord[2])
	return x,y,z


def get_coords(N):
	x,y,z=healpixify(N)

	theta=random_theta()
	## update mjoyce 10/31/18
	phi=random_theta()
	psi=random_theta()

	# these are already in unit by the time they get here
	xd, yd, zd = rotate_shell(x,y,z,theta,"about_z")
	xe, ye, ze = rotate_shell(xd,yd,zd,phi,"about_y")
	xf, yf, zf = rotate_shell(xe,ye,ze,psi,"about_x")

	xf = to_array(xf)
	yf = to_array(yf)
	zf = to_array(zf)

	return xf.flatten(), yf.flatten(), zf.flatten()



def rotate_shell(x_array, y_array, z_array, theta, direction, **kwargs):
    vec=np.array( [x_array, y_array, z_array])

    Rx=np.matrix( [\
    [1.0, 0.0, 0.0],\
    [0, np.cos(theta), -np.sin(theta)],\
    [0, np.sin(theta), np.cos(theta)]\
    ])

    Ry=np.matrix( [\
    [np.cos(theta), 0.0, np.sin(theta)],\
    [0.0, 1.0, 0.0],\
    [-np.sin(theta), 0.0, np.cos(theta)]\
    ])

    Rz=np.matrix( [\
    [np.cos(theta), -np.sin(theta), 0.0],\
    [np.sin(theta), np.cos(theta), 0.0],\
    [0.0, 0.0, 1.0]\
    ])

    if str(direction) == "about_z":
        new=Rz*vec
    elif str(direction) == "about_y":
        new=Ry*vec
    else:
        new=Rx*vec

    new_x=np.array(new[0].transpose()).astype(float)
    new_y=np.array(new[1].transpose()).astype(float)
    new_z=np.array(new[2].transpose()).astype(float)

    return new_x, new_y, new_z


###########################################################################
#
# basic math functions
#
###########################################################################
def R_to_solar_f(r):
	return r*R_to_solar

def to_Rsun(R_in_cm):
	try:
		Rs=np.float(R_in_cm)/R_to_solar
		return Rs
	except:
		Rs=[]
		for  i in range(len(R_in_cm)):
			Rs.append(float(R_in_cm[i])/R_to_solar)
		Rs= np.array(Rs)
		return Rs 

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx


def to_log(xq):
	logged=[]
	for i in range(len(xq)):
		try:
			q=np.log10(xq[i])
			logged.append(q)
		except TypeError:
			pass
	logged = np.array(logged).astype(float)
	return logged 


def unlog(xq):
	try:
		unlogged=10.0**xq
	except TypeError:
		unlogged=10.0**float(xq)
	except:
		print("error in unlog (MESA2HYDRO/lib/converge_funcs.py")
		sys.exit()
	return unlogged 


def random_theta():
	#.random gives random float between 0 and 1
    theta=rand.random()*2.0*np.pi 
    return theta


def to_rad(theta):
    theta=theta*np.pi/180.0
    return theta


def to_array(array):
    return np.array(array).astype(np.float)


def min_index(array):
	mn,idx = min( (array[i],i) for i in xrange(len(array)) )
	return idx

def volume(r):
	r=float(r)
	vol=(4.0/3.0)*np.pi*r**3.0
	return float(vol)


###########################################################################
#
# Recovery tools
#
###########################################################################
def one_over_r(xdata,A,B,C,D):
	return A*( 1.0/( (D*xdata)-B) ) + C

def get_curve(r_array,rho_array,guess_a,guess_b,guess_c,guess_d,functional_form):
	import scipy.optimize 
	p=[guess_a,guess_b,guess_c, guess_d]
	params, cov = curve_fit(functional_form, r_array, rho_array, p0=p)
	a, b, c,d = params

	return a,b,c,d


def poly_curve(xdata,ydata,degree):
	import scipy.optimize 
	p=np.polyfit(xdata,ydata,degree)
	d=degree
	if d==2:
		y=p[0]*(xdata**d) + p[1]*(xdata**(d-1)) +  p[d] 
	elif d==3:
		y=p[0]*(xdata**d) + p[1]*(xdata**(d-1))+ p[2]*(xdata**(d-2)) +  p[d]
	elif d==4:
		y=p[0]*(xdata**d) + p[1]*(xdata**(d-1))+ p[2]*(xdata**(d-2))+ p[3]*(xdata**(d-3)) +  p[d]
	elif d==5:
		y=p[0]*(xdata**d) + p[1]*(xdata**(d-1))+ p[2]*(xdata**(d-2))\
		+ p[3]*(xdata**(d-3))+ p[4]*(xdata**(d-4)) +  p[d]
	else:
		print('fits with degree > 5 not supported')
		sys.exit()
	return y

# end module converge_funcs
