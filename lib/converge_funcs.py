#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
#import pygadgetreader as pgr # works- credit this person
import MESAhandling as MJ
import mainlib as mn

import datetime as dt 
import random as rand
import healpy as hp



M_to_solar=1.988e33 #*10.0**33.0 ## g/Msolar
R_to_solar=6.957e10 #*10.0**10.0 ## cm/Rsolar

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


#######################################################################################
#
# data scanning functions
#
#######################################################################################
def rho_r(r, MESA_file, masscut, *args, **kwargs):
	############################################################
	#
	# WARNING! FIXING cgs units!!!
	#
	############################################################
	

	#### REEVALUATE MASSCUT BASED ON NEW R
	# masscut = mn.masscut_from_r(r, MESA_file)


	fit_region_R=mn.MESA_r(MESA_file,masscut)

	# #print "inside rho_r: ", fit_region_R.max()

	fit_region_rho=mn.MESA_rho(MESA_file,masscut)
	r0,idx=find_nearest(fit_region_R,r)
	# WARNING! THIS RELIES ON LOADED DATA BEING SORTED! DO NOT TAMPER!

	#try:
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
		print "WARNING! r not found between r0 and r1....end of MESA density values"
		#print "WARNING! rrho forced to rho_1"
		return  #float(rho1)
	# except IndexError:
	# 	rrho_r= fit_region_rho[idx]
	# 	#print "WARNING! index error in try/catch block cf.rho_r()"
	# 	return float(rrho_r)		

		#pass
		#return 
		
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
		print "could not compute m(r) for point"
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
	    #print "r, m: ", r, m
	Mshell=m
	return Mshell


from scipy import integrate
from scipy.integrate import romberg, quad 
def Mshell_from_Romberg(rl, rmax, MESA_file, masscut,  *args,**kwargs): #step, MESA_file,masscut,
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
	#print "Mshell from quad: ", Mshell
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
	#print "in RK routine"

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
			#print "local rho activated"
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

	    #print "(Romberg location 4) step j=", j

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
	#print "(Rombger location 5) ", r, mass_val

	### should not be returing b! onyl mass_val
	return b, mass_val ## WARNING returning b not r

##################################################################################

# def get_placement_radii(rl, rmax, RKstep, TOL, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):
# 	#print "loc 1"
# 	#scale the tolerance by the maximum stellar mass--- units must be preserved
# 	# tolerance = float(TOL)*(MJ.get_quantity(MESA_file,'mass').astype(np.float)*M_to_solar).max() ##(1.0 - float(TOL))
# 	fit_region_R   =mn.MESA_r(MESA_file, masscut)
# 	fit_region_E   =mn.MESA_E(MESA_file, masscut) 
# 	rmin = fit_region_R.min()
# 	rmax = fit_region_R.max()

# 	Mshell_target = target_Mshell(force_N, mp)

# 	#tolerance = float(TOL)*Mshell_target  ## random tolerance adjustment
# 	tolerance =TOL*Mshell_target

# 	## initial guess
# 	rl = rmin
# 	ru_mass_loop = rl + (rmax-rl)/1000.0 #rl + RKstep
# 	while ru_mass_loop <= rmax:	

# 		def f(r_slide):
# 			#print "loc 5"
# 			# Mshell_target = target_Mshell(force_N, mp)
# 			Mshell_integral = Mshell_from_RK(rl, r_slide, RKstep, MESA_file, masscut)

# 			quantity= Mshell_integral - Mshell_target
# 			#print "rl, ru: ",rl/rmax, r_slide/rmax,"   value of Mshell_from_RK: ", Mshell_integral, "   Mshell_target: ", Mshell_target,  
# 			#print "         mass difference: ", quantity
# 			return  quantity

# 		def dfdr(r_slide):
# 			#return 4.0*np.pi*r_slide**2.0*rho_r(r_slide, MESA_file, masscut)
# 			#print "loc 6"
# 			return density_integral_numeric(r_slide,rho_r(r_slide, MESA_file, masscut))

# 		def dx(f, x):
# 			#print "loc 4"
# 			return np.abs(0.0-f(x))

# 		def Newton(f, dfdr, r0, e):
# 			#print "loc 3"
# 			delta = dx(f,r0)

# 			#print "delta, e: ", delta, e
			
# 			while delta > e:
# 			## MODIFIED!
# 				#if r0 <= rmin:
# 					#print "attempted to use ru < rmin...refining guess for ru..."
# 					#ru_mass_loop = 0.5*ru_mass_loop
# 					#r0 = 2.0*r0
# 					#r0 = rmin + 1e-9*rmin
# 					#r0 = r0 - f(r0)/dfdr(r0)
# 					#delta = dx(f, r0)
# 					#return rmin ##go back to the minimum in event of boundary crossing
# 				#else:
# 				r0 = r0 - f(r0)/dfdr(r0)
# 				print "(Newton) r=", r0/fit_region_R.max()
# 				delta = dx(f, r0)
# 			#print "convergence found at: ", r0	
# 			return r0

# 		##print "loc 2"
# 		#try:
# 		ru_found = Newton(f, dfdr, ru_mass_loop, tolerance)
# 		# except IndexError:
# 		# 	print "attempted to use ru < rmin...refining guess for ru..."
# 		# 	ru_mass_loop = 0.5*ru_mass_loop
# 		# 	print "changing value of ru.."
# 		# 	break



# 		#print "convergence to within ", TOL, " found at r=", "%.5f"%(ru_found/rmax)

# 		## find central point between upper and lower integral limits for converged value
# 		Mshell_integral_at_ru= Mshell_from_RK(rl, ru_found, RKstep, MESA_file, masscut)
		
# 		r_print=(ru_found + rl)/2.0
# 		r_nearest,rdex=find_nearest(fit_region_R,r_print)
# 		u_local=fit_region_E[rdex]
		
# 		#print >> outf, force_N, r_print, Mshell_integral_at_ru, u_local


# 		print "Mshell agreement ",  ('%.3f'%(100.0*Mshell_integral_at_ru/Mshell_target))\
# 		       ," at ", "%.5f"%(ru_mass_loop/rmax), r"%radius  ...current RKstep=","%1.5e"%RKstep#

# 		rl = ru_found
# 		ru_mass_loop = ru_found + RKstep

		
# 	return 


###################################################################################


Romberg=False
def get_placement_radii(rl, rmax, RKstep, TOL, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):
	lower_convergence_limit = 1.0 - float(TOL)
	upper_convergence_limit = 1.0 + float(TOL)	

	Romberg=kwargs.get("Romberg", False)
	input_RKstep=RKstep
	#rtot=(mn.MESA_r(MESA_file, masscut)).max()
	fit_region_R   =mn.MESA_r(MESA_file, masscut)
	fit_region_E   =mn.MESA_E(MESA_file, masscut) 
	
	#rmax=ru#fit_region_R.max()
	Mshell_target=target_Mshell(force_N,mp)
	Mshell=0
	#rdiff= rmax

	#rl=fit_region_R.min()
	temp=rl + RKstep
	ru_mass_loop = temp
	Mshell_integral = 0.0	
	while ru_mass_loop <= rmax:	
		Mshell_target = 12.0*force_N**2.0*mp

		#print "loc 1 value of RKstep: ", RKstep
		while ( (Mshell_integral/Mshell_target) <= lower_convergence_limit)\
		   or ( (Mshell_integral/Mshell_target) >= upper_convergence_limit): 

			try:		
				if Romberg:
					#print "WARNING!!! ROMBERG SWITCH ON!"
					Mshell_integral = Mshell_from_Romberg(rl, ru_mass_loop, MESA_file, masscut, steps=1)
					#print "loc 2"
				else:
					Mshell_integral = Mshell_from_RK(rl, ru_mass_loop, RKstep, MESA_file, masscut)#, load_unlogged=use_unlog)
					#Mshell_integral= Mshell_from_quad(rl, ru_mass_loop, MESA_file, masscut)
				#print "Mshell_integral, ru_mass_loop, stepsize:             ",\
			 	#	  ('%.3f'%(100.0*Mshell_integral/Mshell_target)), "    ", ru_mass_loop, "   ", RKstep
			
			except TypeError:
				print "TypeError in cf.get_placement_radii()"
				break  	 		 

			### adapative step size???
			#print "ru before reset: ", ru_mass_loop
			if (Mshell_integral/Mshell_target) >= upper_convergence_limit:
				RKstep=RKstep-0.5*RKstep
				#print "too high, new step=", RKstep
				ru_mass_loop = ru_mass_loop - RKstep
				Mshell_integral =0.0
			
			elif (Mshell_integral/Mshell_target) <= lower_convergence_limit:
				#
				# experimental
				#print "too low, new step=", RKstep

				#RKstep=RKstep-0.5*RKstep  #john said to do this but I don't buy it 
				RKstep = RKstep + RKstep   #my version 
				
				ru_mass_loop = ru_mass_loop + RKstep
				Mshell_temp=Mshell_integral
				Mshell_integral =0.0
			else:
				pass
				

		## find central point between upper and lower integral limits for converged value
		r_print=(ru_mass_loop + rl)/2.0
		r_nearest,rdex=find_nearest(fit_region_R,r_print)
		u_local=fit_region_E[rdex]
		print >> outf, force_N, r_print, Mshell_integral, u_local

		RKtemp = RKstep	
		Mshell_temp=Mshell_integral
		print "Mshell agreement ",  ('%.3f'%(100.0*Mshell_temp/Mshell_target))\
		       ," at ", "%.5f"%(ru_mass_loop/rmax), r"%radius  ...current RKstep=","%1.5e"%RKtemp,\
		        "%1.5e"%input_RKstep,"  from ", "%1.5e"%RKtemp

		## reset
		print "not resetting!"
		#RKstep=input_RKstep

		#RKstep=RKstep + RKstep ## doubling for basically no reason
		rl = ru_mass_loop
		ru_mass_loop = ru_mass_loop + RKstep
		Mshell_integral = 0

	return #



###########################################################################
#
# load MESA data in correct format
#
###########################################################################

#############################l######################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################


def get_MESA_profile_edge(MESA_file,**kwargs):#, strip):
	#print MJ.show_allowed_MESA_keywords(MESA_file)
	strip=bool(kwargs.get('strip',False))
	keyword=str(kwargs.get('quantity','zone'))
	masscut=float(kwargs.get('masscut',0.65))

	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]
	
	try:
		quantity=MJ.get_quantity(MESA_file,keyword)
	except:
		print "Quantity keyword not found. Allowed keywords are:"
		print MJ.show_allowed_MESA_keywords(MESA_file)
		sys.exit()

	masses = MJ.get_quantity(MESA_file,'mass').astype(np.float)#*M_to_solar ## WARNING
	Mtot=masses[0]

	## this is OK-- non-solar units
	#print "in get_MESA_profile_edge, Mtot=", Mtot

	bound = masscut*Mtot
	fit_region_indices=np.where( (masses>=bound) )[0]	
	fit_region  = quantity[fit_region_indices]

	if keyword =='mass':
		fit_region = outer_mass(Mtot, fit_region)#[Mtot-p for p in mf]

	## this is OK-- non-solar units
	#print "in get_MESA_profile_edge, fit_region=", fit_region

	return np.array(fit_region).astype(float)



def outer_mass(Mtot,fit_region):
	mf=fit_region
	fit_region = [Mtot-p for p in mf]
	#print "fit_region in outer_mass: ", fit_region
	fit_region=np.array(fit_region)

	select=np.where( (fit_region>0.0) )[0]
	nfit_region=fit_region[select]

	return np.array(nfit_region).astype(float)



###########################################################################
#
# healpix
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
# basic
#
###########################################################################
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
	#print type(xq)
	try:
		unlogged=10.0**xq
	except TypeError:
		unlogged=10.0**float(xq)
	except:
		print "error in unlog (MESA2HYDRO/lib/converge_funcs.py"
		sys.exit()
	return unlogged 


def random_theta():
    theta=rand.random()*2.0*np.pi #.random gives random float between 0 and 1
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
# Recovery
#
###########################################################################
import scipy.optimize 

def one_over_r(xdata,A,B,C,D):
	#print "xdata: ", xdata
	return A*( 1.0/( (D*xdata)-B) ) + C


#------------------------------ analytic ----------------------------------------

def get_curve(r_array,rho_array,guess_a,guess_b,guess_c,guess_d,functional_form):

	p=[guess_a,guess_b,guess_c, guess_d]
	params, cov = curve_fit(functional_form, r_array, rho_array, p0=p)
	a, b, c,d = params

	return a,b,c,d#,d#,y 


def poly_curve(xdata,ydata,degree):
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
		print 'fits with degree > 5 not supported'
		sys.exit()
	return y
