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



M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar


#######################################################################################
#
# Numerical integration
#
#######################################################################################

def RK1(r, m, fx, h,MESA_file,masscut, *args, **kwargs):
	#need to pass the value of rho roughly at r but don't update it, just need it for calculation
	#use_unlog=bool(kwargs.get('load_unlogged',False))
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

def density_integral_numeric(r, rho_r): 
	# will need to write an interpolation function to go between points in order to have
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


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx



def rho_r(r,MESA_file,masscut, *args, **kwargs):
	############################################################
	#
	# WARNING! FIXING cgs units!!!
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
		print "no."
		return 
		
def m_r(r,MESA_file,masscut, *args, **kwargs):
	############################################################
	#
	# WARNING! FIXING cgs units!!!
	#
	############################################################
	fit_region_R=mn.MESA_r(MESA_file,masscut)
	fit_region_m=mn.MESA_m(MESA_file,masscut)
	r0,idx=find_nearest(fit_region_R,r)
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
		print "no."
		return 



def get_logE(r,MESA_file,masscut):

	fit_region_R=mn.MESA_r(MESA_file,masscut)
	fit_region_E=mn.MESA_internalE(MESA_file,masscut)

	r0,idx=find_nearest(fit_region_R,r)

	if r0 <= r:
		E0=fit_region_E[idx]
		r1=fit_region_R[idx-1]
		E1=fit_region_E[idx-1]
	else:
		r0=fit_region_R[idx+1]
		E0=fit_region_E[idx+1]
		r1=fit_region_R[idx]
		E1=fit_region_E[idx]

	E= (  (r1-r)*E0 + (r-r0)*E1 ) /(r1-r0)
	if (r0 <= r <= r1):
		return E
	else:
		print "Error obtaining internal energy E."
		return 



def target_Mshell(N,mp):
	Mshell=N**2.0*12.0*mp
	return Mshell


def get_placement_radii(rl, ru, RKstep, force_N, mp, MESA_file, masscut, *args, **kwargs):
	############### FINE
	#use_unlog=bool(kwargs.get('load_unlogged',True))
	SO=bool(kwargs.get('suppress_output',False))


	Mshell_target=target_Mshell(force_N,mp)
	Mshell=0
	Mshell_temp=0
	oldru=rl
	while Mshell <= Mshell_target:
		Mshell=Mshell_temp+Mshell_from_RK(oldru, ru, RKstep, MESA_file,masscut)#, load_unlogged=use_unlog)
		oldru=ru
		Mshell_temp=Mshell
		ru=ru+RKstep

	if SO:
		pass
	else:
		print "target_Mshell",('%.3E'%Mshell_target)," convergence at Mshell:", ('%.3E'%Mshell),\
		 '   ', ('%.3e'%(Mshell/Mshell_target)),r'x target',\
		"  rl, ru:", ('%1.3e'%rl),('%1.3e'%ru), "  RKstep: ",  ('%1.3e'%RKstep)
	r_place=(ru+rl)/2.0
	return r_place,Mshell



###########################################################################
#
# load MESA data in correct format
#
###########################################################################
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

	masses = MJ.get_quantity(MESA_file,'mass').astype(np.float)
	Mtot=masses[0]
	bound = masscut*Mtot
	fit_region_indices=np.where( (masses>=bound) )[0]	
	fit_region  = quantity[fit_region_indices]

	#print "keyword: ",keyword
	if keyword =='mass':
		#xq_unlogged = [10.0**p for p in xq]
		#mf = fit_region
		fit_region = outer_mass(Mtot, fit_region)#[Mtot-p for p in mf]

	return np.array(fit_region).astype(float)



def outer_mass(Mtot,fit_region):
	mf=fit_region
	fit_region = [Mtot-p for p in mf]
	return np.array(fit_region).astype(float)



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
	# these are already in unit by the time they get here
	xd, yd, zd = rotate_shell(x,y,z,theta,"about_z")
	xe, ye, ze = rotate_shell(xd,yd,zd,theta,"about_y")
	xf, yf, zf = rotate_shell(xe,ye,ze,theta,"about_x")

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
	unlogged=[]
	for i in range(len(xq)):
		#print xq
		try:
			q=10.0**float(xq[i])
			unlogged.append(q)
		except ValueError:
			pass
	unlogged = np.array(unlogged).astype(float)
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

# #------------------------------------------------------------------------------------
# def density_integral(rl,ru, A,B,C):
# 	rdiff=(ru - rl)
# 	f1=0.5*A*rdiff**2.0  + A*B*rdiff  + A*B**2.0 * np.log10( abs(rdiff - B) ) + (1.0/3.0)*C*rdiff**3.0
# 	return f1

# ### is it necessary to do this density integral or should I just fit the mass profile and pass that?

# def Mshell_from_integral(rl, ru, A,B,C):
# 	Mshell=4.0*np.pi*density_integral(rl, ru, A, B, C)
# 	return Mshell


# def approximate_mp(rl, ru, n_p, A, B, C):
# 	Mshell=4.0*np.pi*density_integral(ru,rl,A,B,C)
# 	mp=Mshell/n_p
# 	return mp