#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
#import pygadgetreader as pgr # works- credit this person
import MESAlibjoyce as MJ
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
	use_unlog=bool(kwargs.get('load_unlogged',False))

	rho_k0=rho_r(r,MESA_file,masscut,load_unlogged=use_unlog)
	#print ">>>>>>>>>>>>>>>>>r, m, rho in start RK:", r, m, rho_k0#, ('%1.5e'%h)
	rho_k12=rho_r(r+0.5*h,MESA_file,masscut,load_unlogged=use_unlog)
	rho_k3=rho_r(r+h,MESA_file,masscut,load_unlogged=use_unlog)

	k0=fx(r, rho_k0 )*h
	#print "k0:", k0
	k1=fx(r + 0.5*h, rho_k12 )*h  #+ 0.5*k0
	k2=fx(r + 0.5*h, rho_k12 )*h  #+ 0.5*k1
	k3=fx(r + h, rho_k3 )*h  #+ k2
	r     = r     + h
	m = m + (k0 + 2.0*k1 + 2.0*k2 + k3)/6.0
	#print ">>>>>>>>>>>>>>>>>r, m, rho in end RK:", r, m, rho_k3
	return r, m

def density_integral_numeric(r, rho_r): 
	# will need to write an interpolation function to go between points in order to have
	# a rho_r value for any possible r
	dmdr=4.0*np.pi*r**2.0*rho_r
	#print "value of r, rho_r, dmdr:", r, rho_r, dmdr
	return dmdr


def Mshell_from_RK(rl, rmax, step, MESA_file,masscut, *args,**kwargs):
	use_unlog=bool(kwargs.get('load_unlogged',False))
	#print "value for use_unlog sent to Mshell_from_RK:", use_unlog
	### rmax is NOT MODIFIED BY THIS ROUTINE, it is simply a STOPPING CRITERION
	r = rl      
	m=0
	while r<rmax:
	    r, m= RK1(r,m, density_integral_numeric, step, MESA_file,masscut, load_unlogged=use_unlog) #rho_r(r,MESA_file,masscut,load_unlogged=use_unlog)
	    #print "r being used in Mshell_from_RK:",r,m
	    #sys.exit()
	Mshell=m
	return Mshell#, r#, rho_r

def find_nearest(array,value):
	# this isn't quite right, use bisect??
    idx = (np.abs(array-value)).argmin()
    return array[idx],idx



def rho_r(r,MESA_file,masscut, *args, **kwargs):
	############################################################
	#
	# WARNING! FIXING cgs units!!!
	#
	############################################################
	use_unlog=bool(kwargs.get('load_unlogged',False))
	#print "value of use_unlog in rho_r:", use_unlog


	fit_region_R = get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
	fit_region_R=unlog(fit_region_R)

	if use_unlog:
		fit_region_rho = get_MESA_profile_edge(MESA_file, quantity='rho', masscut=masscut ,strip=False)
	else:
		fit_region_rho = get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
		fit_region_rho=unlog(fit_region_rho)

	## convert unlogged radii from Rsolar units to cgs
	fit_region_R=fit_region_R*R_to_solar
	#print "r, fit_region_R after conversion:", r, fit_region_R[-1]

	r0,idx=find_nearest(fit_region_R,r)
	#
	# WARNING! THIS RELIES ON LOADED DATA BEING SORTED! DO NOT TAMPER!
	#

	if r0 <= r:
		rho0=fit_region_rho[idx]
		r1=fit_region_R[idx-1]
		rho1=fit_region_rho[idx-1]
	else:
		r0=fit_region_R[idx+1]
		rho0=fit_region_rho[idx+1]
		r1=fit_region_R[idx]
		rho1=fit_region_rho[idx]
	#np.polyfit()

	rrho_r= (  (r1-r)*rho0 + (r-r0)*rho1 ) /(r1-r0)
	#print "rrho_r in rho_r(r) routine:", rrho_r

	#print idx, r0,rho0, r1, rho1
	if (r0 <= r <= r1):
		#print r0, r, r1
		return rrho_r
	else:
		print "no."
		return #rho_r



def target_Mshell(N,mp):
	Mshell=N**2.0*12.0*mp
	return Mshell


def get_placement_radii(rl, ru, RKstep, force_N, mp, MESA_file, masscut, *args, **kwargs):
	############### FINE
	use_unlog=bool(kwargs.get('load_unlogged',True))
	SO=bool(kwargs.get('suppress_output',False))

	##print "value of use_unlog:", use_unlog

	Mshell_target=target_Mshell(force_N,mp)
	#print "target_Mshell", Mshell_target
	#sys.exit()

	Mshell=0
	Mshell_temp=0
	oldru=rl
	while Mshell <= Mshell_target:
		Mshell=Mshell_temp+Mshell_from_RK(oldru, ru, RKstep, MESA_file,masscut, load_unlogged=use_unlog)
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


####################################################################
#
# functions for validating density profile recovered from GADGET IC file
#
#####################################################################
# def get_density_from_loaded_particles(r, r_bin, r_recovered, p_mass):
# 	##unsure about r and r_bin in these
# 	r=float(r)
# 	r_bin=float(r_bin)
# 	mshell=Mshell_r(r,r_bin,r_recovered, p_mass)
# 	vol_annulus=volume(r+r_bin)-volume(r)
# 	rho = float(mshell)/float(vol_annulus)
# 	#print "   rho", rho
# 	#(1.0/(4.0*np.pi * (r+r_bin)**3.0 )) #-( 3.0/(4.0*np.pi * (r)**3.0 ))
# 	return rho



##############################################
#
# THIS FUNCTION IS THE PROBLEM 
#
# #################################################
# def Mshell_r(r, r_bin,r_recovered,p_mass):
# 	######################################
# 	# warning! requires a fixed particle mass, can't handle changing m
# 	######################################
# 	# total number of particles in the binned region'
# 	region = np.where( (r_recovered>r) & (r_recovered <(r+r_bin) )  )
# 	N=len(r_recovered[region])
# 	mshell=p_mass*N
# 	#print 'recovered_region: ',r_recovered[region]
 

# 	#print "\nr",r,"r_bin",r_bin,'   len(region)=N', N,'      mshell',mshell	
# 	return mshell





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
	#vol=4.0*np.pi*r**2.0 ## this is actually the surface area but fuck it
	return float(vol)
