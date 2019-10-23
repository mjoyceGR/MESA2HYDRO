#!/usr/bin/env python

#************************************************************************
#
#   Copyright (C) 2019  M. Joyce, L. Lairmore, D. J. Price
#
#   See MESA2HYDRO/LICENSE
#
#************************************************************************

'''

Contains: main subroutines for run_conversion.py 

'''

import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import sys 
import hdf5lib as hdf5lib
import time

import converge_funcs as cf
import MESAhandling as MJ
import io_lib as rw
import constants as const

const.logo()

M_to_solar=const.Msun 
R_to_solar=const.Rsun 

### FOR DEBUGGING ONLY
use_normalized=False


def masscut_from_r(r, MESA_file):
    fit_region_R = (10.0**(MJ.get_quantity(MESA_file, 'logR').astype(np.float)))*R_to_solar 
    fit_region_M = MJ.get_quantity(MESA_file,'mass').astype(np.float)*M_to_solar 
    r_bound = r 
    assoc_region=np.where( fit_region_R >= r_bound)[0]
    mass_assoc = fit_region_M[assoc_region]
    masscut = mass_assoc.min()/mass_assoc.max()
    return masscut


def MESA_r(MESA_file, masscut):
    ## returns R in cm
    fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
    fit_region_R=cf.unlog(fit_region_R)*R_to_solar
    return fit_region_R


def MESA_rho(MESA_file,masscut):
    ## returns rho in cgs
    fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
    fit_region_rho=cf.unlog(fit_region_rho)
    return fit_region_rho

def MESA_m(MESA_file,masscut):
    ### returns mass in g

    ## WARNING!!!!
    ## THIS MIGHT BE DEFINITELY WRONG! 7/30/19

    fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)*M_to_solar
    return fit_region_M


def MESA_P(MESA_file, masscut):
    ## returns pressure in cgs
    fit_region_P=cf.get_MESA_profile_edge(MESA_file, quantity='logP', masscut=masscut, strip=False)
    fit_region_P=cf.unlog(fit_region_P)
    return fit_region_P

def MESA_E(MESA_file, masscut):
    ## cgs
    fit_region_E=cf.get_MESA_profile_edge(MESA_file, quantity='logE', masscut=masscut, strip=False)
    fit_region_E=cf.unlog(fit_region_E)
    return fit_region_E


def central_mass(MESA_file, masscut):
    ## define the mass to be cast as the sink particle/gravity well/stellar core
    masses = MJ.get_quantity(MESA_file,'mass').astype(np.float)*M_to_solar 
    radii = 10.0**(MJ.get_quantity(MESA_file, 'logR').astype(np.float))
    Mtot=masses.max()
    central_M = masscut*Mtot
    return central_M, Mtot


def check_MESA(MESA_file, masscut,uselog=True,save=True):
    fit_region_R=MESA_r(MESA_file,masscut)
    fit_region_rho=MESA_rho(MESA_file,masscut)
    import matplotlib.pyplot as plt
   
    if uselog:
        fit_region_R=cf.to_log(fit_region_R)
        fit_region_rho=cf.to_log(fit_region_rho)

    plt.plot(fit_region_R, fit_region_rho,'b.')
    plt.xlabel('Radius (cm)')
    plt.ylabel(r'$\rho$ (g/cm$^3$)')
    if uselog:
        plt.xlabel('log R (cm)')
        plt.ylabel(r'log $\rho$ (g/cm$^3$)')
    plt.title(MESA_file)
    if save:
        plt.savefig('latest_loaded_MESA_profile.png')
    plt.show()
    plt.close()

def get_sink_mass(MESA_file, masscut):
    sink_mass=central_mass(MESA_file, masscut)[0]
    return sink_mass


##########################################################
#
# Make an NR file from MESA data
#
###########################################################
def make_NR_file(MESA_file,masscut,N,mp, RKstep, TOL, NR_file, *args, **kwargs):
    Romberg=kwargs.get("Romberg", False)    

    start_time=time.time()
    fit_region_R   =MESA_r(MESA_file, masscut)

    rmin=fit_region_R.min()
    rmax=fit_region_R.max()

    outf=NR_file
    outf.write('## fname',  MESA_file ,\
               ' masscut',  masscut,\
               '   N',      N,\
               '  mp (Ms)', mp/M_to_solar,\
               '  mp (g)',  mp,\
               '  stepsize',  ('%1.7e'% RKstep), "\n")
    
    outf.write('#N    (ru+rl)/2 (cm)    Mcontained in shell ru-rl     u(rmid)\n')


    cf.get_placement_radii(rmin, rmax, RKstep, TOL,  N, mp, MESA_file,masscut, outf)

    print('runtime: ', "%.1f"%(time.time()-start_time), "seconds")
    outf.write('#\n#\n# runtime: ', time.time()-start_time, " seconds\n")


##########################################################
#
# Make the initial conditions file from an NR file
#
###########################################################
def get_IC(MESA_file, masscut, NR_file_name,output_filename,mp,which_dtype='f', *args, **kwargs): #temp remove rmax
    filetype=str(kwargs.get('format_type','binary'))
    lognorm=kwargs.get('lognorm',False)
    phantom_rescale=kwargs.get('phantom_rescale',True)

    ## do not actually need to load E here; better to search MESA profile directly 
    N,rmid=np.loadtxt(NR_file_name, usecols=(0,1), unpack=True) 

    super_x=[]
    super_y=[]
    super_z=[]

    fit_region_rho=MESA_rho(MESA_file,masscut)
    fit_region_R=MESA_r(MESA_file, masscut)
    fit_region_P=MESA_P(MESA_file, masscut)
    fit_region_E=MESA_E(MESA_file, masscut)
    super_rho=[]
    super_P=[]
    super_E=[]
    
    for i in range(len(N)):
        NSIDE= N[i]
        r_mid= rmid[i]

        r_nearest,r_idx=cf.find_nearest(fit_region_R,r_mid)
        local_MESA_rho=fit_region_rho[r_idx]
        local_MESA_P=fit_region_P[r_idx]
        local_MESA_E=fit_region_E[r_idx]

        radius=float(rmid[i])
        if use_normalized:
            print("WARNING! unit changed from [radius (cm)] to [radius/Rsolar (cm)]")
            radius = radius/R_to_solar
        else:
            radius=radius

        ###############################################################
        #
        # Arbitrarily rotate shells
        #
        ###############################################################
        x,y,z=cf.get_coords(NSIDE)
        x=x*radius
        y=y*radius
        z=z*radius
 
        rho_array = x*0 + local_MESA_rho
        P_array   = x*0 + local_MESA_P
        E_array   = x*0 + local_MESA_E

        super_x  =np.concatenate((super_x,x),axis=0)
        super_y  =np.concatenate((super_y,y),axis=0)
        super_z  =np.concatenate((super_z,z),axis=0)

        super_rho =np.concatenate((super_rho,rho_array), axis=0) 
        super_P   =np.concatenate((super_P, P_array), axis=0)
        super_E   =np.concatenate((super_E, E_array), axis=0)


    ############################################################    
    #
    # add in the central mass point at coordiante 0,0,0
    #
    ############################################################
    zero_space=np.array([0.0])

    super_x=np.concatenate((super_x,zero_space),axis=0)
    super_y=np.concatenate((super_y,zero_space),axis=0)
    super_z=np.concatenate((super_z,zero_space),axis=0)
    
    central_point_mass, Mstar=central_mass(MESA_file, masscut)

    super_x=cf.to_array(super_x)
    super_y=cf.to_array(super_y)
    super_z=cf.to_array(super_z)

    super_rho=cf.to_array(super_rho)
    super_P= cf.to_array(super_P)
    super_E=cf.to_array(super_E)


    ############################################################    
    #
    # Renormalization of physical units
    #
    ############################################################
    if phantom_rescale:
        print("IC WARNING: Renormalization of [mp] to [(Mstar - Mcore)/Np]\n")
        mp = (Mstar-central_point_mass)/len(super_x) 
    else:
        print("IC WARNING: Renormalization of mp is OFF\n")

    if use_normalized:
        print("IC WARNING: unit changed from [mass (g)] to [mass/Msolar (g)]\n")
        mp = mp/M_to_solar
        central_point_mass= central_point_mass/M_to_solar
    else:
        mp = mp
        central_point_mass= central_point_mass

    #############################################
    #
    # 5/22/19
    #
    # WARNING! 
    # only gadget_binary writer is up to date!!!!
    #
    # 5/22/19
    #
    #############################################
    if filetype=='hdf5':
        var=rw.make_IC_hdf5(str(output_filename)+ '.hdf5',\
        mp, central_point_mass,\
        super_x, super_y, super_z,super_E) #, userho=False
        #svar=rw.make_IC_hdf5_old_way(hdf5file, mp, super_x, super_y, super_z,super_E)


    elif filetype=='phantom_binary':
        if lognorm:
            print("Phantom IC WARNING: log normalization rho --> log10(rho) for phantom_binary format IC is ON!\n")
            logged_super_rho = np.log10(np.array(super_rho))

            var = rw.make_IC_Phantom(str(output_filename),\
                    mp, central_point_mass,\
                    super_x, super_y, super_z,\
                    logged_super_rho, super_P, super_E, lognorm=True)
        else:
            var = rw.make_IC_Phantom(str(output_filename),\
                    mp, central_point_mass,\
                    super_x, super_y, super_z,\
                    super_rho, super_P, super_E, lognorm=False)
    
        import subprocess                    
        subprocess.call("mv ../work/star_00000.tmp  " + str(var) +"_00000.tmp", shell=True)  


    elif filetype=='gadget_binary':
        var=rw.make_IC_binary(str(output_filename)+ '.bin',\
        mp, central_point_mass,\
        super_x, super_y, super_z,\
        super_rho, super_P, super_E,
        which_dtype=which_dtype)  #central mass not handled (??? is this true)

    else:
        var=rw.make_IC_text(str(output_filename)+ '.txt',\
        (super_x*0. + mp), central_point_mass,\
        super_x, super_y, super_z,\
        super_rho, super_P, super_E,\
        which_dtype=which_dtype)


    print(var, type(var))
    return



#############################################
#
# plotting routines---experimental; not currently set in run_conversion.py
#
#############################################
def reload_IC( IC_file, format_type, which_dtype='f'):
    filetype=str(format_type)

    print("\n\n\n<----testing mn.reload_IC()---->")

    if filetype=='hdf5':
        hdf5_file=IC_file+'.hdf5'
        PartType=0
        masses=rw.read_block_single_file(hdf5_file,'Masses',PartType)[0][:]
        x=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,0]
        y=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,1]
        z=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,2]
    elif filetype=='gadget_binary':   
        f=open(IC_file+'.bin','r')
        ptype=0
        header=rw.load_gadget_binary_header(f)
        attribute_dictionary=rw.load_gadget_binary_particledat(f, header, ptype, skip_bh=0, which_dtype=which_dtype)

        positions=attribute_dictionary['Coordinates']
        masses=attribute_dictionary['Masses'] ##all of these are the same value, mp
        x=positions[:,0]
        y=positions[:,1]
        z=positions[:,2]

    else:
        masses,x,y,z,E=np.loadtxt(IC_file+'.txt',usecols=(0,1,2,3,4), unpack=True)

    r_recovered= np.sqrt(x**2.0 + y**2.0 + z**2.0)#*float(rmax)

    return r_recovered, masses



def binned_r_rho(r_array,mp,nbin):
    rmin=r_array.min()
    rmax=r_array.max()
    binsize=(rmax-rmin)/nbin
    r_set=np.arange(rmin,rmax,binsize)
    r_b=[]
    rho_b=[]
    for i in range(len(r_set)-1):
        r1=r_set[i]
        r2=r_set[i+1]
        region=np.where( (r1<=r_array) &(r2>r_array))  #size of this should be ~12N^2
        if len(r_array[region])==0:
            break
        r_b.append(r2)
        rho_b.append( len(r_array[region])*mp/(cf.volume(r2)-cf.volume(r1))  )
    return cf.to_array(r_b), cf.to_array(rho_b)


def bins_from_NR(NR_file_name, r_array, mp):
    N,r_set=np.loadtxt(NR_file_name,usecols=(0,1),unpack=True)
    r=[]
    rho=[]
    for i in range(len(r_set)-1):
        r1=r_set[i]
        r2=r_set[i+1]
        try:
            region=np.where( (r1<=r_array) &(r2>r_array))  #size of this should be ~12N^2
        except RuntimeWarning:
            print("Data types in generated vs recovered IC files do not match")
            print("Update 'which_dtype' value in config file")
            sys.exit()

        if len(r_array[region])==0:
            break
        r.append(r2)
        rho.append( len(r_array[region])*mp/(cf.volume(r2)-cf.volume(r1))  )
    return cf.to_array(r), cf.to_array(rho)

## end module mainlib
