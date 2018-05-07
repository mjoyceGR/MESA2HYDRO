#!/usr/bin/env python
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import sys 
import converge_funcs as cf

######
import hdf5lib as hdf5lib

from yanked import *

###############################################################
#
# hdf5 writing routine
#
###############################################################
def make_IC_hdf5(out_fname, mp, x, y, z,E, **kwargs):

    print "\n\n\nusing modified hdf5 writer 5/3/18\n"

    file = h5py.File(out_fname,'w') 
    Ngas = len(x)
    npart = np.array([len(x),0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others

    h = file.create_group("Header");
    h.attrs['NumPart_ThisFile'] = npart; 
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.array([mp,0,0,0,0,0]) #np.zeros(6); 
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = 1.0; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs

    particles = file.create_group("PartType0")

    rho_desired = 1.0           # box average initial gas density
    P_desired = 1.0             # initial gas pressure
    gamma_eos = 5./3.           # polytropic index of ideal equation of state the run will assume

    IDs=np.arange(1,Ngas+1)
    masses=mp + 0.*x 
    hsml=0.*x + (-1)
    #U=P_desired/((gamma_eos-1.)*rho_desired) 
    #internalE = U + 0.*x
    internalE=E
    dens=0.*x 

    pos=np.zeros((Ngas,3))
    pos[:,0]=x
    pos[:,1]=y
    pos[:,2]=z

    vel=np.zeros((Ngas,3))
    vel[:,0]=0.*x
    vel[:,1]=0.*x
    vel[:,2]=0.*x

    print "masses being used in hdf5: ", masses
    print 'positions being used in hdf5: ', pos
   #internal energy 
    #particle IDs

    particles.create_dataset("ParticleIDs",data=IDs)
    particles.create_dataset("Masses",data=masses)#,dtype=np.dtype('d'))
    particles.create_dataset("Coordinates",data=pos)#,dtype=np.dtype('d'))
    particles.create_dataset("Velocities",data=vel)#,dtype=np.dtype('d'))
    particles.create_dataset("InternalEnergy",data=internalE)#,dtype=np.dtype('d'))
    particles.create_dataset("HSML", data=hsml)#,dtype=np.dtype('d'))
    particles.create_dataset("Density",data=dens)#,dtype=np.dtype('d'))

    file.close()
    return file


def make_IC_hdf5_old_way(out_fname, mp, x, y, z,E, **kwargs):    #rho out
    print "\n\n\n\nusing sketchy version of hdf5 writer\n\n\n\n"

    ## from kwargs choose either hdf5 or binary, when I feel like doing this
    userho=kwargs.get("userho",False)

    fname=out_fname
    Lbox = 1.0                  # box side length
    rho_desired = 1.0           # box average initial gas density
    P_desired = 1.0             # initial gas pressure
    vgrainrms=0.0               # mjoyce- for me, all particles stationary bc star
    dust_to_gas_ratio = 0.01    # mass ratio of collisionless particles to gas
    gamma_eos = 5./3.           # polytropic index of ideal equation of state the run will assume
    #temp = 

    Ngas = len(x)#xv_g.size

    x=x*Lbox                    #positions
    y=y*Lbox 
    z=z*Lbox; 


    vx=0.*x                     #3D velocities
    vy=0.*y                     
    vz=0.*z

    bx=0.*x                     #magnetic field components
    by=0.*y
    bz=0.*z

    mv_g=mp + 0.*x              #particle mass
    U=P_desired/((gamma_eos-1.)*rho_desired)    #internal energy 
    id_g=np.arange(1,Ngas+1)    #particle IDs

    xv_d=x
    yv_d=y
    zy_d=z

    Ngrains=int(len(xv_d)) #length of arrays loaded from healpix
    print Ngrains, type(Ngrains)
    ###only pick one of Ngrains / Ngas because we only want one particle type

    # set the IDs: these must be unique, so we start from the maximum gas ID and go up
    id_d = np.arange(Ngas+1,Ngrains+Ngas+1)

    vx_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vy_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vz_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)

    ############# mjoyce ########################################
    # if we had other particle types besides gas, set their masses separately 
    # mv_d = mp + 0.*xv_d

    file = h5py.File(fname,'w') 
    npart = np.array([Ngas,0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others
    h = file.create_group("Header");
    h.attrs['NumPart_ThisFile'] = npart; 
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); 
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = 1.0; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 1; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs

    p = file.create_group("PartType0")

    q=np.zeros((Ngas,3)); q[:,0]=x; q[:,1]=y; q[:,2]=z;
    p.create_dataset("Coordinates",data=q)
    
    q=np.zeros((Ngas,3)); q[:,0]=vx; q[:,1]=vy; q[:,2]=vz;
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_g)
    p.create_dataset("Masses",data=mv_g)
    #rho_desired = 1.0 
    uv_g = U + 0.*xv_d#
    p.create_dataset("InternalEnergy",data=E)

    q=np.zeros((Ngas,3)); q[:,0]=bx; q[:,1]=by; q[:,2]=bz;
    p.create_dataset("MagneticField",data=q)
    ### force density to be zero?
    rho = rho_desired + 0.*xv_d
    p.create_dataset("Density",data=rho)

    # no PartType1 for this IC
    # no PartType2 for this IC
    # no PartType4 for this IC
    # no PartType5 for this IC
    file.close()
    return file




###############################################################
#
# binary writing routine
#
###############################################################
def make_IC_binary(fname, mp, x, y, z, E, which_dtype='f',**kwargs):
    central_mass=float(kwargs.get('central_mass', 10e6)) #<-----WARNING!! not properly handled!!!
    print "mp in make_IC_binary is", mp
 
    import pyIC as pygadgetic
    total_number_of_particles=len(x)
    npart=[total_number_of_particles,0,0,0,0,0]

    my_header=pygadgetic.Header()
    my_body=pygadgetic.Body(npart)
    my_header.NumPart_ThisFile = np.array(npart)
    my_header.NumPart_Total = np.array(npart)

    ##fill the body with minimal information
    my_body.id[:]=np.arange(0,total_number_of_particles) #id_g
    my_body.mass[:]=0.*x + mp #particle_masses
    my_body.hsml[:]=0.*x + (-1)  #shazrene says -1 but that's not working at all

    my_body.pos[:,0]=x
    my_body.pos[:,1]=y
    my_body.pos[:,2]=z

    my_body.vel[:,0]=0.*x
    my_body.vel[:,1]=0.*x
    my_body.vel[:,2]=0.*x

    rho_desired = 1.0
    P_desired = 1.0 
    gamma_eos = 5./3.
  #  U=P_desired/((gamma_eos-1.)*rho_desired) # internal energy 
  #U + 0.*x#internal_energy ### <--------- load this from MESA directly
    ## trying to load this from MESA!
    my_body.u[:]=E   ##this is in ergs maybe??? unclear??


    pygadgetic.dump_ic(my_header,my_body,fname, which_dtype=which_dtype)
    return fname
