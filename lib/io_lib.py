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
def make_IC_hdf5(out_fname, mp, central_point_mass,\
         x, y, z,E, **kwargs):

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
    masses[-1]=central_point_mass

    ###########################################################################
    hsml=0.*x + (-1)
    print "5/21/19"
    print "hsml--smoothing length (loc 1): ", hsml
    ###########################################################################


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


    particles.create_dataset("ParticleIDs",data=IDs)
    particles.create_dataset("Masses",data=masses)#,dtype=np.dtype('d'))
    particles.create_dataset("Coordinates",data=pos)#,dtype=np.dtype('d'))
    particles.create_dataset("Velocities",data=vel)#,dtype=np.dtype('d'))
    particles.create_dataset("InternalEnergy",data=internalE)#,dtype=np.dtype('d'))
    particles.create_dataset("HSML", data=hsml)#,dtype=np.dtype('d'))
    particles.create_dataset("Density",data=dens)#,dtype=np.dtype('d'))

    file.close()
    return file

###############################################################
#
# binary writing routine
#
###############################################################




def make_IC_binary(fname, mp, central_point_mass,\
                   x, y, z,\
                   local_MESA_rho, local_MESA_P, local_MESA_E,\
                   which_dtype='f',**kwargs):
                   # E removed!!!

    #central_mass=float(kwargs.get('central_mass', 10e6)) #<-----WARNING!! not properly handled!!!
    print "mp in make_IC_binary is", mp
    import pyIC as pygadgetic

    gas_particles=len(x)-1#+1  ## mjoyce 11/2/2018  for central mass
    total_number_of_particles=len(x)

    #gas particles = npart[0] = gas_particles
    npart=[gas_particles,1,0,0,0,0] # central mass should be Type 1 according to The Phil

    my_header=pygadgetic.Header()


    my_body=pygadgetic.Body(npart)
    my_header.NumPart_ThisFile = np.array(npart)
    my_header.NumPart_Total = np.array(npart)

    my_body.pos[:,0]=x
    my_body.pos[:,1]=y
    my_body.pos[:,2]=z

    my_body.vel[:,0]=0.*total_number_of_particles
    my_body.vel[:,1]=0.*total_number_of_particles
    my_body.vel[:,2]=0.*total_number_of_particles


    ##fill the body with minimal information
    my_body.id[:]=np.arange(0,total_number_of_particles) #id_g

    #########################
    #
    # WARNING! confused here!
    #
    ##########################
    my_body.mass[:]=0.*total_number_of_particles + mp #particle_masses 

    ############################
    #
    # GADGET "knows" that the last index'd particle is meant to be the Type 1 particle
    #
    ############################
    print "value of central_point_mass: ", central_point_mass
    my_body.mass[-1]=central_point_mass 


    #########################
    #
    # WARNING! confused here!
    #
    ##########################
    h=1.2*(mp/local_MESA_rho)**(1.0/3.0)
    h=2.0*h ### GADGET defines hsml as the kernel radius, so use 2H here
    
    # Q: this should... include the central particle? yes no? 5/21/19
    # A: maybe doesn't matter, but gas_particles/total_num_particles here
    #    needs to MATCH the definition in the pyIC binary writer

    my_body.hsml[:]=h[0:gas_particles]  
    #my_body.hsml[:]=0.*gas_particles + (-1)  


    ### give hydro MESA's density estimate directly, to ease starting guess 
    my_body.rho[:]=local_MESA_rho[0:gas_particles]

    ## give hydro some energy per unit mass from 
    ## an ideal gas equation of state gamma=5/3 (and then hydro code has to use that EOS too)
    gamma_eos    = 5.0/3.0     
    internal_E   = local_MESA_P/((gamma_eos-1.0)*local_MESA_rho)
    my_body.u[:] = internal_E[0:gas_particles]

    # sanity check that MESA and ideal EOS look the same
    #my_body.u[:]=local_MESA_E[0:gas_particles]

    
    pygadgetic.check_header(my_header)
    pygadgetic.dump_ic(my_header,my_body,fname, which_dtype=which_dtype)
    return fname


###############################################################
#
# Safe Dumb Text, no fancy formatting
#
###############################################################
def make_IC_text(\
    fname, mp, central_point_mass,\
    x, y, z, E,\
    which_dtype='f',**kwargs):
    ################################################################
    #
    # last edited 4/26/19 by Mjoyce
    #
    #################################################################

    print "central mass: ", central_point_mass, "  mp: ", mp

    format_str="%.12f"

    outf=open(fname, "w")
    lines = [( (format_str%mp[i])\
          +' '+(format_str%x[i])\
          +' '+(format_str%y[i])\
          +' '+(format_str%z[i])\
          +' '+(format_str%E[i]))\
          for i in range(len(x)-1)]

    print >> outf, "\n".join(lines)

    print >> outf, (format_str%central_point_mass)\
          +' '+(format_str%0)\
          +' '+(format_str%0)\
          +' '+(format_str%0)\
          +' '+(format_str%0)

    outf.close()
    return fname


