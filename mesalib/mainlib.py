#!/usr/bin/env python
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import sys 
import converge_funcs as cf
import io_lib as rw
import hdf5lib as hdf5lib

import time

print "\n\nusing git version\n\n"

###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################



def MESA_r(MESA_file, masscut):
    fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
    fit_region_R=cf.unlog(fit_region_R)*R_to_solar
    return fit_region_R


def MESA_rho(MESA_file,masscut):
    fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
    fit_region_rho=cf.unlog(fit_region_rho)
    return fit_region_rho



##########################################################
#
# Make an NR file from MESA data
#
###########################################################
def make_NR_file(MESA_file,masscut,N,mp, RKstep,NR_file, *args, **kwargs):
    ## NR_file must be FILE OBJECT, not STRING
    # mp must be passed in units of solar masses, e.g. 1e-7 


    check_MESA=kwargs.get('check_MESA',False)
    start_time=time.time()

    #RKstep=RKstepsize
    # outer_step=outer_step

    fit_region_R=MESA_r(MESA_file, masscut)
    fit_region_rho=MESA_rho(MESA_file, masscut)

    fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)*M_to_solar


    if check_MESA:
        import matplotlib.pyplot as plt
        plt.plot(fit_region_R, fit_region_rho,'c.')
        plt.title(MESA_file)
        plt.show()
        plt.close()

    mp=mp*M_to_solar
    rl=fit_region_R.min()
    rmax=fit_region_R.max()

    #if make_NR_file:
    #saveNR=str(NR_file_name)
    #outf=open(saveNR,"w")
    outf=NR_file
    print >> outf, '## fname',MESA_file ,' masscut',masscut,'   N', N, '  mp', mp/M_to_solar,\
       'mp_solar', mp,'  RKstep',('%1.3e'% RKstep)#, '   outer_step', ('%1.3e'%outer_step)
    print >> outf, '#N    (ru+rl)/2 (cm)  M(g) contained in shell ru-rl'

    ru=rl
    while ru <= rmax:
        try:
            ru, Mshell=cf.get_placement_radii(rl, ru, RKstep, N, mp, MESA_file,masscut, suppress_output=False, load_unlogged=False)
            print >> outf, N, ru, Mshell
            rl=ru
        except TypeError:
            print 'reached', ('%1.5f'% (ru*100.0/rmax)), r'% of outer radius' 
            break

    print 'runtime: ', time.time()-start_time
    print >> outf, '#\n#\n# runtime: ', time.time()-start_time, " seconds"
    #outf.close()

    return



##########################################################
#
# Make the initial conditions file from an NR file
#
###########################################################
def get_IC(NR_file_name,output_filename,mp, *args, **kwargs): #temp remove rmax
    filetype=str(kwargs.get('format_type','binary'))
    
    print 'WARNING!! sending physical radius!!!!\nmp IS multipled by Msolar'

    #print "\n\nWARNING! mp not multiplied by M_solar!! \n\n"
    mp=mp*M_to_solar 

    #print "\n\n\nmp:", mp, "\n\n\n"

    #tag=str(output_filename)
    hdf5file=str(output_filename)+ '.hdf5'
    binaryfile=str(output_filename)+ '.bin'
    #saveNR=str(NR_file_name)

    N,rmid=np.loadtxt(NR_file_name, usecols=(0,1), unpack=True) 
    #if make_IC_file:
    super_x=[]
    super_y=[]
    super_z=[]
    for i in range(len(N)):
        NSIDE= N[i]
        r_mid= rmid[i]

        # print 'WARNING!! sending physical radius!!!!\nmp IS multipled by Msolar'
        radius=float(rmid[i])
        #radius=float(rmid[i])/float(rmax) ##normalized coordinates
        ###############################################################
        #
        # Arbitrarily rotate shells
        #
        ###############################################################
        x,y,z=cf.get_coords(NSIDE)
        x=x*radius
        y=y*radius
        z=z*radius
        super_x=np.concatenate((super_x,x),axis=0)
        super_y=np.concatenate((super_y,y),axis=0)
        super_z=np.concatenate((super_z,z),axis=0)
        r_super=np.sqrt(super_x**2.0+ super_y**2.0 + super_z**2.0)
    super_x=cf.to_array(super_x)
    super_y=cf.to_array(super_y)
    super_z=cf.to_array(super_z)
    if filetype=='hdf5':
        var=rw.make_IC_hdf5(hdf5file, mp, super_x, super_y, super_z) #, userho=False
    else:
        var=rw.make_IC_binary(binaryfile, mp, super_x, super_y, super_z, central_mass=1) 
    print var, type(var)
    return


def estimate_stepsize(MESA_file, masscut, Nshells):

    # load density, mass, radius from MESA and convert to cgs
    # fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
    # fit_region_R=cf.unlog(fit_region_R)*R_to_solar
    # fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)*M_to_solar
    # fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
    # fit_region_rho=cf.unlog(fit_region_rho)
    fit_region_R=MESA_r(MESA_file, masscut)
    fit_region_rho=MESA_rho(MESA_file, masscut)

    rl=fit_region_R.min()
    rmax=fit_region_R.max()


    return (float(rmax)-float(rl))/float(Nshells)



def reload_IC( IC_file, format_type): #rmax #NR_file
    #
    filetype=str(format_type)

    print "\n\n\n<----testing reload---->"
    # N,r_set,m_cont=np.loadtxt(NR_file, usecols=(0,1,2), unpack=True)
    if filetype=='hdf5':
        hdf5_file=IC_file+'.hdf5'
        PartType=0
        masses=rw.read_block_single_file(hdf5_file,'Masses',PartType)[0][:]
        x=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,0]
        y=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,1]
        z=rw.read_block_single_file(hdf5_file,"Coordinates",PartType)[0][:][:,2]
    else:   
        f=open(IC_file+'.bin','r')
        ptype=0
        header=rw.load_gadget_binary_header(f)
        attribute_dictionary=rw.load_gadget_binary_particledat(f, header, ptype, skip_bh=0)

        positions=attribute_dictionary['Coordinates']
        masses=attribute_dictionary['Masses'] ##all of these are the same value, mp
        x=positions[:,0]
        y=positions[:,1]
        z=positions[:,2]

    r_recovered= np.sqrt(x**2.0 + y**2.0 + z**2.0)#*float(rmax)
    # p_mass=masses[5]

    # r1=r_recovered.min() 
    # r_temp=[]
    # rho_temp=[]
    # for i in range(len(r_set)-1):
    #     r1=r_set[i]
    #     r2=r_set[i+1]
    #     region=np.where( (r1<=r_recovered) &(r2>r_recovered))

    #     if len(r_recovered[region])==0:
    #         break
    #     #print "length of r_recovered*mp", len(r_recovered[region])*p_mass
    #     #print "volume of shell r2-r1: ", cf.volume(r2)-cf.volume(r1)
    #     #print "r_rec[ r1 to r2 ]*mp / vol(r2-r1): ", len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))
    #     r_temp.append(r2)
    #     rho_temp.append( len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))  )


    return r_recovered, masses
    #cf.to_array(r_temp), cf.to_array(rho_temp)