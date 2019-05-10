#!/usr/bin/env python
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import sys 
import converge_funcs as cf
import MESAhandling as MJ
import io_lib as rw
import hdf5lib as hdf5lib
import time
########################################################
#
# this module contains the main subroutines for run.py 
# and subordinate minor subroutines
#
#########################################################

########### FOR DEBUGGING ONLY
use_normalized=False


print "\n\nusing git version\n\n"

###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################



def MESA_r(MESA_file, masscut):
    ## returns R in cm
    fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
    fit_region_R=cf.unlog(fit_region_R)*R_to_solar
    return fit_region_R


def MESA_rho(MESA_file,masscut):
    fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
    fit_region_rho=cf.unlog(fit_region_rho)
    return fit_region_rho

def MESA_m(MESA_file,masscut):
    ### RETURNS M ARRAY IN cgs UNITS!!!!
    fit_region_M=cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut,strip=False)*M_to_solar
    return fit_region_M

def MESA_internalE(MESA_file,masscut):
    fit_region_E=cf.get_MESA_profile_edge(MESA_file, quantity='logE', masscut=masscut,strip=False)
    fit_region_E=cf.unlog(fit_region_E)
    return fit_region_E



################################################### 

def central_mass(MESA_file, masscut):
    ################################################################
    #
    # last edited 4/26/19 by Mjoyce
    #
    #################################################################

    masses = MJ.get_quantity(MESA_file,'mass').astype(np.float)*M_to_solar ## WARNING
    radii = 10.0**(MJ.get_quantity(MESA_file, 'logR').astype(np.float))

    #Mentry=masses.max()#masses[0]
    Mtot=masses.max()#sum(masses)#masses[0]
    #print "Mentry at (loc 1) = ", Mentry
    print "Mtot at (loc 1) = ", Mtot

    # import matplotlib.pyplot as plt
    # plt.plot(radii, masses, "g-", label="mass vs radius")
    # plt.xlabel("radius Rsolar")
    # plt.ylabel("mass Msolar")
    # plt.show()
    # plt.close()
    # sys.exit()

    #print "masscut: ", masscut

    central_M = masscut*Mtot
    #print "central_M at (loc 2) = ", central_M

    # rest = 1.0e-7*(221185)#(473088) 
    # other_tot = rest + central_M
    # print "rest of mass at (loc (3) = ", rest
    # print "total mass at (loc 3) = ", other_tot

    # print "actual mass percentage represented in atm vs prescribed rate:", (1.0-rest/Mtot), "   ", masscut

    #sys.exit()


    return central_M#.astype(float)
################################################### mjoyce 11/2/2018




##########################################################
#
# Make an NR file from MESA data
#
###########################################################
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
    return 


def make_NR_file(MESA_file,masscut,N,mp, RKstep,NR_file, *args, **kwargs):
    ## NR_file must be FILE OBJECT, not STRING
    # mp must be passed in units of solar masses, e.g. 1e-7 

    start_time=time.time()
    fit_region_R=MESA_r(MESA_file, masscut)
    fit_region_rho=MESA_rho(MESA_file, masscut)

    fit_region_E=MESA_internalE(MESA_file,masscut)

    #mp=mp*M_to_solar
    rl=fit_region_R.min()
    rmax=fit_region_R.max()

    outf=NR_file
    print >> outf, '## fname',MESA_file ,' masscut',masscut,'   N', N, '  mp (Ms)', mp/M_to_solar,\
       '  mp (g)', mp,'  RKstep',('%1.3e'% RKstep)
    print >> outf, '#N    (ru+rl)/2 (cm)    M(g) contained in shell ru-rl   internal energy E (unlog)...ergs?'

    ru=rl
    while ru <= rmax:
        try:
            ru, Mshell=cf.get_placement_radii(rl, ru, RKstep, N, mp, MESA_file,masscut, suppress_output=False, load_unlogged=False)
            print >> outf, N, ru, Mshell, cf.get_logE(ru,MESA_file,masscut)
            rl=ru
        except TypeError:
            print 'reached', ('%1.5f'% (ru*100.0/rmax)), r'% of outer radius' 
            break

    # ################# mjoyce 11/2/18       
    # print >> outf, 

    print 'runtime: ', time.time()-start_time
    print >> outf, '#\n#\n# runtime: ', time.time()-start_time, " seconds"

    return



##########################################################
#
# Make the initial conditions file from an NR file
#
###########################################################
def get_IC(MESA_file, masscut, NR_file_name,output_filename,mp,which_dtype='f', *args, **kwargs): #temp remove rmax
    filetype=str(kwargs.get('format_type','binary'))
    #print 'WARNING!! sending physical radius!!!!\nmp IS multipled by Msolar'
    #print "\n\nWARNING! mp not multiplied by M_solar!! \n\n"
    #mp=mp*M_to_solar 

    N,rmid,E=np.loadtxt(NR_file_name, usecols=(0,1,3), unpack=True) 

    super_x=[]
    super_y=[]
    super_z=[]

    super_E=[]

    for i in range(len(N)):
        NSIDE= N[i]
        r_mid= rmid[i]
        E_val=E[i]
        # print 'WARNING!! sending physical radius!!!!\nmp IS multipled by Msolar'

        radius=float(rmid[i])

        if use_normalized:
            print "\n\nWARNING! NORMALIZED >>> radius <<< COORDINATES NECESSARY FOR BINARY FORMAT!"
            radius = radius/R_to_solar
        else:
            radius=radius
        #print "normalized radius: ", radius

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
        #same size as x, but same value of internal E for points at radius r
        E_array=x*0.+E_val

        super_x=np.concatenate((super_x,x),axis=0)
        super_y=np.concatenate((super_y,y),axis=0)
        super_z=np.concatenate((super_z,z),axis=0)
        #r_super=np.sqrt(super_x**2.0+ super_y**2.0 + super_z**2.0)
        super_E=np.concatenate((super_E,E_array),axis=0)

    ############################################################    
    #
    ### add in the central mass point at coordiante 0,0,0
    #
    ############################################################
    ##
    ## maybe don't need to do this here, because we can do this in io_lib directly?

    zero_space=np.array([0.0])

    super_x=np.concatenate((super_x,zero_space),axis=0)
    super_y=np.concatenate((super_y,zero_space),axis=0)
    super_z=np.concatenate((super_z,zero_space),axis=0)
    super_E=np.concatenate((super_E,zero_space),axis=0) ## don't know what to do about this

    central_point_mass=central_mass(MESA_file, masscut)
    #central_point_mass/M_to_solar

    ############################################################
    super_x=cf.to_array(super_x)
    super_y=cf.to_array(super_y)
    super_z=cf.to_array(super_z)
    super_E=cf.to_array(super_E)


    if use_normalized:
        print "\n\nWARNING! NORMALIZED >>> mass <<< COORDINATES NECESSARY FOR BINARY FORMAT!"
        mp = mp/M_to_solar
        central_point_mass= central_point_mass/M_to_solar
    else:
        mp = mp
        central_point_mass= central_point_mass


    if filetype=='hdf5':
        var=rw.make_IC_hdf5(str(output_filename)+ '.hdf5',\
        mp, central_point_mass,\
        super_x, super_y, super_z,super_E) #, userho=False
        #svar=rw.make_IC_hdf5_old_way(hdf5file, mp, super_x, super_y, super_z,super_E)



    elif filetype=='gadget_binary':
        #print "mp, ex z before norm: ", mp, super_z[5]
        # print "\n\nWARNING! NORMALIZED >>> mass <<< COORDINATES NECESSARY FOR BINARY FORMAT!"
        # #radius = radius/R_to_solar
        # mp = mp/M_to_solar
        # central_point_mass= central_point_mass/M_to_solar

        super_x=super_x #/(R_to_solar)**(0.5)
        super_y=super_y #/(R_to_solar)**(0.5)
        super_z=super_z #/(R_to_solar)**(0.5)
       # print "mp, ex z after norm: ",mp, super_z[5]
       #print super_z[5]

        var=rw.make_IC_binary(str(output_filename)+ '.bin',\
        mp, central_point_mass,\
        super_x, super_y, super_z, super_E, which_dtype=which_dtype)  #central mass not handled (??? is this true)

    

    else:
        var=rw.make_IC_text(str(output_filename)+ '.txt',\
        (super_x*0. + mp), central_point_mass,\
        super_x, super_y, super_z, super_E, which_dtype=which_dtype)

    print var, type(var)
    return



def estimate_stepsize(MESA_file, masscut, Nshells):
    fit_region_R=MESA_r(MESA_file, masscut)
    fit_region_rho=MESA_rho(MESA_file, masscut)
    rl=fit_region_R.min()
    rmax=fit_region_R.max()

    return (float(rmax)-float(rl))/float(Nshells)



def reload_IC( IC_file, format_type, which_dtype='f'): #rmax #NR_file
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
            print "Data types in generated vs recovered IC files do not match"
            print "Update 'which_dtype' value in config file"
            sys.exit()

        if len(r_array[region])==0:
            break
        r.append(r2)
        rho.append( len(r_array[region])*mp/(cf.volume(r2)-cf.volume(r1))  )
    return cf.to_array(r), cf.to_array(rho)






def quick_plot(MESA_file, masscut, r_reloaded,rho_reloaded,IC_format_type,png_tag='latest'):
    import matplotlib.pyplot as plt
    fit_region_R=MESA_r(MESA_file, masscut)
    fit_region_rho=MESA_rho(MESA_file, masscut)

    # plt.plot(r_reloaded, rho_reloaded,'r.', markersize=6, label='GADGET data')
    # plt.plot(fit_region_R, fit_region_rho, "b.", markersize=4, label='MESA data') #cf.to_log()
    # plt.xlabel("R")
    # plt.ylabel("test density")
    # plt.legend(loc=1)
    # #if IC_format_type=='hdf5':
    # #    plt.savefig('lin_'+png_tag+'_hdf5.png')
    # #else:
    # plt.savefig('lin_'+png_tag+'.png')
    # plt.close()

    plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=4, label='MESA data') #cf.to_log()
    plt.plot(r_reloaded, cf.to_log(rho_reloaded),'r.', markersize=6, label='GADGET data')
    plt.xlabel("R")
    plt.ylabel("log(test density)")
    plt.legend(loc=1)
    #if IC_format_type=='hdf5':
    #    plt.savefig('log_'+png_tag+'_hdf5.png')
    #else:
    plt.savefig('log_'+png_tag+'.png')
    plt.close()
