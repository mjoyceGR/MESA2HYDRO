#!/usr/bin/env python
#try:
import numpy as np
import codecs, re
import subprocess, os
import h5py
import pygadgetreader as pyg
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter 
import healpy as hp
import random as rand
#except:
#    print 'Missing module!\nThe following are required: '
#    print 'healpy\n'
#   exit(0)


def plotter(xmaj, xmin, ymaj, ymin,**kwargs):
    xf = kwargs.get('xf', '%1.1f')
    yf = kwargs.get('yf', '%1.1f')

    majorLocator_x  = MultipleLocator(xmaj)     # I want a major tick every "number"
    majorFormatter_x = FormatStrFormatter(xf)#('%1.1f')     # label these (the major ones) with a 1.2f format string
    minorLocator_x  = MultipleLocator(xmin)     # I want a minor tick every "number"

    majorLocator_y  = MultipleLocator(ymaj)     # now for the y axis...
    majorFormatter_y = FormatStrFormatter(yf)#('%1.1f')     # 
    minorLocator_y  = MultipleLocator(ymin)     #

    fig, ax = plt.subplots()

    ax.xaxis.set_major_locator(majorLocator_x)
    ax.xaxis.set_major_formatter(majorFormatter_x)
    ax.xaxis.set_minor_locator(minorLocator_x)

    ax.yaxis.set_major_locator(majorLocator_y)
    ax.yaxis.set_major_formatter(majorFormatter_y)
    ax.yaxis.set_minor_locator(minorLocator_y)

    return fig,ax

###
# ideally need dictionary or some hash table containing map of i for p[i] to physical quantity name 
### even better than this: read in all the column titles from the MESA output in one shot, shove them in a list,
# iterate over that list using enumerate or whatever, then BUILD a dictionary out of names AS LOADED FROM THE THING and assign an index to that 
# simply by iteration


def submit_MESA_job(mesa_work_path,qsub_script):
    mesa_work_path=str(mesa_work_path)
    os.chdir(mesa_work_path)
    #	subprocess.call('cd '+mesa_work_path, shell=True)
    print "entered"
    subprocess.call('pwd', shell=True)
    compile_MESA()
    run_MESA()
    return 

def compile_MESA():
	try:
		subprocess.call('./mk', shell=True)
	except:
		print 'MESA star compilation failed'
	return

def execute_MESA():
	subprocess.call('cd '+mesa_work_path, shell=True)
	try:
		subprocess.call('./rn', shell=True)
	except:
		print 'MESA star execution failed'	
	return


def update_MESA_inlist_value(MESA_inlist, field, value):
    inlist_dict=grab_fields(MESA_inlist)
    old_value=inlist_dict.get(field)
    new_value=str(value)
    oldstr=r'{}\s*=\s*{}'.format(field, old_value)
    newstr='{} = {}'.format(field, value)
    #print("oldstr: {}\nnewstr: {}".format(oldstr, newstr))
    f = codecs.open(MESA_inlist)
    contents = f.read()
    newcontents=re.sub(oldstr, newstr, contents) 
    f.close()
    outf=open(MESA_inlist,"w")
    print >>outf, newcontents
    outf.close()
    return 

def grab_fields(MESA_inlist):
    inf=open(MESA_inlist,'r')
    inlist_dict={}
    for line in inf:
        if line and "=" in line:
            try:
                p=line.split("=")
                try: 
                    p[1]=p[1].split('\n')[0]
                except:
                    pass
                if "!" in p[0]:
                    del p    
                inlist_dict[str(p[0]).strip()]=str(p[1]).strip()
            except:
                pass 
    return inlist_dict


#########################################################################
#
# MESA output functions
#
#########################################################################
def strip_MESA_header(in_filename, out_filename, *args, **kwargs):

    ### add hook to prevent double removes

	#[0]-- returns the file object
	#[1]-- returns the string/name of the reformatted text file
	n = int(kwargs.get('n', 5))
	num_delete=int(n) #the number of line to be read and deleted
	outfn=out_filename

	with open(in_filename) as f:
	    mylist = f.read().splitlines()
	newlist = mylist[:]
	#os.remove("bigFile.txt")

	thefile = open(outfn, 'w')
	del mylist[:num_delete]
	for item in mylist:
		thefile.write("%s\n" % item)
	return  thefile, outfn

def get_MESA_output_fields(filename):
    inf=open(filename,'r')
    for line in inf:
        #print "for line in inf"
        # not robust but whatever for now
        #try:
            #print "in try catch"
        if "zone" in line and "num_zones" not in line:
            #try:
            p=line.split()
            #except:
            #    strip_MESA_header(filename,filename)
            #    p=line.split()
        #except:
        else:
            #print "history"
            if 'model_number' in line:                    
                #try:
                p=line.split()
                #except:
                #    strip_MESA_header(filename,filename)
                #    p=line.split()
    phys_dict={}		 	
    for i,v in enumerate(p):
        phys_dict[str(v)]=i#-1 	#wow fucking damn if this was the problem 
        #yes, it was. 
    return phys_dict


def get_columns(filename,keyname_list):
    phys_dict=get_MESA_output_fields(filename)
    column_dict={}
    for i in range(len(keyname_list)):
        try:
            keyname=str(keyname_list[i])
            column_dict[ keyname ]=(get_key(filename, keyname)[0])
        except:
            print show_allowed_MESA_keywords(readfile)
    return column_dict


def get_key(filename,keyname):
    inf=open(filename,'r')
    phys_dict=get_MESA_output_fields(filename)
    indx=phys_dict.get(str(keyname))
    param=[]
    for line in inf:
      #avoid 'zone' instead of '#'
      if line and 'zone' not in line:
               #print line
          try: 
              p=line.split()
              param.append(p[indx])
          except:
              pass
    return param, type(param)

def get_quantity(readfile,keyname):
    keyname=str(keyname)
    keyname_list=get_MESA_output_fields(readfile).keys()
    column_dict=get_columns(readfile,keyname_list)
    quantity=np.array(column_dict.get(keyname)).astype(float)
    return quantity


def show_allowed_MESA_keywords(readfile):
    fstr=""
    for i in get_MESA_output_fields(readfile).keys():
        fstr=fstr+str(i) +'\n'
    return fstr#get_MESA_output_fields(readfile).keys()

###########################################################################
#
# healpix
#
############################################################################
def get_coords(N,r_mid, rmax):#,file_index):
    r_mid=float(r_mid); rmax=float(rmax)
    ##Need N to be (4) 8 or 16, and must change mp per shell to maintain that, probably
    NSIDE = N #closest power of 2 to N
    NSIDE=int(NSIDE)

    theta=random_theta()
    ipix_array=np.arange(hp.nside2npix(NSIDE)) #this is just a straight up array of particle IDs
    x=[]
    y=[]
    z=[]
    for i in range(len(ipix_array)):
        ipix=ipix_array[i]
        coord=hp.pixelfunc.pix2vec(NSIDE, ipix, nest=True)
        #print >> outf, coord[0], coord[1], coord[2]
        x.append(coord[0])
        y.append(coord[1])
        z.append(coord[2])

    # no fam this rotate shell thing is some serious math
    xd, yd, zd = rotate_shell(x,y,z,theta,"about_z")
    xe, ye, ze = rotate_shell(xd,yd,zd,theta,"about_y")
    xf, yf, zf = rotate_shell(xe,ye,ze,theta,"about_x")

    xf = np.array(xf).astype(float)
    yf = np.array(yf).astype(float)
    zf = np.array(zf).astype(float)

    xf = to_physical(xf, r_mid/rmax)
    yf = to_physical(yf, r_mid/rmax)
    zf = to_physical(zf, r_mid/rmax)

    xf = np.array(xf).astype(float)
    yf = np.array(yf).astype(float)
    zf = np.array(zf).astype(float)

    #print "x: ",xf, "\n\ny: ", yf, "\n\nz: ",zf
    return xf.flatten(), yf.flatten(), zf.flatten()


def to_physical(xq, r_mid):
    xf = [ (r_mid*x_i) for x_i in xq]
    return xf

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

def random_theta():
    theta=rand.random()*2.0*np.pi #.random gives random float between 0 and 1
    return theta


# x_array=np.array([3,3,3,3])
# y_array = 2.0+0.0*x_array
# z_array = 5.0+0.0*x_array

# x,y,z=rotate_shell(x_array,y_array,z_array, random_theta())

# print "\none column? x: ", x
# print "\nx[0]:", x[0]
# print "\nwhat? x[0][0]:", x[0][0]

# for i in range(len(x)):
#     print "array[i] etc: ", np.sqrt(x_array[i]**2.0+y_array[i]**2.0+z_array[i]**2.0)
#     print "x[i] etc", np.sqrt(x[i]**2.0 + y[i]**2.0 + z[i]**2.0)


def to_rad(theta):
    theta=theta*np.pi/180.0
    return theta




#######################################################################
#
# Gadget Input Functions
#
#######################################################################

def write_IC_binary(out_fname):
    print "pyg: ",pyg
    return


#def make_IC_hdf5(out_fname, mp, coord_file, **kwargs):
#x,y,z=np.loadtxt(coord_file, usecols=(0,1,2), unpack=True)
def make_IC_hdf5(out_fname, mp, x, y, z, rho,**kwargs):    
    ## from kwargs choose either hdf5 or binary, when I feel like doing this
    fname=out_fname
    Lbox = 1.0                  # box side length
    rho_desired = 1.0           # box average initial gas density
    P_desired = 1.0             # initial gas pressure
    vgrainrms=0.0               # mjoyce- for me, all particles stationary bc star
    dust_to_gas_ratio = 0.01    # mass ratio of collisionless particles to gas
    gamma_eos = 5./3.           # polytropic index of ideal equation of state the run will assume
    
    Ngas = len(x)#xv_g.size

    x=x*Lbox                    #positions
    y=y*Lbox 
    z=z*Lbox; 

    rho=rho

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
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs

    p = file.create_group("PartType0")

    q=np.zeros((Ngas,3)); q[:,0]=x; q[:,1]=y; q[:,2]=z;
    p.create_dataset("Coordinates",data=q)
    
    q=np.zeros((Ngas,3)); q[:,0]=vx; q[:,1]=vy; q[:,2]=vz;
    p.create_dataset("Velocities",data=q)
    
    p.create_dataset("ParticleIDs",data=id_g)
    p.create_dataset("Masses",data=mv_g)
    uv_g = U + 0.*xv_d#
    p.create_dataset("InternalEnergy",data=uv_g)
    q=np.zeros((Ngas,3)); q[:,0]=bx; q[:,1]=by; q[:,2]=bz;
    
    p.create_dataset("MagneticField",data=q)

    #q=
    p.create_dataset("Density",data=rho)

    # no PartType1 for this IC
    # no PartType2 for this IC
    # no PartType4 for this IC
    # no PartType5 for this IC
    file.close()
    return file
