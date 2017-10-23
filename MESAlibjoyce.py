#!/usr/bin/env python
import numpy as np
import codecs, re
import subprocess, os
import h5py
import pygadgetreader as pyg
import matplotlib.pyplot as plt


def plotter(xmaj,xmin,ymaj,ymin, xf, yf):
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

# def myfunc(a,b, *args, **kwargs):
#       c = kwargs.get('c', None)
#       d = kwargs.get('d', None)


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
		# not robust but whatever for now
		if "zone" in line and "num_zones" not in line:
			try:
				p=line.split()
			except:
			 	pass
	phys_dict={}		 	
	for i,v in enumerate(p):
		phys_dict[str(v)]=i#-1 	#wow fucking damn if this was the problem 
        #yes, it was. 
	return phys_dict


# fname='profile175.data'#'profile175.data' #140, 32
# MESA_file="{}".format(fname)
# print get_MESA_output_fields(MESA_file)



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


##########################################################################

def fit_MESA_density_profile(readfile):
    Rsolar=6.955*10.0**10.0 #cm
    try:
        logr=get_quantity(readfile,'logR')
        logrho=get_quantity(readfile,'logRho')
    except:
        print show_allowed_MESA_keywords(readfile)      

    radius=[]
    for i,p in enumerate(logr):
        try:
            radius.append(10.0**float(p))
        except:
            #logrho=list(logrho)
            #del logrho[i]
            radius.append(0.0)
    rho=[]
    for j,q in enumerate(logrho):
        try:
            rho.append(10.0**float(q))
        except:
            #logr=list(logr)
            #del logr[i]
            rho.append(0.0)
    radius = [Rsolar*r for r in radius]         
    return rho, radius



#######################################################################
#
# Gadget Input Functions
#
#######################################################################

def write_IC_binary(out_fname):
    print "pyg: ",pyg
    return




def make_IC_box_hdf5(out_fname):
    DIMS=2; # number of dimensions 
    N_1D=32; # 1D particle number (so total particle number is N_1D^DIMS)
    fname=out_fname#'box_3d_r32.hdf5'; # output filename 

    Lbox = 1.0 # box side length
    rho_desired = 1.0 # box average initial gas density
    P_desired = 1.0 # initial gas pressure
    vgrainrms = 1.0 # rms velocity of collisionless particles
    dust_to_gas_ratio = 0.01 # mass ratio of collisionless particles to gas
    gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
    
    # first we set up the gas properties (particle type 0)
    
    # make a regular 1D grid for particle locations (with N_1D elements and unit length)
    x0=np.arange(-0.5,0.5,1./N_1D); x0+=0.5*(0.5-x0[-1]);
    # now extend that to a full lattice in DIMS dimensions
    if(DIMS==3):
        xv_g, yv_g, zv_g = np.meshgrid(x0,x0,x0, sparse=False, indexing='xy')
    if(DIMS==2):
        xv_g, yv_g = np.meshgrid(x0,x0, sparse=False, indexing='xy'); zv_g = 0.0*xv_g
    if(DIMS==1):
        xv_g=x0; yv_g = 0.0*xv_g; zv_g = 0.0*xv_g; 

    # the gas particle number is the lattice size: this should be the gas particle number
    Ngas = xv_g.size

    # flatten the vectors (since our ICs should be in vector, not matrix format): just want a
    #  simple list of the x,y,z positions here. Here we multiply the desired box size in
    xv_g=xv_g.flatten()*Lbox; yv_g=yv_g.flatten()*Lbox; zv_g=zv_g.flatten()*Lbox; 

    # set the initial velocity in x/y/z directions (here zero)
    vx_g=0.*xv_g; vy_g=0.*xv_g; vz_g=0.*xv_g;

    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx_g=0.*xv_g; by_g=0.*xv_g; bz_g=0.*xv_g;

    # set the particle masses. Here we set it to be a list the same length, with all the same mass
    #   since their space-density is uniform this gives a uniform density, at the desired value
    mv_g=rho_desired/((1.*Ngas)/(Lbox*Lbox*Lbox)) + 0.*xv_g

    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    uv_g=P_desired/((gamma_eos-1.)*rho_desired)

    # set the gas IDs: here a simple integer list
    id_g=np.arange(1,Ngas+1)
    
    # now we set the properties of the collisionless particles: we will assign these to particle type '3', 
    #   but (barring special compile-time flags being set) GIZMO will treat all collisionless particle types
    #   the same. so the setup would be identical for any of the particle types 1,2,3,4,5

    # set the desired number of particles (here to about twice as many as the gas particles, because we feel like it)
    Ngrains = np.round(2. * (1.*N_1D)**DIMS)
    Ngrains=int(Ngrains) #mjoyce

    # set the x/y/z positions: again a simple list for each: here to random numbers from a uniform distribution
    print Ngrains, type(Ngrains)
    xv_d = (np.random.rand(Ngrains)-0.5)*Lbox

    yv_d = (np.random.rand(Ngrains)-0.5)*Lbox
    zv_d = (np.random.rand(Ngrains)-0.5)*Lbox

    # set the IDs: these must be unique, so we start from the maximum gas ID and go up
    id_d = np.arange(Ngas+1,Ngrains+Ngas+1)

    # set the velocities. again we will set to a random value, here a Gaussian-distributed one
    vx_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vy_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vz_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)

    # set the masses, again a list with all the same mass
    mv_d = dust_to_gas_ratio * (1.*Ngas)/(1.*Ngrains) * mv_g[0] + 0.*xv_d

    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
    file = h5py.File(fname,'w') 

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    npart = np.array([Ngas,0,0,Ngrains,0,0]) # we have gas and particles we will set for type 3 here, zero for all others

    # now we make the Header - the formatting here is peculiar, for historical (GADGET-compatibility) reasons
    h = file.create_group("Header");

    # here we set all the basic numbers that go into the header
    # (most of these will be written over anyways if it's an IC file; the only thing we actually *need* to be 'correct' is "npart")
    h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs 
                                         #  if we make a multi-part IC file. with this simple script, we aren't equipped to do that.
   
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since 
                                        # we set masses explicitly by-particle this should be zero. that is more flexible anyways, as it 
                                        # allows for physics which can change particle masses 
    

    ## all of the parameters below will be overwritten by whatever is set in the run-time parameterfile if
    ##   this file is read in as an IC file, so their values are irrelevant. they are only important if you treat this as a snapshot
    ##   for restarting. Which you shouldn't - it requires many more fields be set. But we still need to set some values for the code to read
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
    ## ok, that ends the block of 'useless' parameters
    
    # Now, the actual data!
    #   These blocks should all be written in the order of their particle type (0,1,2,3,4,5)
    #   If there are no particles of a given type, nothing is needed (no block at all)
    #   PartType0 is 'special' as gas. All other PartTypes take the same, more limited set of information in their ICs
    
    # start with particle type zero. first (assuming we have any gas particles) create the group 
    p = file.create_group("PartType0")

    # now combine the xyz positions into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=xv_g; q[:,1]=yv_g; q[:,2]=zv_g;

    # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q)

    # similarly, combine the xyz velocities into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=vx_g; q[:,1]=vy_g; q[:,2]=vz_g;

    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=q)

    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)

    # write particle masses to the Masses block
    p.create_dataset("Masses",data=mv_g)

    # write internal energies to the InternalEnergy block
    p.create_dataset("InternalEnergy",data=uv_g)

    # combine the xyz magnetic fields into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=bx_g; q[:,1]=by_g; q[:,2]=bz_g;

    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with 
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can 
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=q)

    # no PartType1 for this IC
    # no PartType2 for this IC

    # now assign the collisionless particles to PartType3. note that this block looks exactly like 
    #   what we had above for the gas. EXCEPT there are no "InternalEnergy" or "MagneticField" fields (for 
    #   obvious reasons). 
    p = file.create_group("PartType3")
    q=np.zeros((Ngrains,3)); q[:,0]=xv_d; q[:,1]=yv_d; q[:,2]=zv_d;
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((Ngrains,3)); q[:,0]=vx_d; q[:,1]=vy_d; q[:,2]=vz_d;
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_d)
    p.create_dataset("Masses",data=mv_d)

    # no PartType4 for this IC
    # no PartType5 for this IC

    # close the HDF5 file, which saves these outputs
    file.close()
    return file


##################################

# better form to put readwritehdf5 in here?
# include readsnap_m_edit.py?
