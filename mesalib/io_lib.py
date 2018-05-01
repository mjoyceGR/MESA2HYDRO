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

###################################################
#
# binary routines
#
####################################################
def load_gadget_binary_header(f):
    import array
    f.read(4)
    Npart = array.array('I')
    Npart.fromfile(f, 6)
    Massarr = array.array('d')
    Massarr.fromfile(f, 6)
    a = array.array('d')
    a.fromfile(f, 1)
    a = a[0]
    z = array.array('d')
    z.fromfile(f, 1)
    z = float(z[0])
    FlagSfr = array.array('i')
    FlagSfr.fromfile(f, 1)
    FlagFeedback = array.array('i')
    FlagFeedback.fromfile(f, 1)
    Nall = array.array('i')
    Nall.fromfile(f, 6)
    FlagCooling = array.array('i')
    FlagCooling.fromfile(f, 1)
    NumFiles = array.array('i')
    NumFiles.fromfile(f, 1)
    BoxSize = array.array('d')
    BoxSize.fromfile(f, 1)
    Omega0 = array.array('d')
    Omega0.fromfile(f, 1)
    OmegaLambda = array.array('d')
    OmegaLambda.fromfile(f, 1)
    h = array.array('d')
    h.fromfile(f, 1)
    h = float(h[0])
    FlagAge = array.array('i')
    FlagAge.fromfile(f, 1)
    FlagMetals = array.array('i')
    FlagMetals.fromfile(f, 1)
    NallHW = array.array('i')
    NallHW.fromfile(f, 6)
    flag_entr_ics = array.array('i')
    flag_entr_ics.fromfile(f, 1)
    f.seek(4+256+4+4)
    return {'NumPart_ThisFile':Npart, 'MassTable':Massarr, 'Time':a, 'Redshift':z, \
    'Flag_Sfr':FlagSfr[0], 'Flag_Feedback':FlagFeedback[0], 'NumPart_Total':Nall, \
    'Flag_Cooling':FlagCooling[0], 'NumFilesPerSnapshot':NumFiles[0], 'BoxSize':BoxSize[0], \
    'Omega0':Omega0[0], 'OmegaLambda':OmegaLambda[0], 'HubbleParam':h, \
    'Flag_StellarAge':FlagAge[0], 'Flag_Metals':FlagMetals[0], 'Nall_HW':NallHW, \
    'Flag_EntrICs':flag_entr_ics[0]}



def load_gadget_binary_particledat(f, header, ptype, skip_bh=0):
    ## load old format=1 style gadget binary snapshot files (unformatted fortran binary)
    ### COULD THIS BE THE PROBLEM?? FORMAT 1 versus format 2???
    ###########################
    #
    # WARNING about recovering density, internal energy, smoothing length
    # if this yells at you, IC may not have been generated properly
    #
    #############################

    import array
    gas_u=0.; gas_rho=0.; gas_ne=0.; gas_nhi=0.; gas_hsml=0.; gas_SFR=0.; star_age=0.; 
    zmet=0.; bh_mass=0.; bh_mdot=0.; mm=0.;
    Npart = header['NumPart_ThisFile']



    Massarr = header['MassTable']


    NpartTot = np.sum(Npart)
    NpartCum = np.cumsum(Npart)
    n0 = NpartCum[ptype] - Npart[ptype]
    n1 = NpartCum[ptype]
    
    ### particles positions. 3*Npart*float.
    pos = array.array('d')
    pos.fromfile(f, 3*NpartTot)
    pos = np.reshape(pos, (NpartTot,3))
    f.read(4+4) # Read block size fields.

    ### particles velocities. 3*Npart*float.
    vel = array.array('d')
    vel.fromfile(f, 3*NpartTot)
    vel = np.reshape(vel, (NpartTot,3))
    f.read(4+4) # Read block size fields.

    ### Particle IDs. # (Npart[0]+...+Npart[5])*int
    ids = array.array('i')
    ids.fromfile(f, NpartTot)
    ids = np.array(ids)
    f.read(4+4) # Read block size fields.
        
    ### Variable particle masses. 
    Npart_MassCode = np.copy(np.array(Npart))
    Npart=np.array(Npart)
    Npart_MassCode[(Npart <= 0) | (np.array(Massarr,dtype='d') > 0.0)] = 0
    NwithMass = np.sum(Npart_MassCode)
    mass = array.array('d')

    #print "\n\nmass loaded in binary reader: ", mass

    mass.fromfile(f, NwithMass)
    f.read(4+4) # Read block size fields.
    if (Massarr[ptype]==0.0):
        Npart_MassCode_Tot = np.cumsum(Npart_MassCode)
        mm = mass[Npart_MassCode_Tot[ptype]-Npart_MassCode[ptype]:Npart_MassCode_Tot[ptype]]

    if ((ptype==0) | (ptype==4) | (ptype==5)):
        if (Npart[0]>0):
            ### Internal energy of gas particles ((km/s)^2).

            f.read(4+4) # Read block size fields.
            ### Density for the gas paraticles (units?).

            gas_u = array.array('d')
            gas_rho = array.array('d')
            gas_hsml = array.array('d')


            ############## EDITING HERE NOW #######################
            #### edit edit edit ###

            try:
                gas_u.fromfile(f, Npart[0]) ### <---- EDIT maybe this is the keyword for internal energy????
                gas_rho.fromfile(f, Npart[0])
                #print 'gas_u.whatever', gas_u.fromfile(f, Npart[0])

                ### Smoothing length (kpc/h). ###
                gas_hsml.fromfile(f, Npart[0])
            except EOFError:
                # gas_rho
                print "WARNING: Error loading internal energy OR density OR smoothing length"
                #print "Npart, Npart[0]: ", Npart, Npart[0]
                pass

            f.read(4+4) # Read block size fields.

            if (header['Flag_Cooling'] > 0):
                ### Electron number density for gas particles (fraction of n_H; can be >1).
                gas_ne = array.array('d')
                gas_ne.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.
                ### Neutral hydrogen number density for gas particles (fraction of n_H).
                gas_nhi = array.array('d')
                gas_nhi.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.

            f.read(4+4) # Read block size fields.

            if (header['Flag_Sfr'] > 0):
                ### Star formation rate (Msun/yr). ###
                gas_SFR = array.array('d')
                gas_SFR.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.

        if (Npart[4]>0):
            if (header['Flag_Sfr'] > 0):
                if (header['Flag_StellarAge'] > 0):
                    ### Star formation time (in code units) or scale factor ###
                    star_age = array.array('d')
                    star_age.fromfile(f, Npart[4])
                    f.read(4+4) # Read block size fields.
        
        if (Npart[0]+Npart[4]>0):
            if (header['Flag_Metals'] > 0):
                ## Metallicity block (species tracked = Flag_Metals)
                if (Npart[0]>0):
                    gas_z = array.array('d')
                    gas_z.fromfile(f, header['Flag_Metals']*Npart[0])
                if (Npart[4]>0):
                    star_z = array.array('d')
                    star_z.fromfile(f, header['Flag_Metals']*Npart[4])
                f.read(4+4) # Read block size fields.
                if (ptype==0): zmet=np.reshape(gas_z,(-1,header['Flag_Metals']))
                if (ptype==4): zmet=np.reshape(star_z,(-1,header['Flag_Metals']))
        
        if (Npart[5]>0):
            if (skip_bh > 0):
                ## BH mass (same as code units, but this is the separately-tracked BH mass from particle mass)
                bh_mass = array.array('d')
                bh_mass.fromfile(f, Npart[5])
                f.read(4+4) # Read block size fields.
                ## BH accretion rate in snapshot
                bh_mdot = array.array('d')
                bh_mdot.fromfile(f, Npart[5])
                f.read(4+4) # Read block size fields.
    
    return {'Coordinates':pos[n0:n1,:], 'Velocities':vel[n0:n1,:], 'ParticleIDs':ids[n0:n1], \
        'Masses':mm, 'Metallicity':zmet, 'StellarFormationTime':star_age, 'BH_Mass':bh_mass, \
        'BH_Mdot':bh_mdot, 'InternalEnergy':gas_u, 'Density':gas_rho, 'SmoothingLength':gas_hsml, \
        'ElectronAbundance':gas_ne, 'NeutralHydrogenAbundance':gas_nhi, 'StarFormationRate':gas_SFR}


###############################################################
#
# binary writing routine
#
###############################################################
def make_IC_binary(fname, mp, x, y, z, E, **kwargs):
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
  
    ## trying to load this from MESA!
    my_body.u[:]=E 

    #U + 0.*x#internal_energy ### <--------- load this from MESA directly

    pygadgetic.dump_ic(my_header,my_body,fname)
    return fname




#####################################################################
#
# hdf5 routines
#
#####################################################################
class snapshot_header:
	def __init__(self, *args, **kwargs):
		if (len(args) == 1):
			filename = args[0]

			if os.path.exists(filename):
				curfilename=filename
			elif os.path.exists(filename+".hdf5"):
				curfilename = filename+".hdf5"
			elif os.path.exists(filename+".0.hdf5"): 
				curfilename = filename+".0.hdf5"
			else:	
				print( "[error] file not found : ", filename)
				sys.stdout.flush()
				sys.exit()

			f=hdf5lib.OpenFile(curfilename)
			self.npart = hdf5lib.GetAttr(f, "Header", "NumPart_ThisFile") 
			self.nall = hdf5lib.GetAttr(f, "Header", "NumPart_Total")
			self.nall_highword = hdf5lib.GetAttr(f, "Header", "NumPart_Total_HighWord") 
			self.massarr = hdf5lib.GetAttr(f, "Header", "MassTable")
			self.time = hdf5lib.GetAttr(f, "Header", "Time") 
			self.redshift = hdf5lib.GetAttr(f, "Header", "Redshift") 
			self.boxsize = hdf5lib.GetAttr(f, "Header", "BoxSize") 
			self.filenum = hdf5lib.GetAttr(f, "Header", "NumFilesPerSnapshot") 
			self.omega0 = hdf5lib.GetAttr(f, "Header", "Omega0")
			self.omegaL = hdf5lib.GetAttr(f, "Header", "OmegaLambda") 
			self.hubble = hdf5lib.GetAttr(f, "Header", "HubbleParam") 
			self.sfr = hdf5lib.GetAttr(f, "Header", "Flag_Sfr") 
			self.cooling = hdf5lib.GetAttr(f, "Header", "Flag_Cooling") 
			self.stellar_age = hdf5lib.GetAttr(f, "Header", "Flag_StellarAge") 
			self.metals = hdf5lib.GetAttr(f, "Header", "Flag_Metals") 
			self.feedback = hdf5lib.GetAttr(f, "Header", "Flag_Feedback") 
			self.double = hdf5lib.GetAttr(f, "Header", "Flag_DoublePrecision") #GADGET-2 change
			f.close()
		else:
			#read arguments
			self.npart = kwargs.get("npart")
			self.nall = kwargs.get("nall")
			self.nall_highword = kwargs.get("nall_highword")
			self.massarr = kwargs.get("massarr")
			self.time = kwargs.get("time")
			self.redshift = kwargs.get("redshift")
			self.boxsize = kwargs.get("boxsize")
			self.filenum = kwargs.get("filenum")
			self.omega0 = kwargs.get("omega0")
			self.omegaL = kwargs.get("omegaL")
			self.hubble = kwargs.get("hubble")
			self.sfr = kwargs.get("sfr")
			self.cooling = kwargs.get("cooling")
			self.stellar_age = kwargs.get("stellar_age")
			self.metals = kwargs.get("metals")
			self.feedback = kwargs.get("feedback")
			self.double = kwargs.get("double")

			#set default values
			if (self.npart == None):
				self.npart = np.array([0,0,0,0,0,0], dtype="int32")
			if (self.nall == None):				
				self.nall  = np.array([0,0,0,0,0,0], dtype="uint32")
			if (self.nall_highword == None):				
				self.nall_highword = np.array([0,0,0,0,0,0], dtype="uint32")
			if (self.massarr == None):
				self.massarr = np.array([0,0,0,0,0,0], dtype="float64")
			if (self.time == None):				
				self.time = np.array([0], dtype="float64")
			if (self.redshift == None):				
				self.redshift = np.array([0], dtype="float64")
			if (self.boxsize == None):				
				self.boxsize = np.array([0], dtype="float64")
			if (self.filenum == None):
				self.filenum = np.array([1], dtype="int32")
			if (self.omega0 == None):
				self.omega0 = np.array([0], dtype="float64")
			if (self.omegaL == None):
				self.omegaL = np.array([0], dtype="float64")
			if (self.hubble == None):
				self.hubble = np.array([0], dtype="float64")
			if (self.sfr == None):	
				self.sfr = np.array([0], dtype="int32")            
			if (self.cooling == None):	
				self.cooling = np.array([0], dtype="int32")
			if (self.stellar_age == None):	
				self.stellar_age = np.array([0], dtype="int32")
			if (self.metals == None):	
				self.metals = np.array([0], dtype="int32")
			if (self.feedback == None):	
				self.feedback = np.array([0], dtype="int32")
			if (self.double == None):
				self.double = np.array([0], dtype="int32")


def read_block_single_file(filename, block_name, dim2, parttype=-1, no_mass_replicate=False,\
 fill_block_name="", slab_start=-1, slab_len=-1, verbose=False):

	if (verbose):
		print( "[single] reading file           : ", filename   )	
		print( "[single] reading                : ", block_name )
		sys.stdout.flush()
	# error here for some reason	
	head = snapshot_header(filename)
	npart = head.npart
	massarr = head.massarr
	nall = head.nall
	filenum = head.filenum
	doubleflag = head.double #GADGET-2 change
	#doubleflag = 0          #GADGET-2 change

	if (parttype!=-1):
		if (head.npart[parttype]==0):
			return [0, False]
	else:
		if (head.npart.sum()==0):
			return [0, False]
	del head


	#construct indices for partial access
	if (slab_start!=-1) & (slab_len!=-1):
		data_slice = slice(slab_start, (slab_start+slab_len))
	else:
		data_slice = slice(None, None, 1)

	if (verbose):
		print( "[single] data slice: ", data_slice)
		sys.stdout.flush() 

	f=hdf5lib.OpenFile(filename)

	#read specific particle type (parttype>=0, non-default)
	if parttype>=0:
		if (verbose):
			print( "[single] parttype               : ", parttype )
			sys.stdout.flush()
		if ((block_name=="Masses") & (npart[parttype]>0) & (massarr[parttype]>0)):
			if (verbose):
				print( "[single] replicate mass block")
				sys.stdout.flush()	
			ret_val=np.repeat(massarr[parttype], npart[parttype])[data_slice]
		else:		
			part_name='PartType'+str(parttype)
			ret_val = hdf5lib.GetData(f, part_name+"/"+block_name)[data_slice]
		if (verbose):
			print( "[single] read particles (total) : ", ret_val.shape[0]/dim2)
			sys.stdout.flush()

	#read all particle types (parttype=-1, default); slab reading not possible here; i.e. data_slice not used below
	if parttype==-1:
		first=True
		dim1=long(0)
		for parttype in range(0,6):
			part_name='PartType'+str(parttype)
			if (hdf5lib.Contains(f,"",part_name)):
				if (verbose):
					print( "[single] parttype               : ", parttype)
					print( "[single] massarr                : ", massarr)
					print( "[single] npart                  : ", npart)
					sys.stdout.flush()

				#replicate mass block per default (unless no_mass_replicate is set)
				if ((block_name=="Masses") & (npart[parttype]>0) & (massarr[parttype]>0) & (no_mass_replicate==False)):
					if (verbose):
						print( "[single] replicate mass block")
						sys.stdout.flush()
					if (first):
						data=np.repeat(massarr[parttype], npart[parttype])
						dim1+=long(data.shape[0])
						ret_val=data
						first=False
					else:	
						data=np.repeat(massarr[parttype], npart[parttype])
						dim1+=long(data.shape[0])
						ret_val=np.append(ret_val, data)
					if (verbose):
						print( "[single] read particles (total) : ", ret_val.shape[0]/dim2)
						sys.stdout.flush()
					if (doubleflag==0):
						ret_val=ret_val.astype("float32")


				#fill fill_block_name with zeros if fill_block_name is set and particle type is present and fill_block_name not already stored in file for that particle type
				if ((block_name==fill_block_name) & (block_name!="Masses") & (npart[parttype]>0) & (hdf5lib.Contains(f,part_name, block_name)==False)):
					if (verbose):
						print( "[single] replicate block name", fill_block_name)
						sys.stdout.flush()
					if (first):
						data=np.repeat(0.0, npart[parttype]*dim2)
						dim1+=long(data.shape[0])
						ret_val=data
						first=False
					else:
						data=np.repeat(0.0, npart[parttype]*dim2)
						dim1+=long(data.shape[0])
						ret_val=np.append(ret_val, data)
					if (verbose):
						print( "[single] read particles (total) : ", ret_val.shape[0]/dim2)
						sys.stdout.flush()
					if (doubleflag==0):
						ret_val=ret_val.astype("float32")

				#default: just read the block
				if (hdf5lib.Contains(f,part_name,block_name)):
					if (first):
						data=hdf5lib.GetData(f, part_name+"/"+block_name)[:]
						dim1+=long(data.shape[0])
						ret_val=data
						first=False
					else:
						data=hdf5lib.GetData(f, part_name+"/"+block_name)[:]		
						dim1+=long(data.shape[0])
						ret_val=np.append(ret_val, data)
					if (verbose):
						print( "[single] read particles (total) : ", ret_val.shape[0]/dim2)
						sys.stdout.flush()

		if ((dim1>0) & (dim2>1)):
			ret_val=ret_val.reshape(dim1,dim2)
			if (verbose):
				print( "[single] reshape done: ", ret_val.shape)
				sys.stdout.flush()
	f.close()
	return [ret_val, True]



###############################################################
#
# hdf5 writing routine
#
###############################################################
def make_IC_hdf5(out_fname, mp, x, y, z,E, **kwargs):

    print "using updated hdf5 writer"

    file = h5py.File(out_fname,'w') 
    Ngas = len(x)
    npart = np.array([len(x),0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others

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
   #internal energy 
    #particle IDs

    particles.create_dataset("ParticleIDs",data=IDs)
    particles.create_dataset("Masses",data=masses,dtype=np.dtype('d'))
    particles.create_dataset("Coordinates",data=pos,dtype=np.dtype('d'))
    particles.create_dataset("Velocities",data=vel,dtype=np.dtype('d'))
    particles.create_dataset("InternalEnergy",data=internalE,dtype=np.dtype('d'))

    particles.create_dataset("HSML", data=hsml.astype(float))
    particles.create_dataset("Density",data=dens.astype(float))

    file.close()
    return file
