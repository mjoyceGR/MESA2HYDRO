#!/usr/bin/env python
import numpy as np
import h5py as h5py
import os.path
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import math
import sys 

##########################################
#
# read_snap_hdf5.py
#
#########################################
## This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO ##
def readsnap(sdir,snum,ptype,
    snapshot_name='snapshot',
    extension='.hdf5',
    h0=0,cosmological=0,skip_bh=0,four_char=0,
    header_only=0,loud=0):
    
    if (ptype<0): return {'k':-1};
    if (ptype>5): return {'k':-1};

    fname,fname_base,fname_ext = check_if_filename_exists(sdir,snum,\
        snapshot_name=snapshot_name,extension=extension,four_char=four_char)
    if(fname=='NULL'): return {'k':-1}
    if(loud==1): print('loading file : '+fname)

    ## open file and parse its header information
    nL = 0 # initial particle point to start at 
    if(fname_ext=='.hdf5'):
        file = h5py.File(fname,'r') # Open hdf5 snapshot file
        header_master = file["Header"] # Load header dictionary (to parse below)
        header_toparse = header_master.attrs
    else:
        file = open(fname) # Open binary snapshot file
        header_toparse = load_gadget_binary_header(file)

    npart = header_toparse["NumPart_ThisFile"]
    massarr = header_toparse["MassTable"]
    time = header_toparse["Time"]
    redshift = header_toparse["Redshift"]
    flag_sfr = header_toparse["Flag_Sfr"]
    flag_feedbacktp = header_toparse["Flag_Feedback"]
    npartTotal = header_toparse["NumPart_Total"]
    flag_cooling = header_toparse["Flag_Cooling"]
    numfiles = header_toparse["NumFilesPerSnapshot"]
    boxsize = header_toparse["BoxSize"]
    omega_matter = header_toparse["Omega0"]
    omega_lambda = header_toparse["OmegaLambda"]
    hubble = header_toparse["HubbleParam"]
    flag_stellarage = header_toparse["Flag_StellarAge"]
    flag_metals = header_toparse["Flag_Metals"]
    print("npart_file: ",npart)
    print("npart_total:",npartTotal)

    hinv=1.
    if (h0==1):
        hinv=1./hubble
    ascale=1.
    if (cosmological==1):
        ascale=time
        hinv=1./hubble
    if (cosmological==0): 
        time*=hinv
    
    boxsize*=hinv*ascale
    if (npartTotal[ptype]<=0): file.close(); return {'k':-1};
    if (header_only==1): file.close(); return {'k':0,'time':time,
        'boxsize':boxsize,'hubble':hubble,'npart':npart,'npartTotal':npartTotal};

    # initialize variables to be read
    pos=np.zeros([npartTotal[ptype],3],dtype=np.float64)
    vel=np.copy(pos)
    ids=np.zeros([npartTotal[ptype]],dtype=int)
    mass=np.zeros([npartTotal[ptype]],dtype=np.float64)
    if (ptype==0):
        ugas=np.copy(mass)
        rho=np.copy(mass)
        hsml=np.copy(mass) 
        #if (flag_cooling>0): 
        nume=np.copy(mass)
        numh=np.copy(mass)
        #if (flag_sfr>0): 
        sfr=np.copy(mass)
        metal=np.copy(mass)
    if (ptype==0 or ptype==4) and (flag_metals > 0):
        metal=np.zeros([npartTotal[ptype],flag_metals],dtype=np.float64)
    if (ptype==4) and (flag_sfr>0) and (flag_stellarage>0):
        stellage=np.copy(mass)
    if (ptype==5) and (skip_bh==0):
        bhmass=np.copy(mass)
        bhmdot=np.copy(mass)

    # loop over the snapshot parts to get the different data pieces
    for i_file in range(numfiles):
        if (numfiles>1):
            file.close()
            fname = fname_base+'.'+str(i_file)+fname_ext
            if(fname_ext=='.hdf5'):
                file = h5py.File(fname,'r') # Open hdf5 snapshot file
            else:
                file = open(fname) # Open binary snapshot file
                header_toparse = load_gadget_binary_header(file)
                
        if (fname_ext=='.hdf5'):
            input_struct = file
            npart = file["Header"].attrs["NumPart_ThisFile"]
            bname = "PartType"+str(ptype)+"/"
        else:
            npart = header_toparse['NumPart_ThisFile']
            input_struct = load_gadget_binary_particledat(file, header_toparse, ptype, skip_bh=skip_bh)
            bname = ''
            
        
        # now do the actual reading
        if(npart[ptype]>0):
            nR=nL + npart[ptype]
            pos[nL:nR,:]=input_struct[bname+"Coordinates"]
            vel[nL:nR,:]=input_struct[bname+"Velocities"]
            ids[nL:nR]=input_struct[bname+"ParticleIDs"]
            mass[nL:nR]=massarr[ptype]
            if (massarr[ptype] <= 0.):
                mass[nL:nR]=input_struct[bname+"Masses"]
            if (ptype==0):
                ugas[nL:nR]=input_struct[bname+"InternalEnergy"]
                rho[nL:nR]=input_struct[bname+"Density"]
                hsml[nL:nR]=input_struct[bname+"SmoothingLength"]
                if (flag_cooling > 0): 
                    nume[nL:nR]=input_struct[bname+"ElectronAbundance"]
                    numh[nL:nR]=input_struct[bname+"NeutralHydrogenAbundance"]
                if (flag_sfr > 0):
                    sfr[nL:nR]=input_struct[bname+"StarFormationRate"]
            if (ptype==0 or ptype==4) and (flag_metals > 0):
                metal_t=input_struct[bname+"Metallicity"]
                if (flag_metals > 1):
                    if (metal_t.shape[0] != npart[ptype]): 
                        metal_t=np.transpose(metal_t)
                else:
                    metal_t=np.reshape(np.array(metal_t),(np.array(metal_t).size,1))
                metal[nL:nR,:]=metal_t
            if (ptype==4) and (flag_sfr>0) and (flag_stellarage>0):
                stellage[nL:nR]=input_struct[bname+"StellarFormationTime"]
            if (ptype==5) and (skip_bh==0):
                bhmass[nL:nR]=input_struct[bname+"BH_Mass"]
                bhmdot[nL:nR]=input_struct[bname+"BH_Mdot"]
            nL = nR # sets it for the next iteration	

	## correct to same ID as original gas particle for new stars, if bit-flip applied
    if ((np.min(ids)<0) | (np.max(ids)>1.e9)):
        bad = (ids < 0) | (ids > 1.e9)
        ids[bad] += (1 << 31)

    # do the cosmological conversions on final vectors as needed
    pos *= hinv*ascale # snapshot units are comoving
    mass *= hinv
    vel *= np.sqrt(ascale) # remember gadget's weird velocity units!
    if (ptype==0):
        rho *= (hinv/((ascale*hinv)**3))
        hsml *= hinv*ascale
    if (ptype==4) and (flag_sfr>0) and (flag_stellarage>0) and (cosmological==0):
        stellage *= hinv
    if (ptype==5) and (skip_bh==0):
        bhmass *= hinv

    file.close();
    if (ptype==0):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'u':ugas,'rho':rho,'h':hsml,'ne':nume,'nh':numh,'sfr':sfr,'z':metal};
    if (ptype==4):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'z':metal,'age':stellage}
    if (ptype==5) and (skip_bh==0):
        return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids,'mbh':bhmass,'mdot':bhmdot}
    return {'k':1,'p':pos,'v':vel,'m':mass,'id':ids}



def check_if_filename_exists(sdir,snum,snapshot_name='snapshot',extension='.hdf5',four_char=0):
    for extension_touse in [extension,'.bin','']:
        fname=sdir+'/'+snapshot_name+'_'
        ext='00'+str(snum);
        if (snum>=10): ext='0'+str(snum)
        if (snum>=100): ext=str(snum)
        if (four_char==1): ext='0'+ext
        if (snum>=1000): ext=str(snum)
        fname+=ext
        fname_base=fname

        s0=sdir.split("/"); snapdir_specific=s0[len(s0)-1];
        if(len(snapdir_specific)<=1): snapdir_specific=s0[len(s0)-2];

        ## try several common notations for the directory/filename structure
        fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is it a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot'?
            fname_base=sdir+'/snap_'+ext; 
            fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot', AND its a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap(snapdir)' instead of 'snapshot'?
            fname_base=sdir+'/snap_'+snapdir_specific+'_'+ext; 
            fname=fname_base+extension_touse;
        if not os.path.exists(fname): 
            ## is the filename 'snap' instead of 'snapshot', AND its a multi-part file?
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory? (we assume this means multi-part files)
            fname_base=sdir+'/snapdir_'+ext+'/'+snapshot_name+'_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## is it in a snapshot sub-directory AND named 'snap' instead of 'snapshot'?
            fname_base=sdir+'/snapdir_'+ext+'/'+'snap_'+ext; 
            fname=fname_base+'.0'+extension_touse;
        if not os.path.exists(fname): 
            ## wow, still couldn't find it... ok, i'm going to give up!
            fname_found = 'NULL'
            fname_base_found = 'NULL'
            fname_ext = 'NULL'
            continue;
        fname_found = fname;
        fname_base_found = fname_base;
        fname_ext = extension_touse
        break; # filename does exist! 
    return fname_found, fname_base_found, fname_ext;



def load_gadget_binary_header(f):
    ### Read header.
    import array
    # Skip 4-byte integer at beginning of header block.
    f.read(4)
    # Number of particles of each type. 6*unsigned integer.
    Npart = array.array('I')
    Npart.fromfile(f, 6)
    # Mass of each particle type. If set to 0 for a type which is present, 
    # individual particle masses from the 'mass' block are used instead.
    # 6*double.
    Massarr = array.array('d')
    Massarr.fromfile(f, 6)
    # Expansion factor (or time, if non-cosmological sims) of output. 1*double. 
    a = array.array('d')
    a.fromfile(f, 1)
    a = a[0]
    # Redshift of output. Should satisfy z=1/a-1. 1*double.
    z = array.array('d')
    z.fromfile(f, 1)
    z = float(z[0])
    # Flag for star formation. 1*int.
    FlagSfr = array.array('i')
    FlagSfr.fromfile(f, 1)
    # Flag for feedback. 1*int.
    FlagFeedback = array.array('i')
    FlagFeedback.fromfile(f, 1)
    # Total number of particles of each type in the simulation. 6*int.
    Nall = array.array('i')
    Nall.fromfile(f, 6)
    # Flag for cooling. 1*int.
    FlagCooling = array.array('i')
    FlagCooling.fromfile(f, 1)
    # Number of files in each snapshot. 1*int.
    NumFiles = array.array('i')
    NumFiles.fromfile(f, 1)
    # Box size (comoving kpc/h). 1*double.
    BoxSize = array.array('d')
    BoxSize.fromfile(f, 1)
    # Matter density at z=0 in units of the critical density. 1*double.
    Omega0 = array.array('d')
    Omega0.fromfile(f, 1)
    # Vacuum energy density at z=0 in units of the critical density. 1*double.
    OmegaLambda = array.array('d')
    OmegaLambda.fromfile(f, 1)
    # Hubble parameter h in units of 100 km s^-1 Mpc^-1. 1*double.
    h = array.array('d')
    h.fromfile(f, 1)
    h = float(h[0])
    # Creation times of stars. 1*int.
    FlagAge = array.array('i')
    FlagAge.fromfile(f, 1)
    # Flag for metallicity values. 1*int.
    FlagMetals = array.array('i')
    FlagMetals.fromfile(f, 1)

    # For simulations that use more than 2^32 particles, most significant word 
    # of 64-bit total particle numbers. Otherwise 0. 6*int.
    NallHW = array.array('i')
    NallHW.fromfile(f, 6)

    # Flag that initial conditions contain entropy instead of thermal energy
    # in the u block. 1*int.
    flag_entr_ics = array.array('i')
    flag_entr_ics.fromfile(f, 1)

    # Unused header space. Skip to particle positions.
    f.seek(4+256+4+4)

    return {'NumPart_ThisFile':Npart, 'MassTable':Massarr, 'Time':a, 'Redshift':z, \
    'Flag_Sfr':FlagSfr[0], 'Flag_Feedback':FlagFeedback[0], 'NumPart_Total':Nall, \
    'Flag_Cooling':FlagCooling[0], 'NumFilesPerSnapshot':NumFiles[0], 'BoxSize':BoxSize[0], \
    'Omega0':Omega0[0], 'OmegaLambda':OmegaLambda[0], 'HubbleParam':h, \
    'Flag_StellarAge':FlagAge[0], 'Flag_Metals':FlagMetals[0], 'Nall_HW':NallHW, \
    'Flag_EntrICs':flag_entr_ics[0]}


def load_gadget_binary_particledat(f, header, ptype, skip_bh=0):
    ## load old format=1 style gadget binary snapshot files (unformatted fortran binary)
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
    pos = array.array('f')
    pos.fromfile(f, 3*NpartTot)
    pos = np.reshape(pos, (NpartTot,3))
    f.read(4+4) # Read block size fields.

    ### particles velocities. 3*Npart*float.
    vel = array.array('f')
    vel.fromfile(f, 3*NpartTot)
    vel = np.reshape(vel, (NpartTot,3))
    f.read(4+4) # Read block size fields.

    ### Particle IDs. # (Npart[0]+...+Npart[5])*int
    id = array.array('i')
    id.fromfile(f, NpartTot)
    id = np.array(id)
    f.read(4+4) # Read block size fields.
        
    ### Variable particle masses. 
    Npart_MassCode = np.copy(np.array(Npart))
    Npart=np.array(Npart)
    Npart_MassCode[(Npart <= 0) | (np.array(Massarr,dtype='d') > 0.0)] = 0
    NwithMass = np.sum(Npart_MassCode)
    mass = array.array('f')
    mass.fromfile(f, NwithMass)
    f.read(4+4) # Read block size fields.
    if (Massarr[ptype]==0.0):
        Npart_MassCode_Tot = np.cumsum(Npart_MassCode)
        mm = mass[Npart_MassCode_Tot[ptype]-Npart_MassCode[ptype]:Npart_MassCode_Tot[ptype]]

    if ((ptype==0) | (ptype==4) | (ptype==5)):
        if (Npart[0]>0):
            ### Internal energy of gas particles ((km/s)^2).
            gas_u = array.array('f')
            gas_u.fromfile(f, Npart[0])
            f.read(4+4) # Read block size fields.
            ### Density for the gas paraticles (units?).
            gas_rho = array.array('f')
            gas_rho.fromfile(f, Npart[0])
            f.read(4+4) # Read block size fields.

            if (header['Flag_Cooling'] > 0):
                ### Electron number density for gas particles (fraction of n_H; can be >1).
                gas_ne = array.array('f')
                gas_ne.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.
                ### Neutral hydrogen number density for gas particles (fraction of n_H).
                gas_nhi = array.array('f')
                gas_nhi.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.

            ### Smoothing length (kpc/h). ###
            gas_hsml = array.array('f')
            gas_hsml.fromfile(f, Npart[0])
            f.read(4+4) # Read block size fields.

            if (header['Flag_Sfr'] > 0):
                ### Star formation rate (Msun/yr). ###
                gas_SFR = array.array('f')
                gas_SFR.fromfile(f, Npart[0])
                f.read(4+4) # Read block size fields.

        if (Npart[4]>0):
            if (header['Flag_Sfr'] > 0):
                if (header['Flag_StellarAge'] > 0):
                    ### Star formation time (in code units) or scale factor ###
                    star_age = array.array('f')
                    star_age.fromfile(f, Npart[4])
                    f.read(4+4) # Read block size fields.
        
        if (Npart[0]+Npart[4]>0):
            if (header['Flag_Metals'] > 0):
                ## Metallicity block (species tracked = Flag_Metals)
                if (Npart[0]>0):
                    gas_z = array.array('f')
                    gas_z.fromfile(f, header['Flag_Metals']*Npart[0])
                if (Npart[4]>0):
                    star_z = array.array('f')
                    star_z.fromfile(f, header['Flag_Metals']*Npart[4])
                f.read(4+4) # Read block size fields.
                if (ptype==0): zmet=np.reshape(gas_z,(-1,header['Flag_Metals']))
                if (ptype==4): zmet=np.reshape(star_z,(-1,header['Flag_Metals']))
        
        if (Npart[5]>0):
            if (skip_bh > 0):
                ## BH mass (same as code units, but this is the separately-tracked BH mass from particle mass)
                bh_mass = array.array('f')
                bh_mass.fromfile(f, Npart[5])
                f.read(4+4) # Read block size fields.
                ## BH accretion rate in snapshot
                bh_mdot = array.array('f')
                bh_mdot.fromfile(f, Npart[5])
                f.read(4+4) # Read block size fields.
    
    return {'Coordinates':pos[n0:n1,:], 'Velocities':vel[n0:n1,:], 'ParticleIDs':id[n0:n1], \
        'Masses':mm, 'Metallicity':zmet, 'StellarFormationTime':star_age, 'BH_Mass':bh_mass, \
        'BH_Mdot':bh_mdot, 'InternalEnergy':gas_u, 'Density':gas_rho, 'SmoothingLength':gas_hsml, \
        'ElectronAbundance':gas_ne, 'NeutralHydrogenAbundance':gas_nhi, 'StarFormationRate':gas_SFR}


##########################################
#
# read_write_HDF5.py
#
#########################################
############ 
#DATABLOCKS#
############
#descriptions of all datablocks -> add new datablocks here!
#format -> TAG:[HDF5_NAME, DIMENSION]
datablocks = { 	"POS ":["Coordinates",3], 
		"VEL ":["Velocities",3],
		"ID  ":["ParticleIDs",1],
		"MASS":["Masses",1],
		"U   ":["InternalEnergy",1],
		"RHO ":["Density",1],
		"VOL ":["Volume",1],
		"CMCE":["Center-of-Mass",3],
		"AREA":["Surface Area",1],
		"NFAC":["Number of faces of cell",1],
		"NE  ":["ElectronAbundance",1],
		"NH  ":["NeutralHydrogenAbundance",1],
		"HSML":["SmoothingLength",1],
		"SFR ":["StarFormationRate",1],
		"AGE ":["StellarFormationTime",1],
		"Z   ":["Metallicity",1],
		"ACCE":["Acceleration",3],
		"VEVE":["VertexVelocity",3],
		"FACA":["MaxFaceAngle",1],              
		"COOR":["CoolingRate",1],
		"POT ":["Potential",1],
		"MACH":["MachNumber",1],
		"DMHS":["DM Hsml",1],
		"DMDE":["DM Density",1],
		"PHKE":["PHKey",1],
		"GROU":["GroupNr",1],
		#SIDM (GADGET)
		"PTSU":["PSum",1],
		"DMNB":["DMNumNgb",1],
		"NTSC":["NumTotalScatter",1],
		"SHSM":["SIDMHsml",1],
		"SRHO":["SIDMRho",1],
		"SVEL":["SVelDisp",1],
		#SIDM (AREPO)
		"PTSU":["SIDM_Psum",1],
		"DMNB":["SIDM_NumNgb",1],
		"NTSC":["SIDM_NumTotalScatter",1],
		"SHSM":["SIDM_Hsml",1],
		"SRHO":["SIDM_Density",1],
		"SVEL":["SIDM_VelDisp",1],
		#TRACER
                "TRCE":["TracerField", 1],
                #TRACERMC
                "TRNT":["NumTracers", 1],       #parttypes: 0,4,5
                "TRFQ":["FluidQuantities", 13], #parttype: 3
	        "TRID":["TracerID", 1],         #parttype: 3
        	"TRPR":["ParentID", 1],         #parttype: 3
		#GFM 
		"GAGE":["GFM_StellarFormationTime",1],
		"GIMA":["GFM_InitialMass",1],
		"GZ  ":["GFM_Metallicity",1],
		"GMET":["GFM_Metals",9],
		"GWHV":["GFM_WindHostVal",1],
		"GCOL":["GFM_CoolingRate",1],
		"GSPH":["GFM_StellarPhotometrics",8], #band luminosities: U, B, V, K, g, r, i, z
		"AGNR":["GFM_AGNRadiation",1],
		#GDE
		"CACO":["CausticCounter",1],
		"STDE":["StreamDensity",1],
		"PSDE":["PhaseSpaceDensity",1],
		"FLDE":["FlowDeterminant",1],
		"TIPS":["TidalTensorPS",9],
		"DIPS":["DistortionTensorPS",36],
		"SHOR":["SheetOrientation", 9],
		"INDE":["InitDensity",1],
		#ISM
		"BRD ":["BlastRadius", 1],
		"CSO ":["CoolShutoffTime", 1],
		"RadP":["RadiationPressureMoment",1], 
		#SUBFIND
		"SFDE":["SubfindDensity", 1],
		"SFHS":["SubfindHsml", 1],
		"SFVD":["SubfindVelDisp", 1],
		#BHs
		"BHMD":["BH_Mdot", 1],
		"BHHS":["BH_Hsml", 1],
		"BHMA":["BH_Mass", 1],
		"REBH":["RefBHCounter", 1],
		"BHHG":["BH_HaloGasMass", 1],
		"BHCL":["BH_CoolingLuminosity", 1],
		"BHMR":["BH_Mdot_Radio", 1], 
		"BHPR":["BH_Progs", 1],
                "BCEQ":["BH_CumEgyInjection_QM",1],
                "BHMB":["BH_Mass_bubbles",1]
             }
		

#####################################################################################################################
#                                                    READING ROUTINES			                            #
#####################################################################################################################


########################### 
#CLASS FOR SNAPSHOT HEADER#
###########################  
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




##############################
#READ ROUTINE FOR SINGLE FILE#
############################## 
def read_block_single_file(filename, block_name, dim2, parttype=-1, no_mass_replicate=False, fill_block_name="", slab_start=-1, slab_len=-1, verbose=False):

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

##############
#READ ROUTINE#
##############
def read_block(filename, block, parttype=-1, no_mass_replicate=False, fill_block="", slab_start=-1, slab_len=-1, verbose=False):
	if (verbose):
		print( "reading block          : ", block)
		sys.stdout.flush()	

	if parttype not in [-1,0,1,2,3,4,5]:
		print( "[error] wrong parttype given")
		sys.stdout.flush()
		sys.exit()

	if os.path.exists(filename):
		curfilename=filename
		multiple_files=False
	elif os.path.exists(filename+".hdf5"):
		curfilename = filename+".hdf5"
		multiple_files=False
	elif os.path.exists(filename+".0.hdf5"):
		curfilename = filename+".0.hdf5"
		multiple_files=True
	else:
		print( "[error] file not found : ", filename)
		sys.stdout.flush()
		sys.exit()

	slabflag=False
	if ((slab_start!=-1) | (slab_len!=-1)):
		slabflag=True
		if (parttype==-1):
			print( "[error] slabs only supported for specific parttype")
			sys.stdout.flush()
			sys.exit()

	head = snapshot_header(curfilename)
	filenum = head.filenum
	del head


	if (datablocks.has_key(block)):
		block_name=datablocks[block][0]
		dim2=datablocks[block][1]
		first=True
		if (verbose):
			print( "Reading HDF5           : ", block_name)
			print( "Data dimension         : ", dim2)
			print( "Multiple file          : ", multiple_files)
			print( "Slab data              : ", slabflag)
			sys.stdout.flush()
	else:
		print( "[error] Block type ", block, "not known!")
		sys.stdout.flush()
		sys.exit()

	fill_block_name=""
	if (fill_block!=""):
		if (datablocks.has_key(fill_block)):
			fill_block_name=datablocks[fill_block][0]
			dim2=datablocks[fill_block][1]
			if (verbose):
				print( "Block filling active   : ", fill_block_name)
				sys.stdout.flush()

	if (multiple_files):	
		if (slabflag==False):
			first=True
			dim1=long(0)
			for num in range(0,filenum):
				curfilename=filename+"."+str(num)+".hdf5"
				if (verbose):
					print( "Reading file           : ", num, curfilename)
					sys.stdout.flush() 
				if (first):
					data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start, slab_len, verbose)
					if (succ == True):
						dim1+=long(data.shape[0])
						ret_val = data
						first = False 
				else:	 
					data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start, slab_len, verbose)
					if (succ == True):
						dim1+=long(data.shape[0])
						ret_val=np.append(ret_val, data)
				if (verbose):
					if (succ):
						print( "Read particles (total) : ", ret_val.shape[0]/dim2)
						sys.stdout.flush()
					else:
						print( "Read particles (total) : none")
						sys.stdout.flush()


		if (slabflag==True):
			off=slab_start
			left=slab_len
			first=True
			dim1=long(0)
			for num in range(0,filenum):
				curfilename=filename+"."+str(num)+".hdf5"
				head = snapshot_header(curfilename)
				nloc=head.npart[parttype]
				if (nloc > off):
					start = off
					if (nloc - off > left):
						count = left
					else:
						count = nloc - off
					if (verbose):
						print( "Reading file           : ", num, curfilename, start, count)
						sys.stdout.flush()
					if (first):
						data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start=start, slab_len=count, verbose=verbose)
						if (succ == True):
							dim1+=data.shape[0]
							ret_val = data
							first = False
					else:
						data, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start=start, slab_len=count, verbose=verbose)
						if (succ == True):
							dim1+=data.shape[0]
							ret_val=np.append(ret_val, data)
					if (verbose):
						if (succ):
							print( "Read particles (total) : ", ret_val.shape[0]/dim2)
							sys.stdout.flush()
						else:
							print( "Read particles (total) : none")
							sys.stdout.flush()

					left -= count
					off += count
				if (left==0):
					break
				off -= nloc



		if (verbose):
			print( "all partial files read in")
			sys.stdout.flush()

		if ((dim1>0) & (dim2>1)):
			ret_val=ret_val.reshape(dim1,dim2)
			if (verbose):
				print( "Reshape done (total): ", ret_val.shape)
				sys.stdout.flush()


	else:	
		ret_val, succ = read_block_single_file(curfilename, block_name, dim2, parttype, no_mass_replicate, fill_block_name, slab_start, slab_len, verbose)

	return ret_val


#############
#LIST BLOCKS#
#############
def list_blocks(filename, parttype=-1, verbose=False):

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
	for parttype in range(0,6):
		part_name='PartType'+str(parttype)
		if (hdf5lib.Contains(f,"",part_name)):
			print( "Parttype contains : ", parttype)
			print( "-------------------")
			sys.stdout.flush()
			iter = it=datablocks.__iter__()
			next = iter.next()
			while (1):
				if (verbose):
					print( "check ", next, datablocks[next][0])
					sys.stdout.flush()
				if (hdf5lib.Contains(f,part_name,datablocks[next][0])):
					print( next, datablocks[next][0])
					sys.stdout.flush()
				try:
					next=iter.next()
				except StopIteration:
					break	
	f.close() 

#################
#CONTAINS BLOCKS#
#################
def contains_block(filename, tag, parttype=-1, verbose=False):

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

	contains_flag=False
	f=hdf5lib.OpenFile(curfilename)
	part_name='PartType'+str(parttype)
	if (hdf5lib.Contains(f,"",part_name)):
		iter = it=datablocks.__iter__()
		next = iter.next()
		while (1):
			if (verbose):
				print( "check ", next, datablocks[next][0])
				sys.stdout.flush()
			if (hdf5lib.Contains(f,part_name,datablocks[next][0])):
				if (next.find(tag)>-1):
					contains_flag=True	
			try:
				next=iter.next()
			except StopIteration:
				break
	f.close() 
	return contains_flag

############
#CHECK FILE#
############
def check_file(filename):
	f=hdf5lib.OpenFile(filename)
	f.close()
		  


#####################################################################################################################
#                                                    WRITING ROUTINES    		                            #
#####################################################################################################################


#######################
#OPEN FILE FOR WRITING#
#######################
def openfile(filename, mode="w"):
	f=hdf5lib.OpenFile(filename, mode = mode)	 
	return f

############
#CLOSE FILE#
############
def closefile(f):
	f.close()

##############################
#WRITE SNAPSHOT HEADER OBJECT#
##############################
def writeheader(f, header):	
	group_header=hdf5lib.CreateGroup(f, "Header")
	hdf5lib.SetAttr(group_header, "NumPart_ThisFile", header.npart)
	hdf5lib.SetAttr(group_header, "NumPart_Total", header.nall)
	hdf5lib.SetAttr(group_header, "NumPart_Total_HighWord", header.nall_highword)
	hdf5lib.SetAttr(group_header, "MassTable", header.massarr)
	hdf5lib.SetAttr(group_header, "Time", header.time)
	hdf5lib.SetAttr(group_header, "Redshift", header.redshift)
	hdf5lib.SetAttr(group_header, "BoxSize", header.boxsize)
	hdf5lib.SetAttr(group_header, "NumFilesPerSnapshot", header.filenum)
	hdf5lib.SetAttr(group_header, "Omega0", header.omega0)
	hdf5lib.SetAttr(group_header, "OmegaLambda", header.omegaL)
	hdf5lib.SetAttr(group_header, "HubbleParam", header.hubble)
	hdf5lib.SetAttr(group_header, "Flag_Sfr", header.sfr)
	hdf5lib.SetAttr(group_header, "Flag_Cooling", header.cooling)
	hdf5lib.SetAttr(group_header, "Flag_StellarAge", header.stellar_age)
	hdf5lib.SetAttr(group_header, "Flag_Metals", header.metals)
	hdf5lib.SetAttr(group_header, "Flag_Feedback", header.feedback)
	hdf5lib.SetAttr(group_header, "Flag_DoublePrecision", header.double)

###############
#WRITE ROUTINE#
###############
def write_block(f, block, parttype, data):   
	##operates on these block objects which are in the dictionary at the top of this thing
	## 9/5/17 meridith editf
	part_name="PartType"+str(parttype)
	if (hdf5lib.Contains(f, "", part_name)==False):
		group=hdf5lib.CreateGroup(f, part_name)
	else:
		group=hdf5lib.GetGroup(f, part_name)	

	if (datablocks.has_key(block)):
		block_name=datablocks[block][0]
		dim2=datablocks[block][1]		
		if (hdf5lib.ContainsGroup(group, block_name)==False):
			hdf5lib.CreateArray(f, group, block_name, data)
		else:
			print( "I/O block already written")
			sys.stdout.flush()
	else:
		print( "Unknown I/O block")
		sys.stdout.flush()		





