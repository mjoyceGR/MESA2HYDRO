############ 
#GADGET HDF5 DATABLOCKS#
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
		