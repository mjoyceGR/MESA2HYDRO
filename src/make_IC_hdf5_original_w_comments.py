def make_IC_box_hdf5(out_fname, mp, healpix_file):
    fname=out_fname#'box_3d_r32.hdf5'; # output filename 

    ###############################################################################################################
    #---------------------------------------------- ASSIGNING -------------------------------------------------
    ###############################################################################################################

    DIMS=2; # number of dimensions 
    N_1D=32; # 1D particle number (so total particle number is N_1D^DIMS)

    Lbox = 1.0 # box side length
    rho_desired = 1.0 # box average initial gas density
    P_desired = 1.0 # initial gas pressure
    #mjoyce- for me, all particles stationary bc star
    vgrainrms=0.0
    #vgrainrms = 1.0 # rms velocity of collisionless particles
    dust_to_gas_ratio = 0.01 # mass ratio of collisionless particles to gas
    gamma_eos = 5./3. # polytropic index of ideal equation of state the run will assume
    
    # first we set up the gas properties (particle type 0)
    
    # make a regular 1D grid for particle locations (with N_1D elements and unit length)
    # x0=np.arange(-0.5,0.5,1./N_1D)
    # x0+=0.5*(0.5-x0[-1]);
    # # now extend that to a full lattice in DIMS dimensions

    # print "x0 grid: ",x0


    # if(DIMS==3):
    #     xv_g, yv_g, zv_g = np.meshgrid(x0,x0,x0, sparse=False, indexing='xy')
    # if(DIMS==2):
    #     xv_g, yv_g = np.meshgrid(x0,x0, sparse=False, indexing='xy')
    #     zv_g = 0.0*xv_g

    # if(DIMS==1):
    #     xv_g=x0
    #     yv_g = 0.0*xv_g #nothing in y
    #     zv_g = 0.0*xv_g #nothing in z

    x,y,z=np.loadtxt(healpix_file, usecols=(0,1,2), unpack=True)

    # the gas particle number is the lattice size: this should be the gas particle number
    #print "xv_g: ", xv_g
    Ngas = len(x)#xv_g.size
    print "Ngas: ",  Ngas
    # flatten the vectors (since our ICs should be in vector, not matrix format): just want a
    #  simple list of the x,y,z positions here. Here we multiply the desired box size in
    # xv_g=xv_g.flatten()*Lbox
    # yv_g=yv_g.flatten()*Lbox 
    # zv_g=zv_g.flatten()*Lbox; 

    x=x*Lbox
    y=y*Lbox 
    z=z*Lbox; 


    # set the initial velocity in x/y/z directions (here zero)
    vx=0.*x 
    vy=0.*y 
    vz=0.*z

    # set the initial magnetic field in x/y/z directions (here zero). 
    #  these can be overridden (if constant field values are desired) by BiniX/Y/Z in the parameterfile
    bx=0.*x
    by=0.*x
    bz=0.*x

    # set the particle masses. Here we set it to be a list the same length, with all the same mass
    #   since their space-density is uniform this gives a uniform density, at the desired value



    ########## mjoyce- assign my particle masses here #############################################
    #mv_g=rho_desired/((1.*Ngas)/(Lbox*Lbox*Lbox)) + 0.*xv_g
    mv_g=mp + 0.*x
    # mjoyce: redundant, see below

    # set the initial internal energy per unit mass. recall gizmo uses this as the initial 'temperature' variable
    #  this can be overridden with the InitGasTemp variable (which takes an actual temperature)
    U=P_desired/((gamma_eos-1.)*rho_desired)

    # uv_g = U + 0.*xv_d#

    # uv_g=np.array(uv_g).astype(float)
    # print "uv_g", uv_g
    # try:
    #     "len(uv_g)", len(uv_g)
    # except:
    #     pass

    # set the gas IDs: here a simple integer list
    id_g=np.arange(1,Ngas+1)
    
    # now we set the properties of the collisionless particles: we will assign these to particle type '3', 
    #   but (barring special compile-time flags being set) GIZMO will treat all collisionless particle types
    #   the same. so the setup would be identical for any of the particle types 1,2,3,4,5

    # set the desired number of particles (here to about twice as many as the gas particles, because we feel like it)
    # Ngrains = np.round(2. * (1.*N_1D)**DIMS)
    # Ngrains=int(Ngrains) #mjoyce

    # set the x/y/z positions: again a simple list for each: here to random numbers from a uniform distribution


    ############## mjoyce READ IN ex) healpix_to_gadget_shell_1.dat HERE ############################
    #print Ngrains, type(Ngrains)
    # xv_d = (np.random.rand(Ngrains)-0.5)*Lbox
    # yv_d = (np.random.rand(Ngrains)-0.5)*Lbox
    # zv_d = (np.random.rand(Ngrains)-0.5)*Lbox


    # print xv_d, '\n', yv_d, '\n', zv_d
    # print "len(xv_d): ",len(xv_d)


    xv_d,yv_d,zv_d=np.loadtxt(healpix_file, usecols=(0,1,2), unpack=True)
    print "len(x_vd) loaded from healpix: ", len(xv_d)

    Ngrains=int(len(xv_d)) #length of arrays loaded from healpix
    print Ngrains, type(Ngrains)
    ###only pick one of Ngrains / Ngas because we only want one particle type


    # set the IDs: these must be unique, so we start from the maximum gas ID and go up
    id_d = np.arange(Ngas+1,Ngrains+Ngas+1)

    # set the velocities. again we will set to a random value, here a Gaussian-distributed one
    #mjoyce: all zeroes
    vx_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vy_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)
    vz_d = np.random.randn(Ngrains) * vgrainrms/np.sqrt(3.)

    # set the masses, again a list with all the same mass
    ############# mjoyce ########################################
    mv_d = mp + 0.*xv_d#
    #mv_d=dust_to_gas_ratio * (1.*Ngas)/(1.*Ngrains) * mv_g[0] + 0.*xv_d
    print "mv_d: ", mv_d, 'len: ', len(mv_d)

    # now we get ready to actually write this out
    #  first - open the hdf5 ics file, with the desired filename
 










    ###############################################################################################################
    #---------------------------------------------- HDF5 WRITING -------------------------------------------------
    ###############################################################################################################
    #
    # PROBABLY DO NOT HAVE TO MODIFY THIS AT ALL
    #

    file = h5py.File(fname,'w') 

    # set particle number of each type into the 'npart' vector
    #  NOTE: this MUST MATCH the actual particle numbers assigned to each type, i.e.
    #   npart = np.array([number_of_PartType0_particles,number_of_PartType1_particles,number_of_PartType2_particles,
    #                     number_of_PartType3_particles,number_of_PartType4_particles,number_of_PartType5_particles])
    #   or else the code simply cannot read the IC file correctly!
    #
    print Ngas, Ngrains
    #npart = np.array([Ngas,0,0,Ngrains,0,0]) # we have gas and particles we will set for type 3 here, zero for all others
    npart = np.array([Ngas,0,0,0,0,0]) # we have gas and particles we will set for type 3 here, zero for all others


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

    ## "create_dataset" acts on p, which creates a group on a file

    # now combine the xyz positions into a matrix with the correct format
    # q=np.zeros((Ngas,3)); q[:,0]=xv_g; q[:,1]=yv_g; q[:,2]=zv_g;
    q=np.zeros((Ngas,3)); q[:,0]=x; q[:,1]=y; q[:,2]=z;
     # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q)
    print "len(coords): ",len(q)



    # similarly, combine the xyz velocities into a matrix with the correct format
    q=np.zeros((Ngas,3)); q[:,0]=vx; q[:,1]=vy; q[:,2]=vz;
    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=q)
    print "len(velocities): ",len(q)

    
    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)
    print "len(id_g): ",len(id_g)

    # write particle masses to the Masses block
    p.create_dataset("Masses",data=mv_g)
    print "len(masses): ",len(mv_g)


    # write internal energies to the InternalEnergy block
    uv_g = U + 0.*xv_d#
    p.create_dataset("InternalEnergy",data=uv_g)
    print "uv_g", uv_g
    try:
        print "len(uv_g)", len(uv_g)
    except:
        pass



    # combine the xyz magnetic fields into a matrix with the correct format

    q=np.zeros((Ngas,3)); q[:,0]=bx; q[:,1]=by; q[:,2]=bz;
    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with 
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can 
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=q)
    print "len(mag field): ",len(q)


    # no PartType1 for this IC
    # no PartType2 for this IC

    # now assign the collisionless particles to PartType3. note that this block looks exactly like 
    #   what we had above for the gas. EXCEPT there are no "InternalEnergy" or "MagneticField" fields (for 
    #   obvious reasons). 

    #mjoyce: NOTE- WE ONLY WANT ONE TYPE OF GAS PARTICLE WITH MASS mp AND WE DID THAT ALREADY

    # p = file.create_group("PartType3")
    # q=np.zeros((Ngrains,3)); q[:,0]=xv_d; q[:,1]=yv_d; q[:,2]=zv_d;
    # p.create_dataset("Coordinates",data=q)
    # q=np.zeros((Ngrains,3)); q[:,0]=vx_d; q[:,1]=vy_d; q[:,2]=vz_d;
    # p.create_dataset("Velocities",data=q)
    # p.create_dataset("ParticleIDs",data=id_d)
    # p.create_dataset("Masses",data=mv_d)

    # no PartType4 for this IC
    # no PartType5 for this IC

    # close the HDF5 file, which saves these outputs
    file.close()
    return file
