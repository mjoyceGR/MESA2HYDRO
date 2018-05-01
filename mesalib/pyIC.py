import numpy as np
import os
import sys
import struct
# from modules.check import *
# from modules.write import write_header, write_body


print "\n\nWARNING: using pyIC for binary readwrite\n\n"

format_number=1 ##want 2, designed for 1


## restructure and cite Long Do Cao #https://github.com/ldocao
class Header:
    def __init__(self):
        self.NumPart_ThisFile    = np.zeros(6) #number of particles of each type in present file

        self.MassTable           = np.zeros(6).astype(float) #mass of each particle type

        self.Time                = 0 #time of output, or expansion factor for cosmo simu
        self.Redshift            = 0 #redshift
        self.Flag_Sfr            = int(0) #flag star formation (unused in GADGET2)
        self.Flag_Feedback       = int(0) #flag feedback (unused in GADGET2)
        self.NumPart_Total       = np.zeros(6,dtype=np.int8) #total number of particles of each type in simulation
        self.Flag_Cooling        = int(0) #flag for cooling
        self.NumFilesPerSnapshot = int(1) #number of files in each snapshot
        self.BoxSize             = 0 #box size if periodic boundary condition
        self.Omega0              = 0 #matter density at z=0
        self.OmegaLambda         = 0 #vaccum energy at z=0
        self.HubbleParam         = 0 #hubble constant
        self.Flag_StellarAge     = int(0) #creation times of stars (unused)
        self.Flag_Metals         = int(0) #flag mettalicity (unused)

        #self.

        ##not implemented yet
        # self.NumPart_Total_HW    = 0 #not implemented yet. assume number of particles  <2^32
        # self.Flag_Entropy_ICs    = 0


class Body:
    def __init__(self, npart):
        npart=np.array(npart) #make sure it is a numpy array
        total_number_of_particles = np.sum(npart,dtype="int") #was int

        gas_particles=npart[0] #number of gas particles

        #if total_number_of_particles != 0. :
        self.pos = np.zeros([total_number_of_particles,int(3)]) #Positions
        self.vel = np.zeros([total_number_of_particles,int(3)]) #Velocities
        self.id  = (np.zeros(total_number_of_particles)).astype(int)     #Particle ID's

        self.mass =  np.zeros(gas_particles)   #Masses
        ### meridith 
        self.u =  np.zeros(total_number_of_particles) #was gas particles
        self.hsml =(np.zeros(total_number_of_particles)).astype(int)
        self.rho = np.zeros(total_number_of_particles)

        #else:
        #    raise ValueError, "There are no particles !"

        # ##initialize only if enable in makefile to save memory. Maybe there is a smarter way to do that.
        # But the idea is that if you need for example acce, you need all the blocks before. Need to be implemented

        #     self.u =  np.zeros(gas_particles)          #Internal energy per unit mass
        #     self.tstp =np.zeros(total_number_of_particles)
        #     self.endt =np.zeros(gas_particles)
        #     self.acce =np.zeros(total_number_of_particles)
        #     self.pot = np.zeros(total_number_of_particles)            
        #     self.hsml =np.zeros(gas_particles) #SPH Smoothing Length
        #     self.nh = np.zeros(gas_particles) #Hydrogen Abundance
        #     self.ne = np.zeros(gas_particles) #Electron Abundance
        #     self.rho =np.zeros(gas_particles) #Density



def dump_ic(header, body, destination_file="ic.dat", format_output=format_number):#, format_output=1):
    ic_file=open(destination_file,'w')
    write_header(header,ic_file,format_output)
    write_body(body,ic_file,format_output)
    print "=== SUMMARY ==="
    print_summary(header,body)
    ic_file.close()        
    return None





def write_header(header, ic_file, format_output=format_number): ## 1
    """
    Please note that NumPart_Total_HW, nor Flag_Entropy_ICs will be written to the header. Need to be implemented
    """

    print "Writing header (little endian)"
    #* Note that we use struct.pack to form a block, whereas we have to use tostring() on a non-block *#
    #write header into file
    ic_file.write(struct.pack('<I',256))                             #dummy
    ic_file.write(struct.pack('<6I',
                             header.NumPart_ThisFile[0],
                             header.NumPart_ThisFile[1],
                             header.NumPart_ThisFile[2],
                             header.NumPart_ThisFile[3],
                             header.NumPart_ThisFile[4],
                             header.NumPart_ThisFile[5]))
    ic_file.write(struct.pack('<6d',
                             header.MassTable[0],
                             header.MassTable[1],
                             header.MassTable[2],
                             header.MassTable[3],
                             header.MassTable[4],
                             header.MassTable[5]))
    ic_file.write(struct.pack('<d',header.Time))                            #a
    ic_file.write(struct.pack('<d',header.Redshift))                        #z
    ic_file.write(struct.pack('<i',header.Flag_Sfr))                         #sfrFlag
    ic_file.write(struct.pack('<i',header.Flag_Feedback))                          #FBFlag
    ic_file.write(struct.pack('<6I',
                             header.NumPart_Total[0],
                             header.NumPart_Total[1],
                             header.NumPart_Total[2],
                             header.NumPart_Total[3],
                             header.NumPart_Total[4],
                             header.NumPart_Total[5]))
    ic_file.write(struct.pack('<i',header.Flag_Cooling))                     #coolingFlag    
    ic_file.write(struct.pack('<i',header.NumFilesPerSnapshot))                               #numfiles
    ic_file.write(struct.pack('<d',header.BoxSize))                              #boxsize
    ic_file.write(struct.pack('<d',header.Omega0))                              #Omega_0
    ic_file.write(struct.pack('<d',header.OmegaLambda))                              #Omega_Lambda
    ic_file.write(struct.pack('<d',header.HubbleParam))                               #HubbleParam
    ic_file.write(struct.pack('<i',header.Flag_StellarAge))
    ic_file.write(struct.pack('<i',header.Flag_Metals))
    ##should add here NumPart_Total_HW. not implemented yet, nor Flag Entropy


    ##fill in empty space
    header_bytes_left = 260 - ic_file.tell()
    for j in range(header_bytes_left):
        ic_file.write(struct.pack('<x'))
    ic_file.write(struct.pack('<I',256))
    if ic_file.tell()-8 != 256:
        raise IOError, "header has wrong format"
    return None






def write_body(body, ic_file, format_output):
    print "Writing body (little endian)"


    def write_block(block, nbytes, ic_file):
        ic_file.write(struct.pack('<I',nbytes)) #dimensions*number of particles
        ic_file.write(block.tostring())
        ic_file.write(struct.pack('<I',nbytes))
        return None

    total_number_of_particles = np.size(body.pos[:,0])
    gas_particles = np.size(body.u)

    ######## TRY TAKING OUT ALL OF THESE 4'S FOR THE PROFILE DISCREPANCY ISSUE
    ######## doesn't do anything but I'm keeping at at 1 anyway
    ndim=4
    #write in binary format
    write_block(body.pos.astype('f'), 3*ndim*total_number_of_particles, ic_file)

    #print "\n\nvalues in body.pos: ", body.pos
    write_block(body.vel.astype('f'), 3*ndim*total_number_of_particles, ic_file)
    write_block(body.id.astype('I'), ndim*total_number_of_particles, ic_file)
    print "\n\nvalues in body.mass: ", body.mass.astype('f')

    write_block(body.mass.astype('f'), ndim*total_number_of_particles, ic_file)
    write_block(body.u.astype('f'), ndim*gas_particles, ic_file)


    # ##need to set conditions to write these blocks. Maybe it's better to do a loop over each block, but it need some more work.
    # write_block(body.rho.astype('f'), 4*gas_particles, ic_file)
    # write_block(body.ne.astype('f'), 4*gas_particles, ic_file)
    # write_block(body.nh.astype('f'), 4*gas_particles, ic_file)


    ## do not know what I'm doing here
    #print "\n\nvalues in body.rho: ", body.rho

    #### adding this eliminated the loading errors for hsml, u, and rho
    write_block(body.hsml.astype('f'), ndim*total_number_of_particles, ic_file)
    print "\n\nvalues in body.hsml: ", body.hsml.astype('f')    
    write_block(body.u.astype('f'),ndim*total_number_of_particles, ic_file )
    write_block(body.rho.astype('f'),ndim*total_number_of_particles, ic_file )
    
    return None





###############################################
#
# modified version of check.py
#
##############################################
def check_dimension(x,dim):
    if np.shape(x) != dim:
        raise Exception, "Unexpected dimensions"
    return None



def check_header(header):
    """Run a series of tests to check the header

    """
    print "Checking header..."

    ##check if there are some particles
    if (np.sum(header.NumPart_ThisFile) == 0) or (np.sum(header.NumPart_Total) == 0):
        raise ValueError,"No particles in header !"

    ##check if mass are positive
    if np.any(header.MassTable < 0):
        raise ValueError, "MassTable contains negative values"

    ##check if NumFilesPerSnapshot is positive
    if header.NumFilesPerSnapshot <= 0:
        raise ValueError, "NumFilesPerSnapshot is less or equal 0"

    ##check if number of particles have good dimensions
    check_dimension(header.NumPart_ThisFile, (6,))
    check_dimension(header.NumPart_Total, (6,))
    check_dimension(header.MassTable, (6,))
    
    return None



def print_summary(header,body):
    """Print a summary of the parameters
    """
    from pprint import pprint

    ##print header attributes
    #pprint (vars(header))

    ##print body summary
    pprint (vars(body))
    return None