#!/usr/bin/env python
import numpy as np
import os
import sys
import math
import matplotlib.pyplot as plt


# moment of truth 9/22/17
try:
	import healpy as hp
	# import hdf5lib ## 9/5/17
	import read_write_HDF5 as rw
	# import h5py
	import pygadgetreader as pgr
	import MESAlibjoyce as MJ

except:
	print 'Missing module!\nThe following are required: '
	print 'healpy\nhdf5lib\nread_write_HDF5\nh5py\npygadgetreader\nMESAlibjoyce'
	exit(0)



NSIDE = 16

n=1
outf=open('healpix_to_gadget_shell_' + str(n) + '.dat','w')

ipix_array=np.arange(hp.nside2npix(NSIDE)) #this is just a straight up array of particle IDs
x=[]
y=[]
z=[]

for i in range(len(ipix_array)):
	ipix=ipix_array[i]
	coord=hp.pixelfunc.pix2vec(NSIDE, ipix, nest=True)
	print >> outf, coord[0], coord[1], coord[2]
	x.append(coord[0])
	y.append(coord[1])
	z.append(coord[2])


#print x, y, z

##just plots
new_particles=hp.pixelfunc.vec2pix(NSIDE, x, y, z, nest=True)
#hp.mollview(particle_IDs, nest=True, title="Mollview image NESTED")
hp.mollview(new_particles, nest=True,title='')#, title='coordinate form')
plt.savefig('healpix.png')
plt.close()


readfile='/home/meridith/UCT_SAAO/detached_shells/profile99.data'

keyname_list=MJ.get_MESA_output_fields(readfile).keys()
column_dict=MJ.get_columns(readfile,keyname_list)


radius=np.array(column_dict.get('radius'))
logRho=np.array(column_dict.get('logRho'))

# star_age=np.array(column_dict.get('star_age'))
# star_mass=np.array(column_dict.get('star_mass')) 
# log_R=np.array(column_dict.get('log_R'))
# log_L=np.array(column_dict.get('log_L'))
# log_Teff=np.array(column_dict.get('log_Teff')) 

# log_surf_cell_density=np.array(column_dict.get('log_surf_cell_density')) 
# ## is this actually the density-profile density though or do I need 'profile' files for that
# ## no, this will give the evolution of density with age; we need the density profile AT the final age



##########################
#
# 	hdf5 reader
#
#
#(filename, block_name, dim2, parttype=-1, no_mass_replicate=False, fill_block_name="", slab_start=-1, slab_len=-1, verbose=False)
#"RHO ":["Density",1],
#
#########################
hdf5_file='snapshot_010.hdf5'
density=rw.read_block_single_file(hdf5_file,"Density",1)[0][:]
#temp = rw.read_block_single_file(hdf5_file, "Temperature",)[0][:]

x=rw.read_block_single_file(hdf5_file,"Coordinates",3)[0][:][:,0]
y=rw.read_block_single_file(hdf5_file,"Coordinates",3)[0][:][:,1]
z=rw.read_block_single_file(hdf5_file,"Coordinates",3)[0][:][:,2]
print x
print len(x)#,y,z
print len(density)

r= np.sqrt(x**2.0 + y**2.0 + z**2.0)

#print r

plt.plot(r, density,'ro')
#plt.xlim(0,1)
plt.xlabel('radius')
plt.ylabel('density')
plt.savefig('rho_v_r_Gadget_hdf5.png')
plt.close()

#######################################################################################
#
# binary reader
#
#####################################
 # ---------------------
  # -  STANDARD BLOCKS  -
  # ---------------------
  #  pos         - (all)         Position data
  #  vel         - (all)         Velocity data code units
  #  pid         - (all)         Particle ids
  #  mass        - (all)         Particle masses
  #  u           - (gas)         Internal energy
  #  rho         - (gas)         Density
  #  ne          - (gas)         Number density of free electrons
  #  nh          - (gas)         Number density of neutral hydrogen
  #  hsml        - (gas)         Smoothing length of SPH particles
  #  sfr         - (gas)         Star formation rate in Msun/year
  #  age         - (stars)       Formation time of stellar particles
  #  z           - (gas & stars) Metallicty of gas/star particles (returns total Z)
  #  pot         - (all)         Potential of particles (if present in output)

binary_file="/home/meridith/Gadget/Gadget-2.0.7/ICs/AGB_from_glass_binary.dat"

particle_type=0

positions_gas=pgr.readsnap(binary_file,'pos',particle_type)
gas_pos_x=positions_gas[:,0]
gas_pos_y=positions_gas[:,1]
gas_pos_z=positions_gas[:,2]
gas_pos_r=np.power(gas_pos_x*gas_pos_x + gas_pos_y*gas_pos_y + gas_pos_z*gas_pos_z,0.5)  
gas_pos_r=gas_pos_r.astype(np.float)

density=pgr.readsnap(binary_file,'rho',particle_type)

#positions_star=pgr.readsnap(binary_file,'pos',4)
#velocities_gas=pgr.readsnap(binary_file,'vel',0)
#particle_IDs_gas=pgr.readsnap(binary_file,'pid',0)
masses_gas=pgr.readsnap(binary_file,'mass',0)
#mass_star=pgr.readsnap(binary_file, 'mass',4)


plt.plot(radius, logRho,'g-')
#plt.xlim(0,100)
plt.xlabel('radius')
plt.ylabel('log density')
plt.savefig('rho_v_r_MESA.png')
plt.close()


plt.plot(gas_pos_r, masses_gas,'ro')
#plt.xlim(0,1)
plt.xlabel('radius')
plt.ylabel('density')
plt.savefig('rho_v_r_GADGET_binary.png')
plt.close()