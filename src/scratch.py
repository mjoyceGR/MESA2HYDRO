#!/usr/bin/env python
import numpy as np
import os
import sys
import math
import matplotlib.pyplot as plt
import hdf5lib ## 9/5/17

import read_write_HDF5 as rw
import h5py
import MESAlibjoyce as MJ

print "testing"

# moment of truth 9/22/17
#try:
import healpy as hp
#print hp

######################################
#
# (1) Operate MESA
#
######################################

MESA_path='/home/meridith/MESA/mesa-r8118/star/work/'


######################################
#
# (2) Extract MESA density profile
#
######################################
in_filename='safecopy.profile99'
out_filename='test_profile_formatted'

readfile=MJ.strip_header(in_filename, out_filename)[1]

#trying to pass output of strip_header function fo these things
keyname_list=MJ.get_fields(readfile).keys()
column_dict=MJ.get_columns(readfile,keyname_list)


radius=np.array(column_dict.get('radius'))  #radius#
r_direct=np.array(column_dict.get('logR'))
#print type(r_direct)#=list(r_direct)
#r_direct=10.0**r_direct  #some clusterfuck going on around here
logRho=np.array(column_dict.get('logRho'))

# star_age=np.array(column_dict.get('star_age'))
# star_mass=np.array(column_dict.get('star_mass')) 
# log_R=np.array(column_dict.get('log_R'))
# log_L=np.array(column_dict.get('log_L'))
# log_Teff=np.array(column_dict.get('log_Teff')) 

# log_surf_cell_density=np.array(column_dict.get('log_surf_cell_density')) 
# ## is this actually the density-profile density though or do I need 'profile' files for that
# ## no, this will give the evolution of density with age; we need the density profile AT the final age

plt.plot(radius, logRho,'g.')
#plt.xlim(0,100)
plt.xlabel('radius log(R/Rstar)')
plt.ylabel('density log_rho g/cm^3')
plt.savefig('rho_v_r_MESA.png')
plt.close()



######################################
#
# (3) Map MESA density profile to a some kind of spherical representation using healpix/healpy
#
######################################

NSIDE = 32
mollview_test = np.arange(hp.nside2npix(NSIDE))
fig1=hp.mollview(mollview_test, title="Mollview image RING")
plt.xlabel('??')
plt.ylabel('mollview_test')
plt.savefig("mollview.png")
plt.close()

fig2=hp.mollview(mollview_test, nest=True, title="Mollview image NESTED")
plt.xlabel('??')
plt.ylabel('nested_mollview')
plt.savefig('nested_mollview.png')
plt.close()

#print type(mollview_test)
#print type(fig1)

######################################
#
# (4) write the Gadget/GIZMO initial conditions file:
# 	see write_GIZMO_ICs.py
#
######################################
out_fname='box_3d_r32.hdf5'
var=MJ.make_IC_box_hdf5(out_fname)
print var
print type(var)

######################################
#
# (5) Run gadget/qsub
#
######################################

# run run_sph.sh 


######################################
#
# (6) Read snapshot/HDF5 file
#
######################################
hdf5_file='snapshot_010.hdf5'

#(filename, block_name, dim2, parttype=-1, no_mass_replicate=False, fill_block_name="", slab_start=-1, slab_len=-1, verbose=False)
#"RHO ":["Density",1],
density=rw.read_block_single_file(hdf5_file,"Density",1)[0][:]  #_single_file
#temp = rw.read_block_single_file(hdf5_file, "Temperature",)[0][:]
x=rw.read_block_single_file(hdf5_file,"Coordinates",3)[0][:][:,0]
y=rw.read_block_single_file(hdf5_file,"Coordinates",3)[0][:][:,1]
z=rw.read_block_single_file(hdf5_file,"Coordinates",3)[0][:][:,2]

#print x
print len(x)#,y,z
print len(density)
r= np.sqrt(x**2.0 + y**2.0 + z**2.0)

plt.plot(r, density,'r.')
#plt.xlim(0,1)
plt.xlabel('radius')
plt.ylabel('density')
plt.savefig('rho_v_r_GADGET.png')
plt.close()


#outf=open('state_of_the_thing.dat','w')
# # for i in range(len(density)):
# # 	print >> outf, density[i] 
# #outf.close()
#plt.clear()




print "done"




