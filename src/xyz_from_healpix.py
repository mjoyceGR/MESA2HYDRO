#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ

NSIDE=16
r_mid=

x,y,z=MJ.get_coords(NSIDE,1)

for i in range(len(x)):
	print x[i]**2.0 + y[i]**2.0  + z[i]**2.0

####################################################
#
# multiply all of these things by the physical radii with which they're associated
# then divide the final ball by outershell's rmax
#
#######################################

# theta0=to_rad(0.0)
# theta=to_rad(90.0)

# print rotate_shell(x[0],y[0],z[0], theta0)
# print rotate_shell(x[0],y[0],z[0], theta)

###############################
#
# WARNING! 	DO NOT ROTATE SHELL PIECEWISE. PICK ONE RANDOM THETA 
# AND MOVE EVERY ENTRY IN THE ARRAY BY THAT SAME VALUE
#
#############################
fixed_random_theta=MJ.random_theta()
# for i in range(len(x)):
# 	print MJ.rotate_shell(x[i],y[i],z[i], fixed_random_theta), "\n"



##just plots
# new_particles=hp.pixelfunc.vec2pix(NSIDE, x, y, z, nest=True)
# #hp.mollview(particle_IDs, nest=True, title="Mollview image NESTED")
# hp.mollview(new_particles, nest=True, title='coordinate form')
# plt.savefig('mollview_nest_w_coords.png')
# plt.close()