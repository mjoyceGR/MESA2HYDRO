#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import random as rand

try:
	import healpy as hp
	#import MESAlibjoyce as MJ

except:
	print 'Missing module!\nThe following are required: '
	print 'healpy\n'
	exit(0)


def get_coords(NSIDE,file_index):
	n=file_index
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
	outf.close()
	print 'file healpix_to_gadget_shell_', str(n), '.dat generated'
	return x, y, z


def rotate_shell(single_x,single_y,single_z, theta):
	vec=np.matrix([ [single_x], [single_y], [single_z]])

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

	xnew=Rx*vec
	ynew=Ry*vec
	znew=Rz*vec
	return xnew,ynew,znew

def to_rad(theta):
	theta=theta*np.pi/180.0
	return theta

def random_theta():
	theta=rand.random()*2.0*np.pi #.random gives random float between 0 and 1
	#theta=rand.randrange(0.0, 2*np.pi, 0.01)
	return theta

NSIDE=16
x,y,z=get_coords(NSIDE,1)

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
fixed_random_theta=random_theta()
for i in range(len(x)):
	print rotate_shell(x[i],y[i],z[i], fixed_random_theta), "\n"
#print x, y, z

##just plots
# new_particles=hp.pixelfunc.vec2pix(NSIDE, x, y, z, nest=True)
# #hp.mollview(particle_IDs, nest=True, title="Mollview image NESTED")
# hp.mollview(new_particles, nest=True, title='coordinate form')
# plt.savefig('mollview_nest_w_coords.png')
# plt.close()