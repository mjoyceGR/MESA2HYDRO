#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf
#import read_write_HDF5 as rw
import io_lib as rw

unlog=True
#file paper was written on:
MESA_file='../data/profile140.data'

# MESA_file='profile_mainsequence.data'
#'../out/profile_AGB_from_mod.data'
#MESA_file='../out/profile_whitedwarf_from_mod.data'
#MESA_file='../out/profile_OB.data'
#MESA_file='../out/profile_mainsequence.data'

#tag=MESA_file.split('../out/profile_')[1].split('.data')[0]
#tag='profile140'
tag='p140'

MJ.show_allowed_MESA_keywords(MESA_file)
#MESA_file='../data/profile140.data' #140, 32
#####################################################
#
# 	allowed keywords:
#
# logT
# logR
# zone
# logP
# z_mass_fraction_metals
# logRho
# mass
# x_mass_fraction_H 
# y_mass_fraction_He
#
####################################################
masscut=0.95#0.95

fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_M = cf.get_MESA_profile_edge(MESA_file, quantity='mass', masscut=masscut ,strip=False)
if unlog:
	fit_region_R=cf.unlog(fit_region_R)
	fit_region_M=cf.unlog(fit_region_M)
plt.plot(fit_region_R, fit_region_M,'k-', color='darkblue')#, label='text')
if unlog:
	plt.xlabel('R (unlog)')
	plt.ylabel('M (unlog)')
else:
	plt.xlabel('logR')
	plt.ylabel('mass')
plt.savefig(tag+'_mass_profile_'+str(masscut)+'.png')
plt.close()