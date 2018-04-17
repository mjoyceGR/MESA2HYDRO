#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf
#import read_write_HDF5 as rw
import io_lib as rw


MESA_file='../data/profile_whitedwarf_from_mod.data'
#'../data/profile_OB.data'
#'../data/profile_mainsequence.data'

tag=MESA_file.split('../data/profile_')[1].split('.data')[0]

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

# logR=MJ.get_quantity(MESA_file,'logR')
# logRho=MJ.get_quantity(MESA_file,'logRho')
#plt.plot(logR, logRho,'m-', label='text')


fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
plt.plot(fit_region_R, fit_region_rho,'m-', label='text')
plt.xlabel('logR')
plt.ylabel('logRho')
plt.savefig('../tests/'+tag+'_density_profile_'+str(masscut)+'.png')
plt.close()