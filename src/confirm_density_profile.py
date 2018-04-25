#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf
#import read_write_HDF5 as rw
import io_lib as rw

###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################

unlogg_data=True
#file paper was written on:
MESA_file='../out/raw_MESA_output/profile_agb_timmes.data'
#profile_mainsequence.data'
#MESA_file2='../data/profile140.data'

#tag='profile140'
tag='agb_timmes'

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
masscut=0.95

fit_region_R = cf.get_MESA_profile_edge(MESA_file, quantity='logR', masscut=masscut,strip=False)
fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='logRho', masscut=masscut ,strip=False)
#fit_region_rho=cf.to_log(fit_region_rho)

if unlogg_data:
	fit_region_R=cf.unlog(fit_region_R)
	fit_region_rho=cf.unlog(fit_region_rho)

# fit_region_rho2 = cf.get_MESA_profile_edge(MESA_file2, quantity='logRho', masscut=masscut ,strip=False)
# fit_region_R2 = cf.get_MESA_profile_edge(MESA_file2, quantity='logR', masscut=masscut,strip=False)

# if unlog:
# 	fit_region_R=cf.unlog(fit_region_R)
# 	fit_region_rho = cf.get_MESA_profile_edge(MESA_file, quantity='rho', masscut=masscut ,strip=False)
# 	# fit_region_R2=cf.unlog(fit_region_R2)
# 	# fit_region_rho2=cf.unlog(fit_region_rho2)

plt.plot(fit_region_R*R_to_solar, fit_region_rho,'m-', label=tag)
#plt.plot(fit_region_R2, fit_region_rho2,'b-', label='p140')
if unlogg_data:
	plt.xlabel('R (unlog)')
	plt.ylabel('Rho (unlog)')
else:
	plt.xlabel('logR')
	plt.ylabel('logRho')
plt.legend(loc=1)
plt.savefig(tag+'_density_profile_'+str(masscut)+'.png')
plt.close()