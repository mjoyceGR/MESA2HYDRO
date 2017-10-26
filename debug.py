#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf
import read_write_HDF5 as rw

###################### debug ###########################################3
m2='history.data'

#print MJ.show_allowed_MESA_keywords(m2)
#exit(0)
hist_L=MJ.get_quantity(m2,"log_L")
hist_r=MJ.get_quantity(m2, "log_R")
hist_teff=MJ.get_quantity(m2,"log_Teff")
age=MJ.get_quantity(m2,'star_age')

region=np.where( 10e5 < hist_L) 

model_num=MJ.get_quantity(m2,'model_number')
print model_num.astype(int)
#his_r=cf.unlog(history_r)/Rsolar
#hist_L=cf.to_log(hist_L)#/Lsolar

#plt.plot(hist_teff[region], hist_L[region], 'r-')
#plt.plot(hist_teff, hist_L, 'r-')
plt.plot(hist_r, hist_L, 'r-')
plt.xlabel('logR')
plt.ylabel('luminostiy')
plt.savefig('whatisABGradius.png')
plt.close()

print "plot made"

