#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import MESAlibjoyce as MJ
import converge_funcs as cf

# fname='../data/profile_mainsequence.data'
which='whitedwarf_from_mod'#'OB','mainsequence'


histfile='../data/history_'+which+'.data'
tag=histfile.split('../data/history_')[1].split('.data')[0]

MJ.show_allowed_MESA_keywords(histfile)

logL=MJ.get_quantity(histfile,'log_L')
logT=MJ.get_quantity(histfile, 'log_Teff')
# L=cf.unlog(logL)
# T=cf.unlog(logT)

print logL
print logT

ax, fig=plt.subplots()
plt.plot(logT, logL,'g-')
plt.xlabel('log T')
plt.ylabel('Log L')
plt.gca().invert_xaxis()
plt.savefig('HRtest_'+str(tag)+'.png')
plt.close()

