#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pygadgetreader as pgr # works- credit this person
import MESAlibjoyce as MJ


mass = MJ.get_quantity(MESA_file,'mass')
logR = MJ.get_quantity(MESA_file,'logR')
logrho = MJ.get_quantity(MESA_file,'logRho')

zone = MJ.get_quantity(MESA_file,'zone')

R = [10.0**p for p in logR]
plt.plot(logR, mass,'g.', label='full profile')
#plt.plot(np.log10(r_region), np.log10(rho_region), 'r.', label=r'outer 5\%')
#plt.xlim(0,100)
plt.xlabel('log radius')
plt.ylabel('mass')
#plt.legend(loc=3, fontsize='small')
#plt.gca().invert_xaxis()
plt.savefig('mass_v_r_MESA.png')
plt.close()


plt.plot(zone, logR, 'm.', label='logR')
plt.plot(zone,mass, 'c.', label='mass')
plt.plot(zone,logrho, 'g.', label='logrho')
plt.xlabel('zone number')
#plt.ylabel('mass')
plt.legend(loc=5, fontsize='small')
plt.savefig('zone3.png')
plt.close()


plt.plot(zone, logR, 'm.', label='logR')
plt.plot(zone,mass, 'c.', label='mass')
#plt.plot(zone,logrho, 'g.', label='logrho')
plt.xlabel('zone number')
#plt.ylabel('mass')
plt.legend(loc=5, fontsize='small')
plt.savefig('zone2.png')
plt.close()