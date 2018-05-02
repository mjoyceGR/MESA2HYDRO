#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/lib/')
try:
    import MESA2GADGET.lib.MESAlibjoyce as MJ
    import MESA2GADGET.lib.converge_funcs as cf
    import MESA2GADGET.lib.io_lib as rw
    import MESA2GADGET.lib.mainlib as mn
    from MESA2GADGET import MESA_PKG_DIR
except ImportError:
    print("Problem with MESA2GADGET installation")
    print("To use this package please run sudo python setup.py install")
    print("or set your PYTHONPATH environment variable to the directory")
    print("MESA2GADGET is in (pointing it directly to MESA2GADGET still causes problems)")
    exit(1)

import time
start_time = time.time()

check_MESA=False
make_NR_file=False
make_IC_file=False


MESA_file=os.path.join(MESA_PKG_DIR, 'out/sample_MESA_output/profile_mainsequence_logE.data')
masscut=0.95
N=8
mp=1e-7 ##IN UNITES OF Msolar!!!
startype='ms'#'wd_from_mod'

saveNR='/home/meridith/UCT_SAAO/detached_shells/MESA2GADGET/work/NR_files/saveNR_ms_logE.dat'
outname='ms_logE_test'


#rough_Nshells=1000.
#stepsize=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
#print "estimated stepsize: ", '%1.5e'%stepsize
if make_NR_file:
	mn.make_NR_file(MESA_file,masscut,N,mp,stepsize,saveNR,check_MESA=check_MESA)

fit_region_R=mn.MESA_r(MESA_file, masscut)
fit_region_rho=mn.MESA_rho(MESA_file, masscut)
rmax=fit_region_R.max()



if make_IC_file:
	mn.get_IC(saveNR,outname,mp,format_type='binary')#temp remove rmax
	mn.get_IC(saveNR,outname,mp,format_type='hdf5')


#r_temp, rho_temp=mn.reload_IC(rmax, saveNR,outname,format_type) #rmax instead of 1.0
r_temp_b,masses_b=mn.reload_IC(outname,'binary') #saveNR
r_temp_h,masses_h=mn.reload_IC(outname,'hdf5') #saveNR

print "recovered array lengths: ", len(masses_b), len(masses_h)
print 'r_temp_h', len(r_temp_h), len(r_temp_b)
print "r_temp_b/r_temp_h", (r_temp_b/r_temp_h), '\n\n'#/rmax,
print "masses_b/masses_h", (masses_b/masses_h)

p_mass=masses_h[5]


nbin=70.
r_b,rho_b=mn.binned_r_rho(r_temp_b, nbin,p_mass)
r_h,rho_h=mn.binned_r_rho(r_temp_h, nbin,p_mass)

print "number of points num_recovered", len(r_b), len(rho_b)



rho_h=2.0*cf.to_array(rho_h)/np.pi


plt.plot(r_h, rho_h,'r^', markersize=5, label='GADGET hdf5 data')
plt.plot(r_b, rho_b,'g.', markersize=4, label='GADGET binary data')
plt.plot(fit_region_R, fit_region_rho, "b.", markersize=4, label='MESA data') #cf.to_log()
plt.xlabel("R")
plt.ylabel("test density")
plt.legend(loc=1)
plt.savefig('lin_'+outname+'_both.png')
plt.close()

plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=4, label='MESA data') #cf.to_log()
plt.plot(r_h, cf.to_log(rho_h),'r^', markersize=5, label='GADGET hdf5 data')
plt.plot(r_b, cf.to_log(rho_b),'g.', markersize=4, label='GADGET binary data')
plt.ylim(-2.5,0.3)
plt.xlabel("R")
plt.ylabel("log(test density)")
plt.legend(loc=1)
plt.savefig('log_'+outname+'_both.png')
plt.close()

print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))











# plt.plot(r_h, num_recovered_h,'r^', markersize=5, label='GADGET hdf5 data')
# plt.plot(r_b, num_recovered_b,'g.', markersize=4, label='GADGET binary data')
# #plt.plot(fit_region_R, fit_region_, "b.", markersize=4, label='MESA data') #cf.to_log()
# plt.xlabel("R")
# plt.ylabel("Particle number recovered")
# plt.legend(loc=1)
# plt.savefig('np_'+outname+'_both.png')
# plt.close()

# plt.plot(r_h, m_h,'r^', markersize=5, label='GADGET hdf5 data')
# plt.plot(r_b, m_b,'g.', markersize=4, label='GADGET binary data')
# #plt.plot(fit_region_R, fit_region_, "b.", markersize=4, label='MESA data') #cf.to_log()
# plt.xlabel("R")
# plt.ylabel("Cummulative Mass")
# plt.legend(loc=1)
# plt.savefig('mass_'+outname+'_both.png')
# plt.close()
