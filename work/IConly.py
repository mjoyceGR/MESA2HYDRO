#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/mesalib/')
try:
    import MESA2GADGET.mesalib.MESAlibjoyce as MJ
    import MESA2GADGET.mesalib.converge_funcs as cf
    import MESA2GADGET.mesalib.io_lib as rw
    import MESA2GADGET.mesalib.mainlib as mn
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


make_IC_file=True
# try_reload=True
format_type='binary' 


MESA_file=os.path.join(MESA_PKG_DIR, 'out/sample_MESA_output/profile_mainsequence.data')
masscut=0.95
N=8
mp=1e-7 ##IN UNITES OF Msolar!!!
startype='ms'#'wd_from_mod'

# tag=startype+'_m'+str(masscut)+'_N'+str(N)+'_'+'mp'+str(mp)
# outname=tag

saveNR=os.path.join(MESA_PKG_DIR, "work/NR_files/saveNR_"+startype+".dat")
outname='ms_test'


rough_Nshells=1000.
stepsize=mn.estimate_stepsize(MESA_file,masscut,rough_Nshells)
print "estimated stepsize: ", '%1.5e'%stepsize
if make_NR_file:
	mn.make_NR_file(MESA_file,masscut,N,mp,stepsize,saveNR,check_MESA=check_MESA)
fit_region_R=mn.MESA_r(MESA_file, masscut)
fit_region_rho=mn.MESA_rho(MESA_file, masscut)



rmax=fit_region_R.max()
print "rmax: ", rmax


#ftypes=['hdf5'] #'binary',
#for i in range(len(ftypes)):
#	format_type=ftypes[i]

if make_IC_file:
	mn.get_IC(saveNR,outname,mp,format_type='binary')#temp remove rmax
	mn.get_IC(saveNR,outname,mp,format_type='hdf5')

#if try_reload:
#put GADGET data in simple units

N,r_set,m_cont=np.loadtxt(saveNR, usecols=(0,1,2), unpack=True)

#r_temp, rho_temp=mn.reload_IC(rmax, saveNR,outname,format_type) #rmax instead of 1.0
r_temp_b,masses_b=mn.reload_IC(outname,'binary') #saveNR
r_temp_h,masses_h=mn.reload_IC(outname,'hdf5') #saveNR

print len(masses_b), len(masses_h)
#print 'r_temp_h', r_temp_h
print "r_temp_b/r_temp_h", (r_temp_b/r_temp_h), '\n\n'#/rmax,
# print "masses_b", masses_b
# print "masses_h", masses_h
print "masses_b/masses_h", (masses_b/masses_h)

# plt.plot((r_temp_b-r_temp_h),(masses_b-masses_h),'g.')
# plt.show()
# plt.close()



p_mass=masses_h[5]

r1=r_temp_b.min() 
r_b=[]
rho_b=[]
m_b=[]
num_recovered_b=[]
for i in range(len(r_set)-1):
    r1=r_set[i]
    r2=r_set[i+1]
    region=np.where( (r1<=r_temp_b) &(r2>r_temp_b))

    if len(r_temp_b[region])==0:
        break
    #print "length of r_recovered*mp", len(r_recovered[region])*p_mass
    #print "volume of shell r2-r1: ", cf.volume(r2)-cf.volume(r1)
    #print "r_rec[ r1 to r2 ]*mp / vol(r2-r1): ", len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))
    num_recovered_b.append((len(r_temp_b[region]))) 
    m_b.append(len(r_temp_b[region])*p_mass)
    r_b.append(r2)
    rho_b.append( len(r_temp_b[region])*p_mass/(cf.volume(r2)-cf.volume(r1))  )



r1=r_temp_h.min() 
r_h=[]
rho_h=[]
m_h=[]
num_recovered_h=[]
for i in range(len(r_set)-1):
    r1=r_set[i]
    r2=r_set[i+1]
    region=np.where( (r1<=r_temp_h) &(r2>r_temp_h))

    if len(r_temp_h[region])==0:
        break
    #print "length of r_recovered*mp", len(r_recovered[region])*p_mass
    #print "volume of shell r2-r1: ", cf.volume(r2)-cf.volume(r1)
    #print "r_rec[ r1 to r2 ]*mp / vol(r2-r1): ", len(r_recovered[region])*p_mass/(cf.volume(r2)-cf.volume(r1))
    num_recovered_h.append((len(r_temp_h[region])) )
    m_h.append( (len(r_temp_h[region])*p_mass) ) 
    r_h.append(r2)
    rho_h.append( len(r_temp_h[region])*p_mass/(cf.volume(r2)-cf.volume(r1))  )












plt.plot(r_h, rho_h,'r.', markersize=4, label='GADGET hdf5 data')
plt.plot(r_b, rho_b,'g.', markersize=4, label='GADGET binary data')
plt.plot(fit_region_R, fit_region_rho, "b.", markersize=4, label='MESA data') #cf.to_log()
plt.xlabel("R")
plt.ylabel("test density")
plt.legend(loc=1)
plt.savefig('lin_'+outname+'_both.png')
plt.close()

plt.plot(r_h, num_recovered_h,'r.', markersize=4, label='GADGET hdf5 data')
plt.plot(r_b, num_recovered_b,'g.', markersize=4, label='GADGET binary data')
#plt.plot(fit_region_R, fit_region_, "b.", markersize=4, label='MESA data') #cf.to_log()
plt.xlabel("R")
plt.ylabel("Particle number recovered")
plt.legend(loc=1)
plt.savefig('np_'+outname+'_both.png')
plt.close()



plt.plot(r_h, m_h,'r.', markersize=4, label='GADGET hdf5 data')
plt.plot(r_b, m_b,'g.', markersize=4, label='GADGET binary data')
#plt.plot(fit_region_R, fit_region_, "b.", markersize=4, label='MESA data') #cf.to_log()
plt.xlabel("R")
plt.ylabel("Cummulative Mass")
plt.legend(loc=1)
plt.savefig('mass_'+outname+'_both.png')
plt.close()



plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=4, label='MESA data') #cf.to_log()
plt.plot(r_h, cf.to_log(rho_h),'r.', markersize=4, label='GADGET hdf5 data')
plt.plot(r_b, cf.to_log(rho_b),'g.', markersize=4, label='GADGET binary data')
#plt.ylim(-2.5,0.3)
plt.xlabel("R")
plt.ylabel("log(test density)")
plt.legend(loc=1)
plt.savefig('log_'+outname+'_both.png')
plt.close()

# renormalize the MESA data
#M_to_solar=1.988*10.0**33.0 ## g/Msolar
#fit_region_R=fit_region_R/rmax
#fit_region_rho=fit_region_rho/M_to_solar


# plt.plot(r_temp, rho_temp,'r.', markersize=6, label='GADGET data')
# plt.plot(fit_region_R, fit_region_rho, "b.", markersize=4, label='MESA data') #cf.to_log()
# plt.xlabel("R")
# plt.ylabel("test density")
# plt.legend(loc=1)
# if format_type=='hdf5':
# 	plt.savefig('lin_'+outname+'_hdf5.png')
# else:
# 	plt.savefig('lin_'+outname+'_bin.png')
# plt.close()

# plt.plot(fit_region_R, cf.to_log(fit_region_rho), "b.", markersize=4, label='MESA data') #cf.to_log()
# plt.plot(r_temp, cf.to_log(rho_temp),'r.', markersize=6, label='GADGET data')
# #plt.ylim(-2.5,0.3)
# plt.xlabel("R")
# plt.ylabel("log(test density)")
# plt.legend(loc=1)
# if format_type=='hdf5':
# 	plt.savefig('log_'+outname+'_hdf5.png')
# else:
# 	plt.savefig('log_'+outname+'_bin.png')
# plt.close()


print "total execution length: "
print("--- %s seconds ---" % (time.time() - start_time))
