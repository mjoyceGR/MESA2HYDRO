#!/usr/bin/env python
import sys
MESA_PKG_DIR = os.path.abspath(
    os.path.join(os.path.abspath(os.path.dirname(__file__)), '..'))
sys.path.insert(0, os.path.join(MESA_PKG_DIR, 'lib'))
import converge_funcs as cf

M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar


Ns=[2,4,8,16,32,64]
mps=[1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12]
#1e6,1e5,1e4,1e3,1e2,1,1e-1,1e-2,1e-3,1e-4,
mps_solar=cf.to_array(mps)*M_to_solar

outf=open("N_mp_combinations.dat","w")
print >> outf, 'N 	 	mp        solar_mp 			Mshell(N,mp)'
for i in range(len(mps)):
	print >> outf, "\n"
	for j in range(len(Ns)):
		print >> outf, Ns[j], '\t\t','%1.1e'%mps[i],'\t\t','%1.1e'%mps_solar[i],'\t\t', cf.target_Mshell(Ns[j],mps_solar[i])

outf.close()
