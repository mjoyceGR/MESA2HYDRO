#!/usr/bin/env python
import argparse
import numpy as np
import matplotlib.pyplot as plt

file ="NR_ms_TOL0.01.dat"
#file= "NR_ms_0.98_Rk4_withmass.dat"

N_rk, R_rk, M_rk, E_rk=np.loadtxt(file, usecols=(0,1,2,3), unpack=True)


#file2="NR_ms_quad_from_scipy.dat"
file2="NR_ms_TOL0.05.dat"
N_q, R_q, M_q, E_q=np.loadtxt(file2, usecols=(0,1,2,3), unpack=True)


N_rm, R_rm, M_rm, E_rm=np.loadtxt("NR_ms_Romberg_from_scipy.dat", usecols=(0,1,2,3), unpack=True)


# plt.plot(R_rk, E_rk, "gD", markersize=10, label="RK")
# plt.plot(R_q,  E_q,  "r*", markersize=7, label="quad")
# plt.plot(R_rm,  E_rm,  "b.",markersize=4, label="Romberg")
# plt.legend(loc=1, fontsize=16)
# plt.xlabel("Radial shell ")
# plt.ylabel("Energy per unit particle")
# plt.show()
# plt.close()



mp = 1e-7
Msolar = 1.98e33
N=8
Mshell_target=12.0*N**2.0*mp*Msolar



plt.plot(R_rk, M_rk, "gD", markersize=10, label="RK TOl=0.01")
plt.plot(R_q, M_q, "m*", markersize=8, label="TK TOL=0.05")
plt.xlabel("Radial shell ")
plt.ylabel("Mass contained (?)")
plt.show()
plt.close()

