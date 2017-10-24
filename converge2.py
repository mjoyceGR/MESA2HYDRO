#!/usr/bin/env python
import converge_funcs as cf

Rsolar = 6.955e10 #[cgs]
Msolar= 2.0e30 #[kg]
Msolar_cgs= 2.0e33 #[cgs]


fname='profile175.data'#'profile175.data' #175, 140, 32
MESA_file="{}".format(fname)

n_p_initial=3000 #72#3072#100,000 
stepsize=0.001#

r_array = cf.get_MESA_profile_edge(MESA_file, quantity='logR', strip=False)
M_array = cf.get_MESA_profile_edge(MESA_file, quantity='mass', strip=False)

cf.do_converge(MESA_file, r_array, M_array, n_p_initial,stepsize,\
outputfile='test_vals_after_overhaul.dat', mp=1e-5)#0.000102319572081 )#)0.1)
## making it pick "mp" with my estimator is definitely breaking this
