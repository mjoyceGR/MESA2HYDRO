#!/usr/bin/env python
import numpy as np
import struct
import matplotlib.pyplot as plt
import pygadgetreader as pgr # works- credit this person

binary_file="AGB_from_glass_binary.dat"


def figure_out_mp_given_gadget_data(binary_file,rl,ru):
	positions_gas=pgr.readsnap(binary_file,'pos',0)
	gas_pos_x=positions_gas[:,0]
	gas_pos_y=positions_gas[:,1]
	gas_pos_z=positions_gas[:,2]

	gas_pos_r=np.power(gas_pos_x*gas_pos_x + gas_pos_y*gas_pos_y + gas_pos_z*gas_pos_z,0.5)  
	gas_pos_r=gas_pos_r.astype(np.float)
	masses_gas=pgr.readsnap(binary_file,'mass',0)

	region=np.where(  (gas_pos_r>rl) & (gas_pos_r<ru) )[0]
	pos_slice=gas_pos_r[region]
	mass_slice=masses_gas[region]

	## possibly replace this with a free parameter
	mp=np.average(mass_slice)
	n_p=len(mass_slice)
	Mshell=n_p*mp

	n1 = np.sqrt(Mshell/(12.0*mp))
	n2 = calc_n2(rl, ru)

	return n1, n2, mp


def Mshell(rl, ru, A,B,C):
	rdiff=(ru - rl)
	f1=0.5*A*rdiff**2.0  + A*B*rdiff  + A*B**2.0 * np.log( abs(rdiff - B) ) + (1.0/3.0)*C*rdiff**3.0
	return 4.0*np.pi*f1


def calc_n2(rl, ru):
	return np.sqrt(np.pi/12.0)*(ru + rl)/(ru-rl)


def get_n1_n2(rl,ru, A, B, C):
	M_shell=Mshell(rl, ru, A, B, C)
	n_p=M_shell/mp

	n1 = np.sqrt(M_shell/(12.0*mp))
	n2 = calc_n2(rl, ru)

	return n1, n2

mp= 5000 #8.266e-12 #fix this to whatever

## these were found by numirically fitting the last 5% of the MESA profile
A= 25795.790860453468
B=-2.9056248746991034
C= 17446.663602504697

rl=3.1 #the starting guess for the lower radius comes from the MESA output


#rl=0.5
# should be tuned according to how radial thing looks

incr=0.0000001
begin=rl+incr
end = 5.5#50*begin

weak_tolerance=10.0
strong_tolerance=0.5

ru=np.arange(begin,end,incr)

rl_list=[]
ru_list=[]
n1_list=[]
n2_list=[]


rl_list_conv=[]
ru_list_conv=[]
n1_list_conv=[]
n2_list_conv=[]


for i in range(len(ru)):
	quant=get_n1_n2(rl,ru[i], A, B, C)
	n1 = quant[0]
	n2 = quant[1]

	# rl_list.append(rl)
	# ru_list.append(ru[i])
	# n1_list.append(n1)
	# n2_list.append(n2)

	#print "n1: ", n1, " n2: ", n2, "\t at ru=", ru[i]
	if abs(n2 -n1) <= weak_tolerance:
		#print "n1: ", n1, " n2: ", n2, "\t at ru=", ru[i]

		if abs(n1 - n2)<=strong_tolerance:
			print "\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru[i]

			# R_shell_vals.append(rl)
			# ru_list.append(ru[i])

			rl=ru[i]
		else:
	 		pass
	 		
i=0
for i in range(len(ru_list)-1):
	print "ru: ", ru_list[i],  "\t rl: ",R_shell_vals[i], "diff: ", ru_list[i]-R_shell_vals[i]



 #from fit range to tail of MESA profile; CAN'T BE THE SAME AS ru[0] or log will freak out

# def approximate_mp(N,ru, rl, A,B,C):
# 	M_shell=Mshell(rl,ru,A,B,C)
# 	mp=Mshell/N**2.0
# 	return mp
# test_N=10.0e7
# test_ru=3.2
# mp = approximate_mp(test_N,test_ru,rl,A,B,C)
# print "value of mp: ", mp
# ru=np.arange(3.0,5.5,0.01)
#print ru


# plt.plot(ru, density_integral(ru,rl,A,B,C), 'g-', label='LHS(ru) (density integral)')
# plt.plot(ru, RHS(ru, rl, mp), 'm-', label='RHS(ru)')
# plt.legend(loc=1)#, fontsize='x-small')
# plt.xlabel("radius ru")
# plt.ylabel("???")
# plt.savefig('find_intersection.png')
# plt.close()