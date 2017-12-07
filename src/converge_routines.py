#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pygadgetreader as pgr # works- credit this person
import MESAlibjoyce as MJ


def get_MESA_profile_edge(MESA_file, strip):
	#print MJ.show_allowed_MESA_keywords(MESA_file)
	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]

	mass = MJ.get_quantity(MESA_file,'mass')
	logR = MJ.get_quantity(MESA_file,'logR')
	logrho = MJ.get_quantity(MESA_file,'logRho')
	logrho = logrho.astype(np.float)
	rho = [10.0**p for p in logrho]

	mass = mass.astype(np.float)
	Mtot=mass[-1]
	# print "last mass entry (total M??): ", Mtot, "\tlargest mass value: ", max(mass),\
	#  "\tmass value at most inner R: ", mass[np.argmin(logR)]
	### warning: argmin returns first occurance only if given a set of co-miminal (maximal) values
	cummulative=[]
	mass_intergral=0
	for i in range(len(mass)):
		mass_intergral = mass_intergral + mass[i]
		cummulative.append(mass_intergral)
	cummulative= np.array(cummulative)
	cummulative=cummulative.astype(np.float)

	fit_region_rho=[]
	fit_region_R=[]
	fit_region_M=[]
	for m in range(len(cummulative)):
		while cummulative[m] <= 0.05*cummulative[np.argmin(logR)]:
			fit_region_rho.append( (10.0**float(logrho[m])))#/100000.0 )) 
			fit_region_R.append(10.0**float(logR[m]))
			fit_region_M.append(float(cummulative[m]))
			m+=1 
	fit_region_R   = np.array(fit_region_R).astype(np.float)
	fit_region_rho = np.array(fit_region_rho).astype(np.float)
	fit_region_M   = np.array(fit_region_M).astype(np.float)
	
	return fit_region_R, fit_region_rho, fit_region_M 


def get_first_mp(MESA_file, n_p_initial):#(r_array, M_array, rl_init, n_p):
	try:
		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, False)
	except:
		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, True) # not sure if this is the right order here

	rl=r_array.min() #3.1
	rl_init=rl

	pos_r=r_array.astype(np.float)  
	Mcumm=M_array.astype(np.float)
	#rl = rl_init

	step_size=0.0000001
	region=[]
	while len(region) ==0: 
		ru = rl_init + step_size
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		pos_slice=pos_r[region]

		step_size=step_size*5.0
		print "trying step size: ", step_size
	else:
		ru = rl_init + step_size
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		pos_slice=pos_r[region]
		#print "while loop broken"
		#break

	print "step size selected: ", step_size
	Mshell=Mcumm[region.argmax()] - Mcumm[region.argmin()] #this should be a single, float value
	n_p=n_p_initial
	mp = Mshell/n_p

	print "selected mass per particle mp=",mp

	return mp



def n1_n2_MESA(r_array, M_array, rl, ru, mp, step_size):
	pos_r=r_array.astype(np.float)  
	Mcumm=M_array.astype(np.float)

	region=[]
	while len(region) ==0: 
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		ru= ru + 0.1*stepsize #magic numbering
		#print "array length: ", len(region)
		#print "step_size in n1_n2_MESA: ", step_size, " ru: ", ru
	else:
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		pos_slice=pos_r[region]
		mass_slice=Mcumm[region]

	Mshell=mass_slice[mass_slice.argmax()] - mass_slice[mass_slice.argmin()]

	n_p=Mshell/mp

	n1 = np.sqrt(Mshell/(12.0*mp))
	n2 = calc_n2(rl, ru)
	return n1, n2, ru, n_p


def calc_n2(rl, ru):
	return np.sqrt(np.pi/12.0)*(ru + rl)/(ru - rl)



######################### in below, assign returned ru value to rl, run again ##############################
def get_N(r_array, M_array, rl, ru, mp, stepsize):#(MESA_file, mp, step_size):
	n1=0.0
	n2=100.0
	weak_tolerance=10.0
	strong_tolerance=0.1
	temp=9999

	while abs(n2 -n1) > 0.0:#weak_tolerance:
		quant=n1_n2_MESA(r_array, M_array, rl, ru, mp, stepsize)
		n1 = quant[0]
		n2 = quant[1]
		ru = ru + stepsize
		n_p = quant[3]

		#if abs(n2 -n1) <= weak_tolerance:
			#print "n1: ", n1, " n2: ", n2, "\t at ru=", ru
		if abs(n1 - n2)<=strong_tolerance:
			#while abs(n2 - n1) <= strong_tolerance:
			print "\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru
			temp = ru					
			#else:
	 		#	break#pass

		if ru > temp:#+3.0*stepsize:
			break
	set_R= (temp + rl)/2.0	

	if np.floor(n1) == np.floor(n2):
		return rl, temp, set_R, np.floor(n1), n_p #, mp
	else:
		print "FAILED TO CONVERGE"
		return 0


# def strip_MESA_header(in_filename, out_filename, *args, **kwargs):
#     ### add hook to prevent double removes
# 	#[0]-- returns the file object
# 	#[1]-- returns the string/name of the reformatted text file
# 	n = int(kwargs.get('n', 5))


def do_converge(MESA_file,n_p_initial,stepsize,*args,**kwargs):
	filename=str(kwargs.get('outputfile','vals.dat'))
	outf=open(filename,"a")
	try:
		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, False)
	except:
		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, True) # not sure if this is the right order here

	rl=r_array.min() #3.1
	rmax = r_array.max() #5.5#50*begin

	#particles per shell; mp is determined based on this
	first_mp=get_first_mp(MESA_file, n_p_initial)
	print "\nfirst mp: ", first_mp, " first np: ", n_p_initial
	mp=first_mp
	ru = rl + stepsize


	print >> outf, 'MESA_file: ', MESA_file
	print >> outf, 'rl				ru				N		n_p				mp			stepsize'
	print >> outf, rl,'\t',ru, "\tN", "\t", n_p_initial, "\t", mp, "\t", stepsize

	while ru < rmax:
		new_vals= get_N(r_array, M_array, ru, ru+stepsize, mp, stepsize) ## not sure about this ru, rl assignment
		try:
			rl=new_vals[0]
		except TypeError:
			print '\n\nRoutine terminated\n\n'
			break	
		ru=new_vals[1]
		rmid=new_vals[2]
		N=new_vals[3]
		n_p=new_vals[4]

		print >> outf, rl,'\t',ru, "\t", N, "\t", n_p, "\t", mp, "\t", stepsize
		print "now: N=", N, " at r_mid=", rmid, ' rl=', rl, ' ru=', ru, 'n_p=',n_p, " and mp=",mp, "\n\n"

	outf.close()
	print 'profile data generated in file', filename
	return 0  


fname='profile175.data'#'profile175.data' #140, 32
MESA_file="{}".format(fname)
n_p_initial=5000#100,000 
stepsize=0.0001

do_converge(MESA_file, n_p_initial,stepsize,outputfile='test_vals.dat')