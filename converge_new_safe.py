#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import pygadgetreader as pgr # works- credit this person
import MESAlibjoyce as MJ


def get_MESA_profile_edge(MESA_file,**kwargs):#, strip):
	#print MJ.show_allowed_MESA_keywords(MESA_file)
	strip=bool(kwargs.get('strip',False))
	keyword=str(kwargs.get('quantity','zone'))

	print "keyword: ",keyword


	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]
	
	try:
		quantity=MJ.get_quantity(MESA_file,keyword)
	except:
		print "Quantity keyword not found. Allowed keywords are:"
		print MJ.show_allowed_MESA_keywords(MESA_file)

	masses = MJ.get_quantity(MESA_file,'mass')
	masses = masses.astype(np.float)
	Mtot=masses[0]

	## THE "mass" VARIABLE ALREADY DESCRIBES THE MASS CONTAINED
	cummulative=masses#mass.sort()

	fit_region_indices=[]
	for m in range(len(cummulative)):
		#FOR EVERYTHING THAT'S NOT MASS, USE THE OUTER REGIONS:
		if cummulative[m] >= 0.95*Mtot:
			fit_region_indices.append(m)	
	fit_region  = quantity[fit_region_indices]

	print "keyword: ",keyword
	if keyword =='mass':
		#xq_unlogged = [10.0**p for p in xq]
		mf = fit_region
		fit_region = [Mtot-p for p in mf]

	fit_region = np.array(fit_region).astype(float)
		
	return fit_region


def get_first_mp(MESA_file, n_p_initial):#(r_array, M_array, rl_init, n_p):
	#try:
	r_array = get_MESA_profile_edge(MESA_file, quantity='logR', strip=False)
	rho_array = get_MESA_profile_edge(MESA_file, quantity='logRho', strip=False)
	M_array = get_MESA_profile_edge(MESA_file, quantity='mass', strip=False)
	# except: #what kind of error happpens if the file isn't stripped? and when?
	# 	r_array = get_MESA_profile_edge(MESA_file, quantity='logR', strip=True)
	# 	rho_array = get_MESA_profile_edge(MESA_file, quantity='logRho', strip=True)
	# 	M_array = get_MESA_profile_edge(MESA_file, quantity='mass', strip=True)

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
	## WARNING!!!!!!!!!!!!!! THIS MIGHT BE BACKWARDS NOW
	Mshell=abs(Mcumm[region.argmin()] - Mcumm[region.argmax()]) #this should be a single, float value
	n_p=n_p_initial
	mp = Mshell/n_p

	print "selected mass per particle mp=",mp

	return mp

# def M_contained(MESA_file, ru):
# 	cumm_M_array=get_MESA_profile_edge(MESA_file,quantity='mass')
# 	r_array=get_MESA_profile_edge(MESA_file,quantity='logR')
# 	idx = where r_array equals ru, then pass that index to cumm_M_arry and print MCA[idx]
# 	#print "cummulative M array: ", cumm_M_array
# 	Mc=MCA[idx] #sum(cumm_M_array)
# 	return Mc


def n1_n2_MESA(r_array, M_array, rl, ru, mp, step_size):
	pos_r=r_array.astype(np.float)  
	Mcumm=M_array.astype(np.float)

	region=[]
	ru=rl
	while len(region) ==0: 
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		ru= ru + 0.00001#0.1*stepsize #magic numbering
		#print "array length: ", len(region)
		#print "step_size in n1_n2_MESA: ", step_size, " ru: ", ru
	#else:
	region=np.where( (pos_r>rl) & (pos_r<temp_ru) )[0]
	pos_slice=pos_r[region]
	mass_slice=Mcumm[region]

	#print "\n\nmass_slice:", mass_slice, 
	print "\nargmin: ", mass_slice.argmin(),"\tslice[argmin]: ",mass_slice[mass_slice.argmin()]
	print 'mass slice.min()', mass_slice.min(),'mass slice.max()', mass_slice.max()

	########  AMBIGUOUS MASS ORIENTATION, MAY NEEED TO REVERSE THISSSSSSSSSSSSS
	Mshell=abs(mass_slice[mass_slice.argmin()] - mass_slice[mass_slice.argmax()])
	n_p=Mshell/mp
	print "Mshell: " ,Mshell, "n_p", n_p


	n1 = np.sqrt(Mshell/(12.0*mp))
	n2 = calc_n2(rl, ru)
	return n1, n2, ru, n_p

#n1_n2_MESA()


def calc_n2(rl, ru):
	return np.sqrt(np.pi/12.0)*(ru + rl)/(ru - rl)



######################### in below, assign returned ru value to rl, run again ##############################
def get_N(r_array, M_array, rl, ru, mp, stepsize, n_p_initial,tolerance,**kwargs):#(MESA_file, mp, step_size):

	strong_tolerance=tolerance#float(kwargs.get('tolerance',0.1))

	n1=0.0
	n2=100.0
	weak_tolerance=10.0
	#strong_tolerance=0.1
	keep_ru=9999

	while abs(n2 -n1) > 0.0:#weak_tolerance:
		#print #'r_array: ',r_array, 
		#print 'M_array', M_array, '\nrl: ',rl, '\tru: ',ru, '\n'
		quant=n1_n2_MESA(r_array, M_array, rl, ru, mp, stepsize)
		#print 'quant: ', quant, '\n\n'
		### R is fine, mass is not

		n1 = quant[0]
		n2 = quant[1]
		ru = ru + stepsize
		print "ru now: ", ru
		n_p = quant[3]

#		if debug:
		#if abs(n2 -n1) <= weak_tolerance:
		print "n1: ", n1, " n2: ", n2, "\t at ru=", ru

		if abs(n1 - n2)<=strong_tolerance:
			#while abs(n2 - n1) <= strong_tolerance:
			print "tolerance met\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru
			keep_ru = ru
		# else:
		# 	#CAREFUL TWO OF THESE
		# 	keep_ru = 0#ru + stepsize					
		# 	#else:
	 # 		#	break#pass

	 	#if ru > r_array.max():
	 	#	break	
	
		#if ru > temp:#+3.0*stepsize:
		#	break

	set_R= (keep_ru + rl)/2.0
	#temp = 'error'	

	#if np.floor(n1) == np.floor(n2):
	return rl, keep_ru, set_R, np.floor(n1), n_p #, mp
	# else:
	# 	print "\n\nFAILED TO CONVERGE with tolerance ", strong_tolerance,\
	# 	 ', initial particle number ', n_p_initial,\
	# 	 ', step size ', stepsize
	# 	return 0


def do_converge(MESA_file,n_p_initial,stepsize,*args,**kwargs):
	filename=str(kwargs.get('outputfile','vals.dat'))
	tolerance=float(kwargs.get('tolerance',0.1))

	outf=open(filename,"a")
	#try:
	r_array = get_MESA_profile_edge(MESA_file, quantity='logR', strip=False)
	rho_array = get_MESA_profile_edge(MESA_file, quantity='logRho', strip=False)
	M_array = get_MESA_profile_edge(MESA_file, quantity='mass', strip=False)
	# except: #what kind of error happpens if the file isn't stripped? and when?
	# 	r_array = get_MESA_profile_edge(MESA_file, quantity='logR', strip=True)
	# 	rho_array = get_MESA_profile_edge(MESA_file, quantity='logRho', strip=True)
	# 	M_array = get_MESA_profile_edge(MESA_file, quantity='mass', strip=True)

	rl=r_array.min() #3.1
	rmax = r_array.max() #5.5#50*begin

	#particles per shell; mp is determined based on this

	#first_mp=get_first_mp(MESA_file, n_p_initial)
	first_mp=10e-12 

	print "\nfirst mp: ", first_mp, " first np: ", n_p_initial
	mp=first_mp
	ru = rl + stepsize

	print >> outf, 'MESA_file: ', MESA_file
	print >> outf, 'rl				ru				N		n_p				mp			stepsize'
	print >> outf, rl,'\t',ru, "\tN", "\t", n_p_initial, "\t", mp, "\t", stepsize

	while ru < rmax:
		new_vals= get_N(r_array, M_array, ru, ru+stepsize, mp, stepsize, n_p_initial, tolerance) 
		## not sure about this ru, rl assignment
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
n_p_initial=100#3072#100,000 
stepsize=0.1#1e-6#0.01

do_converge(MESA_file, n_p_initial,stepsize,outputfile='test_vals_after_overhaul.dat', tolerance=0.1)#)0.1)






#------------------------------- plotting routines ----------------------------------
def plot_region(MESA_file,**kwargs):
	# filename=str(kwargs.get('outputfile','vals.dat'))
	# tolerance=float(kwargs.get('tolerance',0.1))

	x_quantity=str(kwargs.get('x','logR'))
	y_quantity=str(kwargs.get('y','logRho'))
	unlog=bool(kwargs.get('unlog',False))


	xq = MJ.get_quantity(MESA_file,x_quantity)
	yq = MJ.get_quantity(MESA_file,y_quantity)

	x_region= get_MESA_profile_edge(MESA_file, quantity=x_quantity, strip=False)#[0]
	y_region = get_MESA_profile_edge(MESA_file, quantity=y_quantity , strip=False)#[1]

	if unlog:
		xq_unlogged = [10.0**p for p in xq]
		yq_unlogged = [10.0**p for p in yq]
		x_region_un = [10.0**p for p in x_region]
		y_region_un = [10.0**p for p in y_region]
		plt.plot(xq_unlogged, yq_unlogged,'g.', label='full profile')
		plt.plot(x_region_un, y_region_un, 'r.', label='outer 5%')
		plt.xlabel('log ',x_quantity)
		plt.ylabel('log ',y_quantity)
	else:
		plt.plot(xq, yq,'g.', label='full profile')
		plt.plot(x_region, y_region, 'r.', label='outer 5%')
		plt.xlabel(x_quantity)
		plt.ylabel(y_quantity)

	plt.legend(loc=3, fontsize='small')
	plt.savefig(x_quantity +'_v_'+ y_quantity+'_MESA_region.png')
	plt.close()


Rsolar = 6.955e10 #[cgs]
Msolar= 2.0e30 #[kg]
Msolar_cgs= 2.0e33 #[cgs]
#print "mass contained: ", M_contained(MESA_file), "  what unit??"

#plot_region(MESA_file, x='logR', y='logRho', unlog=False)
#plot_region(MESA_file, x='logR', y='mass', unlog=False)

#print "M_contained: ", M_contained(MESA_file)