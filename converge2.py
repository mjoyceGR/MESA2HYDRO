#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
#import pygadgetreader as pgr # works- credit this person
import MESAlibjoyce as MJ
import datetime as dt 

def get_MESA_profile_edge(MESA_file,**kwargs):#, strip):
	#print MJ.show_allowed_MESA_keywords(MESA_file)
	strip=bool(kwargs.get('strip',False))
	keyword=str(kwargs.get('quantity','zone'))
	masscut=float(kwargs.get('masscut',0.65))

	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]
	
	try:
		quantity=MJ.get_quantity(MESA_file,keyword)
	except:
		print "Quantity keyword not found. Allowed keywords are:"
		print MJ.show_allowed_MESA_keywords(MESA_file)

	masses = MJ.get_quantity(MESA_file,'mass').astype(np.float)
	Mtot=masses[0]
	bound = masscut*Mtot
	fit_region_indices=np.where( (masses>=bound) )[0]	
	fit_region  = quantity[fit_region_indices]

	print "keyword: ",keyword
	if keyword =='mass':
		#xq_unlogged = [10.0**p for p in xq]
		mf = fit_region
		fit_region = [Mtot-p for p in mf]

	fit_region = np.array(fit_region).astype(float)		
	return fit_region


def get_first_mp(MESA_file, n_p_initial):#(r_array, M_array, rl_init, n_p):
	## tTHIS FUNCTION IS VERY BAD AND NOT CURRENTLY IMPLEMENTED PLEASE IGNORE

	#try:
	r_array = get_MESA_profile_edge(MESA_file, quantity='logR', strip=False)
	#rho_array = get_MESA_profile_edge(MESA_file, quantity='logRho', strip=False)
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

	step_size=0.01#000001
	region=[]
	while len(region) ==0: 
		ru = rl_init + step_size
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]

		#pos_slice=pos_r[region]
#		step_size=step_size*5.0
		print "trying step size: ", step_size
	else:
		#ru = rl_init + step_size
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




def get_region(pos_r, rl, temp_ru, stepsize):
	region=[]
	########## function A ##################
	while len(region) ==0: 
		region=np.where( (pos_r>rl) & (pos_r<temp_ru) )[0]
		temp_ru= temp_ru + stepsize
	else:
		region=np.where( (pos_r>rl) & (pos_r<temp_ru) )[0]
	########## function A ##################
	return region, temp_ru 


def iterate_MESA(r_array, M_array, rl, ru, mp, step_size):
	### kwargs --> stepsize
    #remove ru? JB: no

	r_array=to_log(r_array)
    #r is logged

	#temp_ru=ru+step_size #this is what updates ru each time in the while loop in get_N
	# put step_size here? what is step_size? where am I
	region, temp_ru = get_region( r_array.astype(np.float) , rl, ru, 0.0001) #ru CHANGES in this function, so it must be returned
	mass_slice=get_mass_slice(M_array, region)
	Mshell=get_Mshell(mass_slice)
	n_p=get_np(Mshell,mp)
	n1 = calc_n1(Mshell, mp)
	n2 = calc_n2(rl, temp_ru)

	temp_ru=unlog(temp_ru)
	#r is unlogged

	return n1, n2, temp_ru, n_p

def get_np(Mshell, mp):
	return float(Mshell)/float(mp)

def get_mass_slice(M_array,region):
	M=M_array.astype(np.float)
	return M[region]
 	#return

def get_Mshell(mass_slice):
	########  AMBIGUOUS MASS ORIENTATION, MAY NEEED TO REVERSE THIS
	return abs(  mass_slice.max() - mass_slice.min() )

def calc_n1(Mshell,mp):
	return np.sqrt(Mshell/(12.0*mp))

def calc_n2(rl, ru):
	return np.sqrt(np.pi/12.0)*(ru + rl)/(ru - rl)


def get_N(r_array, M_array, rl, ru, mp, stepsize, n_p_initial,**kwargs):#(MESA_file, mp, step_size):
	n1, n2, keep_ru = 0.0, 100.0, 9999 
	while (n2 -n1) > 0.0: 
	#this will stop as soon as n1 exceeds n2, which is as close as you can get to making them equal given the discrete data
		quant=iterate_MESA(r_array, M_array, rl, ru, mp, stepsize)
		n1 = quant[0]
		n2 = quant[1]
		ru = quant[2]
		n_p = quant[3]

	 	if ru > r_array.max():
	 		break	
	#the next two lines were in an unnecessary "if" within the "while" loop;
	# it's what you want to happen when the while loop is done
	print "condition satisfied:\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru
	keep_ru = ru
	set_R= (keep_ru + rl)/2.0
	return rl, keep_ru, set_R, np.floor(n1), n_p 


def do_converge(MESA_file,r_array, M_array, n_p_initial,stepsize,*args,**kwargs):
	filename=str(kwargs.get('outputfile','vals.dat'))

	r_array=unlog(r_array)
	#r is unlogged

	rl=r_array.min()
	rmax = r_array.max() 

	first_mp=get_first_mp(MESA_file, n_p_initial)
	first_mp=float(kwargs.get('mp',first_mp))
	#first_mp=1e-8

	mp=first_mp
	print "\nfirst mp: ", first_mp, ' mp: ', mp, " first np: ", n_p_initial
	ru = rl + stepsize

	outf=open(filename,"a")
	print >> outf, 'MESA_file: ', MESA_file, ' on ', dt.date.today()
	print >> outf, 'rl				ru				N		n_p				mp			stepsize'
	print >> outf, rl,'\t',ru, "\tN", "\t", n_p_initial, "\t", mp, "\t", stepsize

	while ru < rmax:
		#these both have to be ru but why?
		new_vals= get_N(r_array, M_array, ru, ru, mp, stepsize, n_p_initial) 
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
	outf.close()
	print 'profile data generated in file', filename
	return 0  

def to_log(xq):
	#xq = np.array(xq).astype(float)
	q=[np.log10(p) for p in xq]
	q = np.array(q).astype(float)
	return q 


def unlog(xq):
	#xq = np.array(xq).astype(float)
	q=[10.0**p for p in xq]
	q = np.array(q).astype(float)
	return q 

def plot_region(MESA_file,x_quantity,y_quantity,**kwargs):
	#x_quantity=str(kwargs.get('x','logR'))
	#y_quantity=str(kwargs.get('y','logRho'))
	unlog=bool(kwargs.get('unlog',False))
	try: 
		xq = MJ.get_quantity(MESA_file,x_quantity)
		yq = MJ.get_quantity(MESA_file,y_quantity)
	except:
		print 'quantity keyword not found!\nallowed keys are: '
		print MJ.show_allowed_MESA_keywords(MESA_file)
		exit(0)

	x_region= get_MESA_profile_edge(MESA_file, quantity=x_quantity, strip=False)#[0]
	y_region = get_MESA_profile_edge(MESA_file, quantity=y_quantity , strip=False)#[1]

	if unlog:
		xq_unlogged = unlog(xq)
		yq_unlogged = unlog(yq)
		x_region_un = unlog(x_region)
		y_region_un = unlog(y_region)

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


fname='profile32.data'#'profile175.data' #175, 140, 32
MESA_file="{}".format(fname)

n_p_initial=3000 #72#3072#100,000 
stepsize=0.01#0.0000000000001#1e-6#0.01

r_array = get_MESA_profile_edge(MESA_file, quantity='logR', strip=False)
M_array = get_MESA_profile_edge(MESA_file, quantity='mass', strip=False)

do_converge(MESA_file, r_array, M_array, n_p_initial,stepsize,\
outputfile='test_vals_after_overhaul.dat', mp=0.000102319572081 )#)0.1)
## making it pick "mp" with my estimator is definitely breaking this
