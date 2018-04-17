#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#import pygadgetreader as pgr # works- credit this person
import MESAlibjoyce as MJ
import datetime as dt 
import random as rand
import healpy as hp

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

	#print "keyword: ",keyword
	if keyword =='mass':
		#xq_unlogged = [10.0**p for p in xq]
		#mf = fit_region
		fit_region = outer_mass(Mtot, fit_region)#[Mtot-p for p in mf]

	return np.array(fit_region).astype(float)




def outer_mass(Mtot,fit_region):
	mf=fit_region
	fit_region = [Mtot-p for p in mf]
	return np.array(fit_region).astype(float)


def get_mp_given_N(r_array, M_array, n_p_initial):#(r_array, M_array, rl_init, n_p):
	## tTHIS FUNCTION IS VERY BAD AND NOT CURRENTLY IMPLEMENTED PLEASE IGNORE
	pos_r=r_array
	rl=r_array.min() #in unlogged form
	ru = rl #rl_init

	step_size=1.5#0.0001#000001
	region=[]
	while len(region) ==0: 
		ru = ru + step_size
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
	else:
		region=np.where( (pos_r>rl) & (pos_r<ru) )[0]
		pos_slice=pos_r[region]
	print "step size selected: ", step_size

	mass_slice=get_mass_slice(M_array, region)
	Mshell=get_Mshell(mass_slice)
	mp = float(Mshell)/float(n_p_initial)
	print "selected mass per particle mp=",mp
	return mp




def iterate_MESA(r_array, M_array, rl, ru, mp, step_size):
	#DO NOT WANT THIS LOGGED!!!!!!!!1

	### kwargs --> stepsize
    #remove ru? JB: no
	#print 'step_size: ', step_size

	# r_array=to_log(r_array)
 #    #r is logged
	# rl=to_log(rl)
	# ru=to_log(ru)

	#temp_ru=ru+step_size #this is what updates ru each time in the while loop in get_N
	# put step_size here? what is step_size? where am I
	region, temp_ru = get_region( r_array.astype(np.float) , rl, ru, step_size)#0.001) #ru CHANGES in this function, so it must be returned

	mass_slice=get_mass_slice(M_array, region)
	Mshell=get_Mshell(mass_slice)
	n_p=get_np(Mshell,mp)
	n1 = calc_n1(Mshell, mp)
	n2 = calc_n2(rl, temp_ru)

	#temp_ru=unlog(temp_ru)
	#r is unlogged

	return n1, n2, temp_ru, n_p


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


def get_N(r_array, M_array, rl, ru, mp, stepsize,**kwargs):#(MESA_file, mp, step_size):
	n1, n2, keep_ru = 0., 100., 9999. 
	#print "before while", n1, n2
	while (n2 -n1) > 0.0: 
	#this will stop as soon as n1 exceeds n2, which is as close as you can get to making them equal given the discrete data
	 	if ru > r_array.max():
	 		break	
		quant=iterate_MESA(r_array, M_array, rl, ru, mp, stepsize)
		n1 = quant[0]
		n2 = quant[1]
		ru = quant[2]
		n_p = quant[3]
	keep_ru = ru
	set_R= (keep_ru + rl)/2.0
	print "converged:\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru, '\tr_mid: ', set_R, '\tn_p', n_p
	return rl, keep_ru, set_R, min(np.floor(n1), np.floor(n2)), n_p 

def get_N_continuous(rl, rmax, A, B, C, mp, stepsize, **kwargs):
	ru_list=np.arange(rl,rmax,stepsize)
	n1, n2, shell_ru=0.,100.,9999.
	ru = rl
	while (n2 > n1):
	#for i in range(len(ru_list)):
		if (ru > rmax): 
			break
		# ru = ru_list[i]
		Mshell=Mshell_from_integral(rl, ru, A,B,C)
		n1=calc_n1(Mshell, mp)
		n2=calc_n2(rl, ru)
		ru = ru + stepsize

	shell_ru=ru
	set_R=(shell_ru+rl)/2.0	
	print "converged:\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru, '\tr_mid: ', set_R#, '\tn_p', n_p	
	return rl, shell_ru, set_R, min(np.floor(n1), np.floor(n2))



# def get_N_continuous(rl, rmax, A, B, C, mp, stepsize,tolerance, **kwargs):
# 	ru_list=np.arange(rl,rmax,stepsize)
# 	shell_ru=9999
# 	for i in range(len(ru_list)):
# 		ru = ru_list[i]
# 		Mshell=Mshell_from_integral(rl, ru, A,B,C)
# 		n1=calc_n1(Mshell, mp)
# 		n2=calc_n2(rl, ru)
# 		if abs(n1 -n2) <= tolerance:#0.1 
# 			#print "tolerance met\tn1: ", n1, "\tn2: ", n2, "\tr_l: ", rl, "\tr_u: ", ru
# 			shell_ru=ru

# 	set_R=(shell_ru+rl)/2.0		
# 	# figure out how to return the intermediate radius here too		
# 	return shell_ru, set_R, min(np.floor(n1), np.floor(n2))



def do_converge(MESA_file,r_array, M_array, n_p_initial,stepsize,*args,**kwargs):
	filename=str(kwargs.get('outputfile','vals.dat'))

	r_array=unlog(r_array)
	#r is unlogged
	#print r_array

	rl=r_array.min()
	rmax = r_array.max() 

	first_mp=get_first_mp(MESA_file, n_p_initial)
	first_mp=float(kwargs.get('mp',first_mp))
	#first_mp=1e-8

	mp=first_mp
	print "\nfirst mp: ", first_mp, ' mp: ', mp, " first np: ", n_p_initial
	ru = rl + stepsize

	outf=open(filename,"a")
	print >> outf, 'MESA_file: ', MESA_file, ' on ', dt.datetime.now()
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
	logged=[]
	for i in range(len(xq)):
		try:
			q=np.log10(xq[i])
			logged.append(q)
		except TypeError:
			pass
	logged = np.array(logged).astype(float)
	return logged 


def unlog(xq):
	unlogged=[]
	for i in range(len(xq)):
		#print xq
		try:
			q=10.0**float(xq[i])
			unlogged.append(q)
		except ValueError:
			pass
	unlogged = np.array(unlogged).astype(float)
	return unlogged 



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

#------------------------------ analytic ----------------------------------------

def get_curve(r_array,rho_array,guess_a,guess_b,guess_c):
	def model_func(x,a,b,c):
			## this is a 1/r density profile
		return a*(1.0/(x+b)) + c

	p=[guess_a,guess_b,guess_c]
	params, cov = curve_fit(model_func, r_array, rho_array, p0=p)
	a, b, c = params

	def y(x):
		a,b,c=params
		return a*(1.0/(x+b)) + c

	return a,b,c,y 
#------------------------------------------------------------------------------------
def density_integral(rl,ru, A,B,C):
	rdiff=(ru - rl)
	f1=0.5*A*rdiff**2.0  + A*B*rdiff  + A*B**2.0 * np.log10( abs(rdiff - B) ) + (1.0/3.0)*C*rdiff**3.0
	return f1

### is it necessary to do this density integral or should I just fit the mass profile and pass that?


def Mshell_from_integral(rl, ru, A,B,C):
	Mshell=4.0*np.pi*density_integral(rl, ru, A, B, C)
	return Mshell


def approximate_mp(rl, ru, n_p, A, B, C):
	Mshell=4.0*np.pi*density_integral(ru,rl,A,B,C)
	mp=Mshell/n_p
	return mp




###########################################################################
#
# healpix
#
############################################################################
def get_coords(N,r_mid, rmax):#,file_index):
    r_mid=float(r_mid); rmax=float(rmax)
    ##Need N to be (4) 8 or 16, and must change mp per shell to maintain that, probably
    NSIDE = N #closest power of 2 to N
    NSIDE=int(NSIDE)

    theta=random_theta()
    ipix_array=np.arange(hp.nside2npix(NSIDE)) #this is just a straight up array of particle IDs
    x=[]
    y=[]
    z=[]
    for i in range(len(ipix_array)):
        ipix=ipix_array[i]
        coord=hp.pixelfunc.pix2vec(NSIDE, ipix, nest=True)
        #print >> outf, coord[0], coord[1], coord[2]
        x.append(coord[0])
        y.append(coord[1])
        z.append(coord[2])

    # no fam this rotate shell thing is some serious math
    xd, yd, zd = rotate_shell(x,y,z,theta,"about_z")
    xe, ye, ze = rotate_shell(xd,yd,zd,theta,"about_y")
    xf, yf, zf = rotate_shell(xe,ye,ze,theta,"about_x")

    xf = np.array(xf).astype(float)
    yf = np.array(yf).astype(float)
    zf = np.array(zf).astype(float)

    xf = to_physical(xf, r_mid/rmax)
    yf = to_physical(yf, r_mid/rmax)
    zf = to_physical(zf, r_mid/rmax)

    xf = np.array(xf).astype(float)
    yf = np.array(yf).astype(float)
    zf = np.array(zf).astype(float)

    #print "x: ",xf, "\n\ny: ", yf, "\n\nz: ",zf
    return xf.flatten(), yf.flatten(), zf.flatten()


def to_physical(xq, r_mid):
    xf = [ (r_mid*x_i) for x_i in xq]
    return xf

def rotate_shell(x_array, y_array, z_array, theta, direction, **kwargs):
    vec=np.array( [x_array, y_array, z_array])

    Rx=np.matrix( [\
    [1.0, 0.0, 0.0],\
    [0, np.cos(theta), -np.sin(theta)],\
    [0, np.sin(theta), np.cos(theta)]\
    ])

    Ry=np.matrix( [\
    [np.cos(theta), 0.0, np.sin(theta)],\
    [0.0, 1.0, 0.0],\
    [-np.sin(theta), 0.0, np.cos(theta)]\
    ])

    Rz=np.matrix( [\
    [np.cos(theta), -np.sin(theta), 0.0],\
    [np.sin(theta), np.cos(theta), 0.0],\
    [0.0, 0.0, 1.0]\
    ])

    if str(direction) == "about_z":
        new=Rz*vec
    elif str(direction) == "about_y":
        new=Ry*vec
    else:
        new=Rx*vec

    new_x=np.array(new[0].transpose()).astype(float)
    new_y=np.array(new[1].transpose()).astype(float)
    new_z=np.array(new[2].transpose()).astype(float)

    return new_x, new_y, new_z

def random_theta():
    theta=rand.random()*2.0*np.pi #.random gives random float between 0 and 1
    return theta


# x_array=np.array([3,3,3,3])
# y_array = 2.0+0.0*x_array
# z_array = 5.0+0.0*x_array

# x,y,z=rotate_shell(x_array,y_array,z_array, random_theta())

# print "\none column? x: ", x
# print "\nx[0]:", x[0]
# print "\nwhat? x[0][0]:", x[0][0]

# for i in range(len(x)):
#     print "array[i] etc: ", np.sqrt(x_array[i]**2.0+y_array[i]**2.0+z_array[i]**2.0)
#     print "x[i] etc", np.sqrt(x[i]**2.0 + y[i]**2.0 + z[i]**2.0)


def to_rad(theta):
    theta=theta*np.pi/180.0
    return theta


def to_array(array):
    return np.array(array).astype(np.float)


def min_index(array):
	mn,idx = min( (array[i],i) for i in xrange(len(array)) )
	return idx