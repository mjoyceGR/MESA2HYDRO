import mainlib as mn
import converge_funcs as cf
import numpy as np
import sys

def get_placement_radii_orig(rmin, rmax, RKstep, TOL, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):
	lower_convergence_limit = 1.0 - float(TOL)
	upper_convergence_limit = 1.0 + float(TOL)	

	Romberg=kwargs.get("Romberg", False)
	input_RKstep=RKstep
	fit_region_R   =mn.MESA_r(MESA_file, masscut)
	fit_region_E   =mn.MESA_E(MESA_file, masscut) 
	
	Mshell_target=cf.target_Mshell(force_N,mp)

	rl = rmin
	ru_mass_loop = rl # RKstep 
	Mshell_integral = 0.0	

	while ru_mass_loop <= rmax:	
		while ( (Mshell_integral/Mshell_target) <= lower_convergence_limit)\
		   or ( (Mshell_integral/Mshell_target) >= upper_convergence_limit): 

			try:		
				Mshell_integral = cf.Mshell_from_RK(rl, ru_mass_loop, RKstep, MESA_file, masscut)
					#Mshell_integral= Mshell_from_quad(rl, ru_mass_loop, MESA_file, masscut)
				#print "Mshell_integral, ru_mass_loop, stepsize:             ",\
			 	#	  ('%.3f'%(100.0*Mshell_integral/Mshell_target)), "    ", ru_mass_loop, "   ", RKstep
			
			except TypeError:
				print "TypeError in cf.get_placement_radii()"
				break  	 		 

			### adapative step size
			if (Mshell_integral/Mshell_target) >= upper_convergence_limit:
				RKstep=RKstep-0.5*RKstep
				#print "too high, new step=", RKstep
				ru_mass_loop = ru_mass_loop - RKstep
				Mshell_integral =0.0
			

			elif (Mshell_integral/Mshell_target) <= lower_convergence_limit:
				# experimental
				#print "WARNING! no RK doubling"
				RKstep = RKstep + RKstep   #my version 
				
				ru_mass_loop = ru_mass_loop + RKstep
				Mshell_temp=Mshell_integral
				Mshell_integral =0.0
			else:
				pass
				
		###############################################################
		#
		# reset parameters for journey to next shell
		#
		###############################################################
		## find central point between upper and lower integral limits for converged value
		r_print=(ru_mass_loop + rl)/2.0
		r_nearest,rdex=cf.find_nearest(fit_region_R,r_print)
		u_local=fit_region_E[rdex]
		print >> outf, force_N, r_print, Mshell_integral, u_local

		RKtemp = RKstep	
		Mshell_temp=Mshell_integral
		print ('%.3f'%(100.0*Mshell_temp/Mshell_target)),"agreement  for Mshell =",Mshell_temp, \
		       " at ", ru_mass_loop, "    ",\
		        "%.5f"%(ru_mass_loop/rmax), r"%radius  ...current RKstep=","%1.5e"%RKtemp,\
		        " from" , "%1.5e"%input_RKstep#, "  from ", "%1.5e"%RKtemp

		## reset
		RKstep=input_RKstep
		rl = ru_mass_loop
		Mshell_integral = 0
	return #


def John_radii(rmin, rmax, RKstep, TOL, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):

	print "   WARNING!! JOHN'S VERSION     "

	###########################################################
	use_J=2.0
	use_delta=1 #turn on with 1
	use_n=0  #requires use_delta
	use_RK_doubling=1

	ud=float(use_delta)
	un=float(use_n)
	urkd=float(use_RK_doubling)
	###########################################################

	lower_convergence_limit = 1.0 - float(TOL)
	upper_convergence_limit = 1.0 + float(TOL)	

	Romberg=kwargs.get("Romberg", False)
	input_RKstep=RKstep
	fit_region_R   =mn.MESA_r(MESA_file, masscut)
	fit_region_E   =mn.MESA_E(MESA_file, masscut) 
	
	Mshell_target=cf.target_Mshell(force_N,mp)

	rl = rmin
	Mshell_integral = 0.0	
	#delta= ((rmax-rl)/2000.0)/3.0
 	delta=RKstep
	if use_J >= 1.0:
		ru_mass_loop = rl+delta*ud
		print "J working"

	if use_J < 1.0:
		ru_mass_loop = rl
	
	ru_mass_loop_previous = rl

	while ru_mass_loop <= rmax:	
		while ( (Mshell_integral/Mshell_target) <= lower_convergence_limit)\
		   or ( (Mshell_integral/Mshell_target) >= upper_convergence_limit): 

			try:		
				Mshell_integral = cf.Mshell_from_RK(rl, ru_mass_loop, RKstep, MESA_file, masscut)
					#Mshell_integral= Mshell_from_quad(rl, ru_mass_loop, MESA_file, masscut)
				#print "Mshell_integral, ru_mass_loop, stepsize:             ",\
			 	#	  ('%.3f'%(100.0*Mshell_integral/Mshell_target)), "    ", ru_mass_loop, "   ", RKstep
			
			except TypeError:
				print "TypeError in cf.get_placement_radii()"
				break  	 		 

			### adapative step size
			if (Mshell_integral/Mshell_target) >= upper_convergence_limit:
				RKstep=RKstep-0.5*RKstep
				#print "too high, new step=", RKstep                  	
				if use_J < 3.0:
					ru_mass_loop = ru_mass_loop - (1.0-ud)*RKstep - ud*delta
				#should be equivalent to next line ru_mass_loop = ru_mass_loop - RKstep
				if use_J >= 3.0:
					ru_mass_loop = ru_mass_loop_previous 
					ru_mass_loop_previous=ru_mass_loop - RKstep
					print "J >=2.0 activated, ru_mass_loop_previous = ", ru_mass_loop_previous

				Mshell_integral =0.0
				delta=delta/2.0
			
			elif (Mshell_integral/Mshell_target) <= lower_convergence_limit:

				ru_mass_loop_previous = ru_mass_loop
				# experimental
				#print "WARNING! no RK doubling"
				#RKstep = RKstep + urkd*RKstep   #my version 
				RKstep = RKstep + RKstep   #my version 

				try:
					n=np.floor(1.0/( Mshell_integral/Mshell_target))
				except ZeroDivisionError:
					n = 1.0
				if n >= 2.0:
					n = n - 1.0
				
				
				if use_J < 2.0:
					ru_mass_loop = ru_mass_loop + RKstep
				if use_J >= 2.0:
					ru_mass_loop = ru_mass_loop + (1.0-ud)*RKstep + ud*(un*float(n)+(1.0-un)*1.0)*delta
				
				Mshell_temp=Mshell_integral
				Mshell_integral =0.0
			else:
				pass
		

		###############################################################
		#
		# reset parameters for journey to next shell
		#
		###############################################################
		## find central point between upper and lower integral limits for converged value
		r_print=(ru_mass_loop + rl)/2.0
		r_nearest,rdex=cf.find_nearest(fit_region_R,r_print)
		u_local=fit_region_E[rdex]
		print >> outf, force_N, r_print, Mshell_integral, u_local

		RKtemp = RKstep	
		Mshell_temp=Mshell_integral
		print ('%.3f'%(100.0*Mshell_temp/Mshell_target)),"agreement  for Mshell =",Mshell_temp, \
		       " at ", ru_mass_loop, "    ",\
		        "%.5f"%(ru_mass_loop/rmax), r"%radius  ...current RKstep=","%1.5e"%RKtemp,\
		        " from" , "%1.5e"%input_RKstep#, "  from ", "%1.5e"%RKtemp

		## reset
		RKstep=input_RKstep
		#delta= ((ru_mass_loop-rl))/3.0
		delta = RKstep

		rl = ru_mass_loop
                ru_mass_loop_previous = ru_mass_loop
		ru_mass_loop = ru_mass_loop + ud*delta 
		Mshell_integral = 0
	return #