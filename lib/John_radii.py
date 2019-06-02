import mainlib as mn
import converge_funcs as cf
import numpy as np

def John_radii(rmin, rmax, RKstep, TOL, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):
	
	lower_convergence_limit = 1.0 - float(TOL)
	upper_convergence_limit = 1.0 + float(TOL)	

	Romberg=kwargs.get("Romberg", False)
	input_RKstep=RKstep
	fit_region_R   =mn.MESA_r(MESA_file, masscut)
	fit_region_E   =mn.MESA_E(MESA_file, masscut) 
	
	Mshell_target=cf.target_Mshell(force_N,mp)
	#deleted Mshell=0 (not used)

	rl = rmin #switched order here and input to rmin, since rl is the one that varies

	ru_mass_loop = rl
	ru_mass_loop_previous = ru_mass_loop
	Mshell_integral = 0.0	

	delta_numerator= ((rmax-rl)/2000.0)/3.0
	precision = 0.0
	delta = delta_numerator


	while ru_mass_loop <= rmax:	
		rl_temp = rl
		precision = 0.0

		while ( (Mshell_integral/Mshell_target) <= lower_convergence_limit) or ((Mshell_integral/Mshell_target) >= upper_convergence_limit): 
			#try:
				#Mshell_integral_previous = Mshell_integral
				#Mshell_integral = Mshell_integral + Mshell_from_RK(rl_temp, ru_mass_loop, RKstep, MESA_file, masscut)
			Mshell_integral =  cf.Mshell_from_RK(rl, ru_mass_loop, RKstep, MESA_file, masscut)

			#except TypeError:
			#	break  	 		 

			if (Mshell_integral/Mshell_target) >= upper_convergence_limit:

				ru_mass_loop = ru_mass_loop_previous 
				
				#precision= precision + 1.0
				delta = delta/2.0


				if delta/RKstep <= 100.0:
					RKstep = RKstep/3.0


			elif (Mshell_integral/Mshell_target) <= lower_convergence_limit:

				try:
					n=np.floor(1.0/( Mshell_integral/Mshell_target))
				except ZeroDivisionError:
					n = 1.0

				if n >= 2.0:
					n = n - 1.0

				ru_mass_loop_previous = ru_mass_loop
				ru_mass_loop = ru_mass_loop + (float(n))*delta 

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

		RKtemp = RKstep	 # these temp vars aren't used (usefully)
		Mshell_temp=Mshell_integral
		print ('%.3f'%(100.0*Mshell_temp/Mshell_target)),"agreement  for Mshell =",Mshell_temp, \
		       " at ", ru_mass_loop, "    ",\
		        "%.5f"%(ru_mass_loop/rmax), r"%radius  ...current RKstep=","%1.5e"%RKtemp,\
		        " from" , "%1.5e"%input_RKstep#, "  from ", "%1.5e"%RKtemp

		RKstep=input_RKstep

		delta_numerator= ((ru_mass_loop-rl))/3.0
		precision = 0.0
		delta = delta_numerator

		rl = ru_mass_loop
                ru_mass_loop_previous = ru_mass_loop
		ru_mass_loop = ru_mass_loop + delta 
		#Mshell_integral = 0.0
	return #