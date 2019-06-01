def get_placement_radii(rl, ru, RKstep, force_N, mp, MESA_file, masscut, outf, *args, **kwargs):
	#print "!!!!!!!!!!!!!!"

	input_RKstep=RKstep
	rtot=(mn.MESA_r(MESA_file, masscut)).max()

	fit_region_R   =mn.MESA_r(MESA_file, masscut)
	fit_region_E   =mn.MESA_E(MESA_file, masscut) 

	rmax=fit_region_R.max()


	Mshell_target=target_Mshell(force_N,mp)
	Mshell=0


	rdiff= rmax



	rl=fit_region_R.min()
	temp=rl + RKstep
	ru_mass_loop = temp
	Mshell_integral = 0.0	
	while ru_mass_loop <= rmax:
		
		# if Romberg:
		# 	print "WARNING!!! ROMBERG SWITCH ON!"
		#  	Mshell_integral= Mshell_from_Romberg(oldru, ru, RKstep, MESA_file, masscut, steps=7)	
		# else: 	
		
		#ru_mass_loop = rl + RKstep	
		Mshell_target = 12.0*force_N**2.0*mp


		#temp_RKstep=RKstep
		### keep pushing the upper limit outward until the integrals agree within 5%


		print "loc 1 value of RKstep: ", RKstep
		while ( (100.0*Mshell_integral/Mshell_target) <= 95.0) or ((100.0*Mshell_integral/Mshell_target) >= 105.0): 
			### what should ru be here?
			#print "in while loop...."

			#RKstep=temp_RKstep
			#print "Mshell integral, top of while: ", Mshell_integral
			#print "loc 2 value of RKstep: ", RKstep

			try:
				Mshell_integral = Mshell_from_RK(rl, ru_mass_loop, RKstep, MESA_file, masscut)#, load_unlogged=use_unlog)
				print "Mshell_integral, ru_mass_loop, RKstep:             ",\
			 	 		 ('%.3f'%(100.0*Mshell_integral/Mshell_target)), "    ", ru_mass_loop, "   ", RKstep
			except TypeError:
				break  	 		 

			### adapative step size???
			#print "ru before reset: ", ru_mass_loop
			if (100.0*Mshell_integral/Mshell_target) >= 105.0:
				RKstep=RKstep-0.5*RKstep
				ru_mass_loop = ru_mass_loop - RKstep
				Mshell_integral =0.0
			
			elif (100.0*Mshell_integral/Mshell_target) <= 95.0:
				#
				# experimental
				#
				RKstep = RKstep + RKstep
				#
				ru_mass_loop = ru_mass_loop + RKstep
				Mshell_integral =0.0
			else:
				pass

		## find central point between upper and lower integral limits for converged value
		r_print=(ru_mass_loop + rl)/2.0
		r_nearest,rdex=find_nearest(fit_region_R,r_print)
		u_local=fit_region_E[rdex]
		print >> outf, force_N, r_print, Mshell, u_local

		## reset things for next integral at shell r_n+1
		## experimentally commenting this out		
		#RKstep=input_RKstep
		RKstep=RKstep + RKstep ## doubling for basically no reason
		#RKstep=unlog( (np.log(RKstep)+np.log(input_RKstep))/2.0 )	

		rl = ru_mass_loop
		ru_mass_loop = ru_mass_loop + RKstep
		Mshell_integral = 0
					# break


		#ru_mass_loop= ru_mass_loop + RKstep
		print "integral at %radius ", ru_mass_loop/rmax, "...resetting RKstep =    ", RKstep
		print ""

	# print '',\
	# 'target_Mshell/true_Mshell ', ('%.3f'%(100.0*Mshell_target/Mshell)),r'%',\
	# "   rl, ru:", ('%.6f'%rl),('%.6f'%ru),'   radius',\
	# ('%1.3f'%(100.0*ru/rtot))+str('%'),'    RKstep:',  ('%.6f'%RKstep) 

	#r_place=ru   #(ru+rl)/2.0 ## WARNING! MODIFIED 5/21/19
	return #r_place, Mshell