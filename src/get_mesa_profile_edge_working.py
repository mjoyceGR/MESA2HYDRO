def get_MESA_profile_edge(MESA_file,**kwargs):#, strip):
	#print MJ.show_allowed_MESA_keywords(MESA_file)
	strip=bool(kwargs.get('strip',False))
	quantity=str(kwargs.get('quantity',zone))

	if strip:
		MESA_file=MJ.strip_MESA_header(MESA_file,MESA_file,n=5)[1]

	mass = MJ.get_quantity(MESA_file,'mass')
	#print "maximum mass: ", mass.max()

	logR = MJ.get_quantity(MESA_file,'logR')
	logrho = MJ.get_quantity(MESA_file,'logRho')
	logrho = logrho.astype(np.float)
	#rho = [10.0**p for p in logrho]

	mass = mass.astype(np.float)
	Mtot=mass[0]
	# print "Mtot: ", Mtot
	# print "last mass entry (total M??): ", Mtot, "\tlargest mass value: ", max(mass),\
	#  "\tmass value at most inner R: ", mass[np.argmin(logR)]


	## THE "mass" VARIABLE ALREADY DESCRIBES THE MASS CONTAINED
	cummulative=mass

	fit_region_rho=[]
	fit_region_R=[]
	fit_region_M=[]
	for m in range(len(cummulative)):
		#print "m value: ",m,'\t',cummulative[m]
		if cummulative[m] <= 0.05*cummulative[np.argmax(logR)]:
			fit_region_rho.append( (10.0**float(logrho[m])))#/100000.0 )) 
			fit_region_R.append(10.0**float(logR[m]))
			fit_region_M.append(float(cummulative[m]))
			#m+=1 
	fit_region_R   = np.array(fit_region_R).astype(np.float)
	fit_region_rho = np.array(fit_region_rho).astype(np.float)
	fit_region_M   = np.array(fit_region_M).astype(np.float)
	
	#print len(fit_region_R), len(fit_region_M), len(fit_region_rho)

	return fit_region_R, fit_region_rho, fit_region_M 
