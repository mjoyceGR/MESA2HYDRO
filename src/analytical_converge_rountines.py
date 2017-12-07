########################################################################################################
#---------------------------------- old routines, maybe not useful? ------------------------------------
########################################################################################################
binary_file="AGB_from_glass_binary.dat"
def n1_n2_gadget(binary_file,rl,ru):
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


def Mshell_analytical(rl, ru, A,B,C):
	rdiff=(ru - rl)
	f1=0.5*A*rdiff**2.0  + A*B*rdiff  + A*B**2.0 * np.log( abs(rdiff - B) ) + (1.0/3.0)*C*rdiff**3.0
	return 4.0*np.pi*f1

def n1_n2_from_integral_fit(rl,ru, A, B, C):
	M_shell=Mshell_analytical(rl, ru, A, B, C)
	n_p=M_shell/mp

	n1 = np.sqrt(M_shell/(12.0*mp))
	n2 = calc_n2(rl, ru)

	return n1, n2
##############################################################################################################
##############################################################################################################
##############################################################################################################


# add optional arguments for strong/weak tolerances
# def get_mp(MESA_file, n_p):
# 	try:
# 		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, False)
# 	except:
# 		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, True) # not sure if this is the right order here

# 	rl=r_array.min() #3.1
# 	# rmax = r_array.max() #5.5#50*begin
# 	mp=get_first_mp(r_array, M_array, rl, n_p)
# 	return mp



#################### ???
# def get_np_per_shell(mp, MESA_file, n_p):
# 	try:
# 		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, False)
# 	except:
# 		r_array, rho_array, M_array = get_MESA_profile_edge(MESA_file, True) # not sure if this is the right order here

# 	rl=r_array.min() #3.1
# 	# rmax = r_array.max() #5.5#50*begin
# 	mp=get_first_mp(r_array, M_array, rl, n_p)
# 	return mp



################# break above line into single thing that returns mp ##################################################