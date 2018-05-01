#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/src/')
import MESA2GADGET.mesalib.MESAlibjoyce as MJ
import MESA2GADGET.mesalib.converge_funcs as cf
###################################################
M_to_solar=1.988*10.0**33.0 ## g/Msolar
R_to_solar=6.957*10.0**10.0 ## cm/Rsolar
###################################################

savefig=False

ftag=str(raw_input('enter MESA history identifier (ex: "redgiant"): '))
histfile='../out/sample_MESA_output/history_'+ftag+'.data'
sf=str(raw_input('save? (y for yes): '))
if sf=='y':
	savefig=True
else:
	savefig=False
tag=ftag


try:
	logL=MJ.get_quantity(histfile,'log_L')
	logT=MJ.get_quantity(histfile, 'log_Teff')
	ax, fig=plt.subplots()
	plt.plot(logT, logL,'g-')
	# plt.xlim(3.51,3.53)
	# plt.ylim(3.79,3.9)
	plt.xlabel('Log log_Teff')
	plt.ylabel('Log L')
	plt.title('HR diagram for history_'+ftag)
	plt.gca().invert_xaxis()
	if savefig:
		plt.savefig('HRtest_'+str(tag)+'.png')
	plt.show()
	plt.close()
except IOError:
	print 'file ../out/sample_MESA_output/history'+ftag+'.data not found\n'
	print 'available files: '
	subprocess.call('ls ../out/sample_MESA_output/history*', shell=True)
