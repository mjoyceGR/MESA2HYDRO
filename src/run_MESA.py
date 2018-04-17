#!/usr/bin/env python
import MESAlibjoyce as MJ
import subprocess
import os

#### solution to sdk issue: set the source line in your bashrc


mesasdk_root=os.environ['MESASDK_ROOT']
m2g_path=os.environ['MESA2GADGET_ROOT']
mesa_dir=os.environ['MESA_DIR']

star_subdir='work_whitedwarf' 	# cp -vr $MESA_DIR/star/work/ to $MESA_DIR/your_work_directory

#print "going???"
template_inlist=m2g_path+ '/data/whitedwarf/inlist_whitedwarf'
target_inlist=mesa_dir+'/'+star_subdir+'/inlist_project' 
#print 'made it here' 

MJ.transfer_inlist(template_inlist, target_inlist) #works
os.chdir(mesa_dir+'/'+star_subdir)

subprocess.call('./mk', shell=True)
subprocess.call('./rn', shell=True)

subprocess.call('cp ./LOGS/*history* '+m2g_path+'/out/', shell=True)
subprocess.call('cp profile_* '+m2g_path+'/out/', shell=True)

