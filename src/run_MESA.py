#!/usr/bin/env python
import MESAlibjoyce as MJ
import subprocess
import os

#### solution to sdk issue: set the source line in your bashrc

mesasdk_root=os.environ['MESASDK_ROOT']
m2g_path=os.environ['MESA2GADGET_ROOT']
mesa_dir=os.environ['MESA_DIR']

# cp -vr $MESA_DIR/star/work/ to $MESA_DIR/your_work_directory
star_subdir='work_AGB' 	

template_inlist=m2g_path+ '/data/AGB/inlist_agb_from_mod'
target_inlist=mesa_dir+'/'+star_subdir+'/inlist_project' 

MJ.transfer_inlist(template_inlist, target_inlist) 
os.chdir(mesa_dir+'/'+star_subdir)

subprocess.call('./mk', shell=True)
subprocess.call('./rn', shell=True)

subprocess.call('cp ./LOGS/*history* '+m2g_path+'/out/', shell=True)
subprocess.call('cp profile_* '+m2g_path+'/out/', shell=True)

