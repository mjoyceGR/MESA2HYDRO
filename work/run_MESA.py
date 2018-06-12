#!/usr/bin/env python
import subprocess
import sys
import os
m2g_path=os.environ['MESA2HYDRO_ROOT']
sys.path.append(m2g_path+'/lib/')
import MESAhandling as MJ


#### solution to sdk issue: set the source line in your bashrc

mesasdk_root=os.environ['MESASDK_ROOT']
m2g_path=os.environ['MESA2HYDRO_ROOT']
mesa_dir=os.environ['MESA_DIR']

print "m2g_path", m2g_path 
#sys.exit(0)

# # cp -vr $MESA_DIR/star/work/ to $MESA_DIR/your_work_directory
# star_subdir='work_redgiant' 	

template_inlist=m2g_path+ '/data/inlists/inlist_mainsequence'
target_inlist=mesa_dir+'/star/work/inlist_project' 

MJ.transfer_inlist(template_inlist, target_inlist) 
os.chdir(mesa_dir+'/star/work/')#+'/'+star_subdir)

subprocess.call('./mk', shell=True)
subprocess.call('./rn', shell=True)

subprocess.call('cp '+mesa_dir+'/star/work/LOGS/*history* '+m2g_path+'/data/', shell=True)
subprocess.call('cp '+mesa_dir+'/star/work/profile_* '+m2g_path+'/data/', shell=True)

