#!/usr/bin/env python
# import numpy as np
# import matplotlib.pyplot as plt
# import MESAlibjoyce as MJ
# import converge_funcs as cf
# import io_lib as rw
import subprocess
import os

mesasdk_root=os.environ['MESASDK_ROOT']
subprocess.call('bash', shell=True)
subprocess.call('source '+mesasdk_root+'/bin/mesasdk_init.sh',shell=True)
m2g_path=os.environ['MESA2GADGET_ROOT']
mesa_dir=os.environ['MESA_DIR']

star_subdir='work_whitedwarf' 	# cp -vr $MESA_DIR/star/work/ to $MESA_DIR/your_work_directory


template_inlist=m2g_path+ '/data/whitedwarf/inlist_whitedwarf' 
target_inlist=mesa_dir+'/'+star_subdir+'/inlist_project' 

MJ.copy_inlist_template(template_inlist, target_inlist) #works


os.chdir(mesa_dir+'/'+star_subdir)

subprocess.call('./mk', shell=True)
subprocess.call('./rn', shell=True)

subprocess.call('cp ./LOGS/*history* '+m2g_path+'/out/', shell=True)
subprocess.call('cp profile_* '+m2g_path+'/out/', shell=True)

