#!/usr/bin/env python
import subprocess
import sys
import os
m2g_path=os.environ['MESA2GADGET_ROOT']
sys.path.append(m2g_path+'/src/')
try:
    import MESA2GADGET.mesalib.MESAlibjoyce as MJ
except ImportError:
    print("Problem with MESA2GADGET installation")
    print("To use this package please run sudo python setup.py install")
    print("or set your PYTHONPATH environment variable to the directory")
    print("MESA2GADGET is in (pointing it directly to MESA2GADGET still causes problems)")
    exit(1)

#### solution to sdk issue: set the source line in your bashrc

mesasdk_root=os.environ['MESASDK_ROOT']
#m2g_path=os.environ['MESA2GADGET_ROOT']
mesa_dir=os.environ['MESA_DIR']

# cp -vr $MESA_DIR/star/work/ to $MESA_DIR/your_work_directory
star_subdir='work_redgiant' 	

template_inlist=m2g_path+ '/data/redgiant/inlist_redgiant'
target_inlist=mesa_dir+'/'+star_subdir+'/inlist_project' 

MJ.transfer_inlist(template_inlist, target_inlist) 
os.chdir(mesa_dir+'/'+star_subdir)

subprocess.call('./mk', shell=True)
subprocess.call('./rn', shell=True)

subprocess.call('cp '+mesa_dir+'/'+star_subdir+'/LOGS/*history* '+m2g_path+'/out/', shell=True)
subprocess.call('cp '+mesa_dir+'/'+star_subdir+'/profile_* '+m2g_path+'/out/', shell=True)

