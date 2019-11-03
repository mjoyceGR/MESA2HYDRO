#!/usr/bin/env python3

import os
import shutil
import sys
import subprocess

################################

LOCAL_BASH = '~/.bashrc'
BACKUP_BASH = '~/.bashrc.bak'

FULL_PATH_BASH = os.path.expanduser(LOCAL_BASH)
FULL_PATH_BASH_BAK = os.path.expanduser(BACKUP_BASH)

PKG_DIR = os.path.abspath(os.path.dirname(__file__))
PYPATH_DIR = os.path.dirname(PKG_DIR)


PYGFUNC_F90_FILE = "lib/pygfunc.f90"
PYGFUNC_PYF_FILE = "lib/pygfunc.pyf"
WRITE_DATA_PHANTOM_F90_FILE = "lib/write_data_phantom.f90"


python_path = os.environ.get('PYTHONPATH') or ''
python_paths = [os.path.abspath(path) for path in python_path.split(':')]

if PYPATH_DIR not in python_paths:
    print("{} not in PYTHONPATH, needed for local development of MESA2HYDRO"
          .format(PYPATH_DIR))
    print("Appending {} to PYTHONPATH in {}".format(PYPATH_DIR, FULL_PATH_BASH))
    shutil.copyfile(FULL_PATH_BASH, FULL_PATH_BASH_BAK)
    with open(FULL_PATH_BASH, 'a') as f:
        f.write('# ADDED BY MESA2HYDRO FOR LOCAL DEVELOPMENT\n')
        f.write('export PYTHONPATH=$PYTHONPATH:{}'.format(PYPATH_DIR))
    print("Restart your terminal to complete local installation")
    
subprocess.run("pip install -r requirements.txt", shell=True)

# install fortran locally
subprocess.run('./install_phantom.sh', shell=True)
