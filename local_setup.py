#!/usr/bin/env python3

import os
import sys
import subprocess

################################

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
    print("Appending {} to PYTHONPATH in ~/.bashrc".format(PYPATH_DIR))
    with open(os.path.expanduser('~/.bashrc'), 'w+') as f:
        f.write('# ADDED BY MESA2HYDRO FOR LOCAL DEVELOPMENT')
        f.write('export PYTHONPATH=$PYTHONPATH:~/')
    
subprocess.run("pip install -r requirements.txt", shell=True)

# install fortran locally
subprocess.run('./install_phantom.sh', shell=True)
