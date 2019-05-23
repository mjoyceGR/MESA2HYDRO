#!/usr/bin/env python
import os
import pip
import sys
import shutil
import subprocess
from setuptools import setup

################################
import subprocess
## lib/setup.py needs to be executable
#subprocess.call("python lib/setup.py build_ext --inplace", shell=True)



PKG_DIR = os.path.abspath(os.path.dirname(__file__))

if len(sys.argv) > 1 and sys.argv[1] == 'local':
    if len(sys.argv) > 2 and sys.argv[2] == "--help":
        print("Sets up the project for local use.")
        print("Downloads and installs needed libraries")
        print("Sets up PYTHONPATH variable")
        exit(0)

    try:
        subprocess.check_call(["{}/install.sh".format(PKG_DIR)], shell=True)
    except:
        print("Failed to install all necessary packages")

    exit(0)

setup(name='MESA2GADGET',
      packages=['MESA2GADGET', 'MESA2GADGET.lib', 'MESA2GADGET.work'],
      package_dir={'MESA2GADGET': '', 'MESA2GADGET.work': 'work',
                   'MESA2GADGET.lib': 'lib',
                   'MESA2GADGET.data': 'data', 'MESA2GADGET.out': 'out'}, 
      package_data={'': ['data/*/*', 'work/N_mp_combinations.dat', 'work/NR_files/*.dat',
                         'work/recovery_images/*.png', 'out/sample_MESA_output/*.data']},
      include_package_data=True,
      version='0.1.0',
      description='TBD',
      author='Meridith Joyce',
      author_email='TBD',
      url='https://github.com/mjoyceGR/MESA2GADGET',
      download_url='https://github.com/mjoyceGR/MESA2GADGET/archive/0.1.0.tar.gz',
      keywords=['mesa', 'hydro', 'astronomy'],
      classifiers=[],
      entry_points={
        'console_scripts': [
            'run_MESA = MESA2GADGET.work.run_MESA:main_func',
            'alg = MESA2GADGET.work.alg:main_func'
      ]},
      requires=required_pkgs,
      zip_safe=False
)


