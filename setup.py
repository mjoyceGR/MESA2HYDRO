#!/usr/bin/env python

import os
import shutil
import subprocess
from setuptools import setup
PKG_DIR = os.path.abspath(os.path.dirname(__file__))
CMD = "cd {} && git clone https://github.com/ldocao/pygadgetic.git".format(PKG_DIR)
print("Try removing any previous packages")
try:
    shutil.rmtree(os.path.join(PKG_DIR, "pygadgetic"))
except OSError:
   print("Package repo didn't need to be removed") 
try:
    shutil.rmtree(os.path.join(PKG_DIR, "mesalib", "pygadgetic"))
except OSError:
   print("No packages needed to be removed") 

print("Pullilng dependent package")
print("Calling " + CMD)
subprocess.check_call(CMD, shell=True, executable='/bin/bash')
print("Moving package files into mesalib")
shutil.copytree(os.path.join(PKG_DIR, "pygadgetic", "pygadgetic"),
                os.path.join(PKG_DIR, "mesalib", "pygadgetic"))
print("Cleaning up git package")
shutil.rmtree(os.path.join(PKG_DIR, "pygadgetic"))

setup(name='MESA2GADGET',
      packages=['MESA2GADGET', 'MESA2GADGET.mesalib', 'MESA2GADGET.work'],
      package_dir={'MESA2GADGET': '', 'MESA2GADGET.work': 'work',
                   'MESA2GADGET.mesalib': 'mesalib',
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
      keywords=['mesa', 'gadget', 'astronomy'],
      classifiers=[],
      entry_points={
        'console_scripts': [
            'run_MESA = MESA2GADGET.work.run_MESA:main_func',
            'alg = MESA2GADGET.work.alg:main_func'
      ]},
      requires=[
        'h5py',
        'scipy',
        'healpy',
        'matplotlib'
      ],
      zip_safe=False
)

