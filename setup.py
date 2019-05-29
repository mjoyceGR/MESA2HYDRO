#!/usr/bin/env python
import os
import pip
import sys
import shutil
import subprocess
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from numpy import get_include
from os import system

################################
import subprocess
## lib/setup.py needs to be executable
#subprocess.call("python lib/setup.py build_ext --inplace", shell=True)



PKG_DIR = os.path.abspath(os.path.dirname(__file__))

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran -fdefault-real-8 write_data_phantom.f90 -c -o write_data_phantom.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)
shared_obj_comp = 'gfortran pygfunc.f90 -c -o pygfunc.o -O3 -fPIC'
print shared_obj_comp
system(shared_obj_comp)

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

ext_modules = [Extension(# module name:
                         'pygfunc',
                         # source file:
                         ['lib/pygfunc.pyx'],
                         libraries=['gfortran'],
                         # other files to link to
                         extra_objects=['lib/write_data_phantom.o', 'lib/pygfunc.o'])]

setup(name='MESA2GADGET',
      cmdclass={'build_ext': build_ext},
      include_dirs=[get_include()],
      ext_modules=cythonize(ext_modules),
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
      #requires=required_pkgs,
      zip_safe=False
)


