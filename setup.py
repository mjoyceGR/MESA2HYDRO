#!/usr/bin/env python
import os
try:
    from pip._internal import main
except ImportError:
    from pip import main
import sys
import shutil
import subprocess
from distutils.core import setup
from distutils.extension import Extension
import subprocess
try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
except ImportError:
    main(['install', 'cython'])
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
from numpy import get_include
from os import system

################################



PKG_DIR = os.path.abspath(os.path.dirname(__file__))

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran -fdefault-real-8 lib/write_data_phantom.f90 -c -o lib/write_data_phantom.o -O3 -fPIC'
print("Compile library: " + fortran_mod_comp)
system(fortran_mod_comp)
shared_obj_comp = 'gfortran lib/pygfunc.f90 -c -o lib/pygfunc.o -O3 -fPIC'
print("Compile interface: " + shared_obj_comp)
system(shared_obj_comp)

if len(sys.argv) > 1 and sys.argv[1] == 'install':
    main(['install', 'numpy', 'h5py', 'scipy', 'healpy', 'matplotlib'])
    subprocess.call("rm {}/lib/*.mod {}/lib/*.o".format(PKG_DIR, PKG_DIR), shell=True)
    subprocess.call("cython -a {}/lib/pygfunc.pyx".format(PKG_DIR), shell=True)
    subprocess.call("python setup.py build_ext --inplace", shell=True)


### only need the things that wrap with python here? AKA not write_data_phantom because it's not pyx
ext_modules = [Extension(# module name:
                         'pygfunc',

                         # source file:
                         ['lib/pygfunc.pyx'],
                         libraries=['gfortran'],

                         # other files to link to
                         extra_objects=['lib/write_data_phantom.o', 'lib/pygfunc.o'])]

## must explicitly include all f90, h, and pyx files
setup(name='MESA2HYDRO',
      scripts=['setup.py',
               'work/confirm_density_profile.py',
               'work/confirm_mass_profile.py',
               'work/run.py'],
      cmdclass={'build_ext': build_ext},
      include_dirs=[get_include()],
      ext_modules=cythonize(ext_modules),
      packages=['MESA2HYDRO', 'MESA2HYDRO.lib', 'MESA2HYDRO.work'],

      package_dir={'MESA2HYDRO': '',
                   'MESA2HYDRO.work': 'work',
                   'MESA2HYDRO.lib': 'lib',
                   'MESA2HYDRO.data': 'data',
                   'MESA2HYDRO.out': 'out',
                   'MESA2HYDRO.DOCUMENTATION':'DOCUMENTATION'},

      package_data={'': ['data/*/*',
                         'work/*.cfg',
                         'out/NR_files/*.dat',
                         'out/sample_MESA_output/*.data',

                         'lib/write_data_phantom.f90',
                         'lib/pygyfunc.f90',
                         'lib/pygyfunc.f90',
                         'lib/pygyfunc.pyx',
                         'lib/pygyfunc.h',]},

      version='0.1.2131',
      description='Convert 1D stellar structure models to 3D particle distributions using the HEALPix spherical tessellation algorithm',
      long_description='Convert 1D stellar structure models to 3D particle distributions using the HEALPix spherical tessellation algorithm\
      \n\
      see academic paper https://iopscience.iop.org/article/10.3847/1538-4357/ab3405/meta\
      \n\
      and User\'s Guide\
      \n\
      https://github.com/mjoyceGR/MESA2HYDRO/blob/master/DOCUMENTATION/MESA2HYDRO_users_guide.pdf\n\
      for more information',
      author='Meridith Joyce',
      author_email='Meridith.Joyce@anu.edu.au',
      url='https://github.com/mjoyceGR/MESA2HYDRO',
      download_url='https://github.com/mjoyceGR/MESA2HYDRO/archive/0.1.0.tar.gz',
      keywords=['MESA', 'stellar evolution', 'stellar structure','astronomy', 'stellar modeling','hydrodynamical ICs','HEALPix'],
      classifiers=[]
)


