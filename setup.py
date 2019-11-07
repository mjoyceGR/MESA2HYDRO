#!/usr/bin/env python3

import os
import sys

from numpy.distutils.core import setup, Extension
from numpy import f2py

################################

VERSION = '1.1.1'
PYGFUNC_F90_FILE = "lib/pygfunc.f90"
PYGFUNC_PYF_FILE = "lib/pygfunc.pyf"
WRITE_DATA_PHANTOM_F90_FILE = "lib/write_data_phantom.f90"

if not os.path.exists(PYGFUNC_PYF_FILE):
    print("Making {} file".format(PYGFUNC_PYF_FILE))
    f2py.run_main([PYGFUNC_F90_FILE, '-m', 'pygfunc', '-h', PYGFUNC_PYF_FILE])


fortran_ext = Extension(
    name = 'pygfunc',
    sources = [PYGFUNC_PYF_FILE, WRITE_DATA_PHANTOM_F90_FILE, PYGFUNC_F90_FILE]
    )

setup(name='MESA2HYDRO',
      version=VERSION,
      description=('Convert 1D stellar structure models to 3D particle '
                   'distributions using the HEALPix spherical tessellation algorithm'),
      long_description=('Convert 1D stellar structure models to 3D particle '
                        'distributions using the HEALPix spherical tessellation algorithm\n'
                        'see academic paper https://iopscience.iop.org/article/10.3847/1538-4357/ab3405/meta\n'
                        'and User\'s Guide\n'
                        'https://github.com/mjoyceGR/MESA2HYDRO/blob/master/DOCUMENTATION/MESA2HYDRO_users_guide.pdf\n'
                        'for more information'),
      author='Meridith Joyce',
      author_email='Meridith.Joyce@anu.edu.au',
      url='https://github.com/mjoyceGR/MESA2HYDRO',
      # Make sure this still exists?
      # download_url='https://github.com/mjoyceGR/MESA2HYDRO/archive/0.1.0.tar.gz',
      keywords=['MESA', 'stellar evolution', 'stellar structure','astronomy',
                'stellar modeling', 'hydrodynamical ICs', 'HEALPix'],
      packages=['MESA2HYDRO',
                'MESA2HYDRO.work',
                'MESA2HYDRO.lib'],
      scripts=['work/confirm_density_profile.py',
               'work/confirm_mass_profile.py',
               'work/run_conversion.py',
               'uninstall_MESA2HYDRO.py'],
      package_dir={'MESA2HYDRO': '',
                   'MESA2HYDRO.work': 'work',
                   'MESA2HYDRO.lib': 'lib',
                   'MESA2HYDRO.data': 'data',
                   'MESA2HYDRO.out': 'out',
                   'MESA2HYDRO.DOCUMENTATION': 'DOCUMENTATION'},
      package_data={'': ['data/*/*',
                         'work/*.cfg',
                         'out/NR_files/*.dat',
                         'out/sample_MESA_output/*.data',
                         'lib/write_data_phantom.f90',
                         'lib/pygfunc.f90',
                         'DOCUMENTATION/*']},
      ext_modules=[fortran_ext]
      )
