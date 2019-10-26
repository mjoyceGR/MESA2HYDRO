#!/usr/bin/env python3

from setuptools import setup, find_packages

#import os
#import sys
#import shutil
#import subprocess
#from numpy import get_include, f2py

################################


# TODO: Add "local" option again that adds MESA2HYDRO_ROOT to python path and installs dependencies

setup(name='MESA2HYDRO',
      version='0.1.2146',
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
      install_requires=['numpy', 'scipy', 'matplotlib', 'healpy', 'h5py', 'setuptools']
      )


"""
PKG_DIR = os.path.abspath(os.path.dirname(__file__))
PYGFUNC_F90_FILE = "{}/lib/pygfunc.f90".format(PKG_DIR)
PYGFUNC_PYF_FILE = "{}/lib/pygfunc.pyf".format(PKG_DIR)
WRITE_DATA_PHANTOM_F90_FILE = "{}/lib/write_data_phantom.f90".format(PKG_DIR)



if len(sys.argv) > 2 and sys.argv[1] == 'install':
    f2py.run_main([PYGFUNC_F90_FILE, '-m', 'pygfunc', '-h', PYGFUCN_PYF_FILE])
    f2py.run_main(['-c', PYGFUNC_PYF_FILE, "--opts='-std=f2003'", WRITE_DATA_PHANTOM_F90_FILE, PYGFUNC_F90_FILE, '-m', 'pygfunc'])
    
"""

"""
## must explicitly include all f90, h, and pyx files
setup(name='MESA2HYDRO',
      scripts=['work/confirm_density_profile.py',
               'work/confirm_mass_profile.py',
               'work/run_conversion.py',
               'uninstall_MESA2HYDRO.py'],
      cmdclass={'build_ext': build_ext},
      include_dirs=[get_include()],
      ext_modules=cythonize(ext_modules, compiler_directives={'language_level': "3"}),
      packages=['MESA2HYDRO',
                'MESA2HYDRO.work',
                'MESA2HYDRO.lib'],
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

      version='0.1.2146',
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
"""
