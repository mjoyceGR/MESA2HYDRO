#!/usr/bin/env python3

import os
import sys
import setuptools


################################


VERSION = '1.1.3'
PYGFUNC_F90_FILE = "lib/pygfunc.f90"
PYGFUNC_PYF_FILE = "lib/pygfunc.pyf"
WRITE_DATA_PHANTOM_F90_FILE = "lib/write_data_phantom.f90"


setup_data = dict(
    name='MESA2HYDRO',
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
    setup_requires=['numpy'],
    install_requires=['numpy', 'scipy', 'matplotlib', 'healpy', 'h5py']
    )


if len(sys.argv) >= 2 and ('--help' in sys.argv[1:] or
        sys.argv[1] in ('--help-commands', 'egg_info', '--version',
                        'clean')):
    # For these actions, NumPy is not required.
    #
    # They are required to succeed without Numpy for example when
    # pip is used to install Scipy when Numpy is not yet present in
    # the system.
    try:
        from setuptools import setup
    except ImportError:
        from distutils.core import setup

else:

    from numpy.distutils.core import setup, Extension
    from numpy import f2py

    if not os.path.exists(PYGFUNC_PYF_FILE):
        print("Making {} file".format(PYGFUNC_PYF_FILE))
        f2py.run_main([PYGFUNC_F90_FILE, '-m', 'pygfunc', '-h', PYGFUNC_PYF_FILE])


    fortran_ext = Extension(
        name = 'pygfunc',
        sources = [PYGFUNC_PYF_FILE, WRITE_DATA_PHANTOM_F90_FILE, PYGFUNC_F90_FILE]
        )

    setup_data['ext_modules'] = [fortran_ext]


setup(**setup_data)
