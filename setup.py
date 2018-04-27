#!/usr/bin/env python

from setuptools import setup
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
