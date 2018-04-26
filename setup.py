#!/usr/bin/env python

from distutils.core import setup
setup(name='MESA2GADGET',
      packages=['MESA2GADGET', 'MESA2GADGET.work'],
      package_dir={'MESA2GADGET': 'src', 'MESA2GADGET.work': 'work'}, 
      package_data={'MESA2GADGET': ['data/*/*', 'work/N_mp_combinations.dat', 'work/NR_files/*.dat',
                                    'work/recovery_images/*.png', 'out/sample_MESA_output/*.data']},
      version='0.1.0',
      description='TBD',
      author='Meridith Joyce',
      author_email='TBD',
      url='https://github.com/mjoyceGR/MESA2GADGET',
      download_url='https://github.com/mjoyceGR/MESA2GADGET/archive/0.1.0.tar.gz',
      keywords=['mesa', 'gadget', 'astronomy'],
      classifiers=[],
      scripts=[
            'work/run_MESA.py',
            'work/alg.py',
            'work/make_IC.py',
            'work/confirm_density_profile.py',
            'work/confirm_mass_profile.py',
            'work/confirm_HR.py'
      ],
      requires=[
        'h5py',
        'scipy',
        'healpy',
        'matplotlib'
      ]
)
