#!/usr/bin/env python

import os
import pip
import sys
import shutil
import subprocess
from setuptools import setup

BUILD_CMDS = [
    'local', # my command
    'build',
    'build_py',
    'build_ext',
    'build_clib',
    'build_scripts',
    'install',
    'install_lib',
    'install_headers',
    'install_scripts',
    'install_data',
    'sdist',
    'register',
    'bdist',
    'bdist_dump',
    'bdist_rmp',
    'upload',
    'check',
    'develop',
    'bdist_egg',
    'easy_install',
    'egg_info',
    'bdist_wheel']
CLEAN_CMDS = BUILD_CMDS+['clean']
PKG_DIR = os.path.abspath(os.path.dirname(__file__))

required_pkgs = [
    'numpy',
    'h5py',
    'scipy',
    'healpy',
    'matplotlib']

if len(sys.argv) < 2:
    print("For local setup\nusage: setup.py local\n\nFor normal setup")
    has_cmd = False
elif len(sys.argv) > 2:
    has_cmd = "help" not in sys.argv[2]
else:
    has_cmd = True

if len(sys.argv) > 1 and sys.argv[1] == 'local':
    if len(sys.argv) > 2 and sys.argv[2] == "--help":
        print("Sets up the project for local use.")
        print("Downloads and installs needed libraries")
        print("Sets up PYTHONPATH variable")
        exit(0)

    install_failed = False
    for pkg in required_pkgs:
        try:
            pip.main(['install', pkg])
        except:
            isntall_failed = True
            print("Couldn't install {}. Please run again with sudo to install it."
                  .format(pkg))

    if not install_failed:
        print("All required packages are installed")
    
    with open(os.path.expanduser("~/.bashrc"), "ra+") as f:
        file_text = f.read()
        if "ADDED BY MESA2GADGET" not in file_text:
            f.write("\n\n# ADDED BY MESA2GADGET\n")
            f.write("export PYTHONPATH=$PYTHONPATH:{}\n\n".format(os.path.dirname(PKG_DIR)))
            print("Updated python path")
            print("Please restart (or source ~/.bashrc) your terminal")
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
      keywords=['mesa', 'gadget', 'astronomy'],
      classifiers=[],
      entry_points={
        'console_scripts': [
            'run_MESA = MESA2GADGET.work.run_MESA:main_func',
            'alg = MESA2GADGET.work.alg:main_func'
      ]},
      requires=required_pkgs,
      zip_safe=False
)

