#!/usr/bin/env python
import sys
import os
import re
import argparse
from MESA2GADGET import MESA_PKG_DIR
######################################################
#
# module CONFIG PARSER
#
#####################################################


def command_line_parser():
##############################################################
#
# Command line argument input--- MUST ALSO EDIT HERE IF ADDING NEW PARAMETER
#
##############################################################


    VALID_CONFIGS = {
        'check_MESA_profile': True,
        'MESA_file': path_from_package('out/sample_MESA_output/profile_mainsequence_logE.data'),
        'make_NR_file': False,
        'loaded_NR_filename': 'work/NR_files/saveNR_ms.dat',
        'new_NR_filename': 'latest_NR.dat',

        'make_IC_file': False,
        'loaded_IC_filename': 'ms_logE',
        'new_IC_filename': 'latest_IC',
        'IC_format_type': 'binary',

        'masscut': 0.95,
        'N': 8,
        'mp': 1e-7,
        'mp_cgs': 1.988e26,
        'stepsize': 2.45e6,

        'try_reload': False,
        'png_tag': 'latest',

        'reload_bin_size': 70.0,
        'use_bins_from_NR_file': False,
        'which_dtype':'f'}



    parser = argparse.ArgumentParser(description='Program for converting MESA output to Gadget simulation files')
    parser.add_argument('--config-file',
                        help='Path to configuration file')
    parser.add_argument('--defaults', action='store_true',
                        help='Lists the scripts default parameters and exits')
    config_args = parser.add_argument_group("Configuration")


    # Boolean option flags set opposite of file default
    config_args.add_argument('--check-MESA',
                             action='store_true' if VALID_CONFIGS['check_MESA_profile'] else 'store_false',
                             help='Sets check-MESA value to {}'
                                  .format(not VALID_CONFIGS['check_MESA_profile']))
    config_args.add_argument('--make-NR-file',
                             action='store_true' if VALID_CONFIGS['make_NR_file'] else 'store_false',
                             help='Sets make-NR-file to {}'
                                  .format(not VALID_CONFIGS['make_NR_file']))
    config_args.add_argument('--make-IC-file',
                             action='store_true' if VALID_CONFIGS['make_IC_file'] else 'store_false',
                             help='Sets make-IC-file value to {}'
                                  .format(not VALID_CONFIGS['make_IC_file']))
    config_args.add_argument('--try-reload',
                             action='store_true' if VALID_CONFIGS['try_reload'] else 'store_false',
                             help='Sets try-reload to {}'
                                  .format(not VALID_CONFIGS['try_reload']))
    config_args.add_argument('--format-type', default=VALID_CONFIGS['IC_format_type'],
                             help='Type of the format (binary...)')
    config_args.add_argument('--MESA-file',
                             default=VALID_CONFIGS['MESA_file'],

                             help='Path to input MESA output')
    config_args.add_argument('--masscut', default=VALID_CONFIGS['masscut'],
                             type=float,
                             help='Sets masscut')
    config_args.add_argument('--N', default=VALID_CONFIGS['N'],
                             type=int,
                             help='Sets N')
    config_args.add_argument('--mp', default=VALID_CONFIGS['mp'],
                             type=float,
                             help='Set the mp value in Msolar units')
    config_args.add_argument('--mp_cgs', default=VALID_CONFIGS['mp_cgs'],
                             type=float,
                             help='Set the mp value in cgs units')
    config_args.add_argument('--step-size', default=VALID_CONFIGS['stepsize'],
                             type=float,
                             help='Sets the step size to the given value. Units in cm')
    # config_args.add_argument('--star-type', default=VALID_CONFIGS['startype'],
    #                          help='Set start type')
    config_args.add_argument('--loaded_NR_filename', default=None,
                             help='Set the NR file')
    config_args.add_argument('--loaded_IC_filename', default=None,
                             help='Set the IC file')
    config_args.add_argument('--usNRbins',default=VALID_CONFIGS['use_bins_from_NR_file'],
                             help='Use bins from NR file')
    return parser
 




def path_from_package(path):
    return os.path.join(MESA_PKG_DIR, path)

def get_path_option(config_file, name, default):
    m = re.search('^\s*{}\s*=\s*(?P<value>[\w\./]+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        value = m.group('value')
        print("Using config file value for {} of {}".format(name, value))
        return value
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default

def get_str_option(config_file,name, default):
    m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        value = m.group('value')
        print("Using config file value for {} of {}".format(name, value))
        return value
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default

def get_float_option(config_file,name, default):
    m = re.search('^\s*{}\s*=\s*(?P<value>[\w\.-]+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        try:
            value = m.group('value')
            value = float(value)
            print("Using config file value for {} of {}".format(name, value))
            return value
        except ValueError:
            print("{} has a misformed value. {} is not a recognized float"
                  .format(name, value))
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default

def get_int_option(config_file,name, default):
    m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        try:
            value = m.group('value')
            value = int(value)
            print("Using config file value for {} of {}".format(name, value))
            return value
        except ValueError:
            print("{} has a misformed value. {} is not a recognized float"
                  .format(name, value))
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default

def get_bool_option(config_file,name, default):
    m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        # case insensitive true
        value = m.group('value')
        if value.upper() == "TRUE":
            print("Using config file value for {} of {}".format(name, True))
            return True
        else:
            if value.upper() == "FALSE":
                print("Using config file value for {} of {}".format(name, False))
                return False
        print("{} needs to have true or false and is {}".format(name, value))
        return default
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default