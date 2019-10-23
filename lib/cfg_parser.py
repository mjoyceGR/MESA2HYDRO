#!/usr/bin/env python

#************************************************************************
#
#   Copyright (C) 2019  M. Joyce, L. Lairmore, D. J. Price
#
#   See MESA2HYDRO/LICENSE
#
#************************************************************************

'''

Contains: command line and config file argument parsers,
          path manipulation functions 

'''

import sys
import os
import re
import argparse
import tempfile

DEBUG = False
MESA_PKG_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
print("Package directory: {}".format(MESA_PKG_DIR))

try:
    m2g_save_path=os.path.abspath(os.environ.get('MESA2HYDRO_ROOT'))
except AttributeError:
    print("Example .cfg files can be found at {}/data/config_file_examples".format(MESA_PKG_DIR))
    print("Error: Use absolute paths in .cfg file or set MESA2HYDRO_ROOT in .bashrc")
    sys.exit(1)

if m2g_save_path is None:
    print("Environmental variable 'MESA2HYDRO_ROOT' isn't set")
    print("Storing output data in default directory root {}"
          .format(MESA_PKG_DIR))
    m2g_save_path = MESA_PKG_DIR
else:
    print("Files will be saved in {}\n\n".format(m2g_save_path))
    if not os.path.exists(m2g_save_path):
        os.makedirs(m2g_save_path)


def relative_to_root(path):
    """If the file doesn't exists as is, make the path relative to 
       MESA2HYDRO_ROOT environment variable"""
    if not os.path.exists(path):
        path = os.path.join(m2g_save_path, path)
    return path



if DEBUG:
    print("PACKAGE PATH: {}".format(MESA_PKG_DIR))

######################################################
#
# module CONFIG PARSER
#
#####################################################
class OptionTypeUnknown:
    pass

class Types:
    STR = "str"
    BOOL = "bool"
    FLOAT = "float"
    INT = "int"

class OptionInputs:

    def __init__(self, options, description=None,
                 cmd_line=True, config_file=True, verbose=True):
        self.verbose = verbose
        self.options = options
        self.config_file = config_file
        self.cmd_line = cmd_line
        self.parser = argparse.ArgumentParser(description=description)
        self.parser.add_argument('--defaults', action='store_true',
                                 help='List file\'s default configuration')
        if self.config_file:
            self.parser.add_argument('config_file', default=None, nargs='?',
                                     help='Path to configuration file. '
                                          'Overrides file defaults. '
                                          'File defaults can be discovered with '
                                          '\'--defaults\'')

        if self.cmd_line:
            self.make_cmd_args()
        self.args = None
        self.file_out = None
        self.out = {}
        for name, default in options.items():
            self.out[name] = default

    def show(self, msg):
        if self.verbose:
            print(msg)

    def init_args(self):
        if self.args:
            return self.args
        self.args = self.parser.parse_args()
        if self.args.defaults:
            print("Default file configuration: \n")
            for config, value in self.options.items():
                print("{} = {}".format(config, value))
            sys.exit(0)
        return self.args

    def init_config(self):
        self.init_args()
        if self.args.config_file is None:
            return {}

        if not os.path.exists(self.args.config_file):
            config_file = path_from_package(self.args.config_file)
        else:
            config_file = self.args.config_file

        if not os.path.exists(config_file):
            print("{} doesn't exist".format(self.args.config_file))
            sys.exit(1)

        with open(config_file, 'r') as cf:
            self.file_out = cf.read()

        lines = self.file_out.splitlines()
        lines = [line for line in lines
                 if not re.match('\s*#', line) and '=' in line]
        self.file_out = "\n".join(lines)
        return self.file_out

    def get_config_file(self):
        if self.file_out is None:
            self.init_config()
            
        return self.file_out

    @staticmethod
    def get_type(option, value):
        if isinstance(value, str):
            return Types.STR
        elif isinstance(value, bool):
            return Types.BOOL
        elif isinstance(value, int):
            return Types.INT
        elif isinstance(value, float):
            return Types.FLOAT
        else:
            raise OptionTypeUnknown("Option {} has default value {}. "
                                    "Unable to determine option's type from default value. "
                                    "Supported add-able option types are str, bool, float, and int."
                                    .format(option, value))

    def make_cmd_args(self):
        config_group = self.parser.add_argument_group(title="Config Options",
                                                      description="Will overwrite config file values.")
        for name, default in self.options.items():
            opt_t = self.get_type(name, default)
            if opt_t == Types.BOOL:
                config_group.add_argument('--{}'.format(name), default=None,
                                          action='store_{}'.format('false' if default else 'true'),
                                          help="Sets {} to {}".format(name, not default))

            elif opt_t == Types.STR:
                config_group.add_argument('--{}'.format(name),
                                          help="Sets {} to value supplied".format(name))

            elif opt_t == Types.INT:
                config_group.add_argument('--{}'.format(name),
                                          type=int,
                                          help="Sets {} to value supplied".format(name))
            elif opt_t == Types.FLOAT:
                config_group.add_argument('--{}'.format(name),
                                          type=float,
                                          help="Sets {} to value supplied".format(name))
                                  

    def read_config_file(self):
        self.init_config()
        for name, default in self.options.items():
            if has_option(self.file_out, name):
                opt_t = self.get_type(name, default)
                if opt_t == Types.STR:
                    self.out[name] = get_path_option(self.file_out, name)
                elif opt_t == Types.BOOL:
                    self.out[name] = get_bool_option(self.file_out, name)
                elif opt_t == Types.INT:
                    self.out[name] = get_int_option(self.file_out, name)
                elif opt_t == Types.FLOAT:
                    self.out[name] = get_float_option(self.file_out, name)

        return self.out

    def read_cmd_args(self):
        self.init_args()
        arguments = vars(self.args)
        for name, value in arguments.items():
            if name not in self.options:
                continue
            if value is not None:
                self.out[name] = value

        return self.out

    def get_configs(self):
        self.init_args()
        if self.config_file and self.args.config_file:
            self.show("Using config values from {}".format(self.args.config_file))
            self.read_config_file()

        if self.cmd_line:
            self.read_cmd_args()

        self.show("Script configurations are:")
        for name, value in sorted(self.out.items()):
            self.show("\t{} = {}".format(name, value))

        return self.out


def path_from_package(path):
    return os.path.join(MESA_PKG_DIR, path)


def has_option(config_file, name):
    return bool(re.search('^\s*{}\s*=\s*(?P<value>[\w\./]+)'.format(name),
                          config_file, re.MULTILINE))


def get_path_str_option(config_file, name, default=None):
    m = re.search('^\s*{}\s*=\s*[\',"](?P<value>[\w\./]+)[\',"]'.format(name),
                  config_file, re.MULTILINE)
    if m:
        value = m.group('value')
        if DEBUG:
            print("Using config file value for {} of {}".format(name, value))
        return value
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default


def get_path_option(config_file, name, default=None):
    m = re.search('^\s*{}\s*=\s*(?P<value>[\w\./]+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        value = m.group('value')
        if DEBUG:
            print("Using config file value for {} of {}".format(name, value))
        return value
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default


def get_str_option(config_file,name, default=None):
    m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        value = m.group('value')
        if DEBUG:
            print("Using config file value for {} of {}".format(name, value))
        return value
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default


def get_float_option(config_file,name, default=None):
    m = re.search('^\s*{}\s*=\s*(?P<value>[\w\.-]+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        try:
            value = m.group('value')
            value = float(value)
            if DEBUG:
                print("Using config file value for {} of {}".format(name, value))
            return value
        except ValueError:
            print("{} has a misformed value. {} is not a recognized float"
                  .format(name, value))
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default


def get_int_option(config_file,name, default=None):
    m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        try:
            value = m.group('value')
            value = int(value)
            if DEBUG:
                print("Using config file value for {} of {}".format(name, value))
            return value
        except ValueError:
            print("{} has a misformed value. {} is not a recognized float"
                  .format(name, value))
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default


def get_bool_option(config_file,name, default=None):
    m = re.search('^\s*{}\s*=\s*(?P<value>\w+)'.format(name),
                  config_file, re.MULTILINE)
    if m:
        # case insensitive true
        value = m.group('value')
        if value.upper() == "TRUE":
            if DEBUG:
                print("Using config file value for {} of {}".format(name, True))
            return True
        else:
            if value.upper() == "FALSE":
                if DEBUG:
                    print("Using config file value for {} of {}".format(name, False))
                return False
        print("{} needs to have true or false and is {}".format(name, value))
        return default
    else:
        if re.search(name, config_file):
            print("{} was malformed in config file".format(name))
        return default

# end module cfg_parser

# test
