#!/usr/bin/env python2
import argparse
import glob
import os
import shutil
import subprocess
import sys
m2g_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.join(m2g_path, 'lib'))
from MESAhandling import transfer_inlist
from cfg_parser import get_path_str_option

def get_environment_var(var):
    """Looks up environment variable. If variable is not in the
       system or doesn't have a value it will write an error to
       stderr and exit with a failure code of 1. Otherwise it
       will return the value of the environment variable"""
    value = os.environ.get(var)
    if value is None:
        sys.stderr.write("Must set {} environment variable "
                         "before running this script\n".format(var))
        sys.exit(1)
    return value

#### solution to sdk issue: set the source line in your bashrc

mesa_dir=get_environment_var('MESA_DIR')


### Gets inlist from user. If none is supplied uses mainsequence
parser = argparse.ArgumentParser()
parser.add_argument("inlist", default="inlist_mainsequence", nargs="?",
                    help="The inlist to run through MESA")
args = parser.parse_args()


inlist_name = args.inlist

if os.path.exists(inlist_name):
    inlist_path = inlist_name
else:
    inlist_path = os.path.join(m2g_path, 'data', 'inlists', inlist_name)
if not os.path.exists(inlist_path):
    sys.stderr.write("The file {} could not be found in {}\n"
                     .format(inlist_name, os.path.dirname(inlist_path)))
    sys.exit(1)

# # cp -vr $MESA_DIR/star/work/ to $MESA_DIR/your_work_directory
# star_subdir='work_redgiant'


target_inlist=mesa_dir+'/star/work/inlist_project'
contents = transfer_inlist(inlist_path, target_inlist, mesa_dir, m2g_path)
profile_name = get_path_str_option(contents,
                               'filename_for_profile_when_terminate')
if profile_name is None:
    sys.stderr.write("Inlist must have a "
                     "'file_name_for_profile_when_terminate' set\n")
    sys.exit(1)

print("Inlist moved to {}".format(target_inlist))

os.chdir(mesa_dir+'/star/work/')#+'/'+star_subdir)
subprocess.call('./mk', shell=True)
subprocess.call('./rn', shell=True)

work_dir = os.path.join(mesa_dir, 'star', 'work')
dest_dir = os.path.join(m2g_path, 'data')
if not os.path.exists(work_dir):
    sys.stderr.write("MESA work directory {} does not exists.\n"
                     .format(work_dir))
    sys.exit(1)

if not os.path.exists(dest_dir):
    sys.stderr.write("MESA2HYDRO data directory {} does not exists.\n"
                     .format(dest_dir))
    sys.exit(1)

if not os.path.exists(os.path.join(work_dir, "LOGS")):
    sys.stderr.write("MESA LOGS directory {} does not exists.\n"
                     .format(os.path.join(work_dir, "LOGS")))
    sys.exit(1)

for f in glob.glob(work_dir+'/LOGS/*history*'):
    shutil.copy(f, dest_dir)

profile_path = os.path.join(work_dir, profile_name)
if not os.path.exists(profile_path):
    sys.stderr.write("MESA failed to save a final profile "
                     " to '{}'.\n".format(profile_path))
    sys.exit(1)
