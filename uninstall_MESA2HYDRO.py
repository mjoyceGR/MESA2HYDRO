#!/usr/bin/python2

import os
import sys
import shutil
from subprocess import check_output, CalledProcessError


scripts = ['confirm_density_profile.py',
           'confirm_mass_profile.py',
           'run_conversion.py']


script_locations = []
for script in scripts:
    try:
        out = check_output("which {} 2>/dev/null".format(script), shell=True)
    except CalledProcessError as err:
        continue
    if out:
        script_locations.append(out)


for location in script_locations:
    location = location[:-1]
    if not os.path.isfile(location):
        continue
    result = raw_input("Remove {}(y/N)?".format(location))
    if result.upper().startswith("Y"):
        try:
            os.remove(location)
        except OSError as err:
            print(err)
            print("Try running again with sudo")
            sys.exit(1)

try:
    import MESA2HYDRO
except:
    print("MESA2HYDRO isn't installed, not removing it")
    sys.exit(0)

package_path = MESA2HYDRO.__path__[0]
print("MESA2HYDRO is installed at {}".format(package_path))

package_dir = os.path.dirname(package_path)
all_packages = os.listdir(package_dir)

related_files = []
for package_file in all_packages:
    if "MESA2HYDRO" in package_file:
        related_files.append(os.path.join(package_dir, package_file))

for location in related_files:
    if not os.path.exists(location):
        print("Skipping {} because doesn't exist".format(location))
        continue

    result = raw_input("Remove {}?(y/N)".format(location))
    if result.upper().startswith("Y"):
        try:
            if os.path.isdir(location):
                shutil.rmtree(location)
            else:
                os.remove(location)
        except OSError as err:
            print(err)
            print("Try running again with sudo")
            sys.exit(1)

file_path = os.path.abspath(__file__)

result = raw_input("Remove this file ({})?(y/N)".format(file_path))
print(result)
if result.upper().startswith("Y"):
    try:
        os.remove(file_path)
    except OSError as err:
        print(err)
        print("Try running again with sudo")
        sys.exit(1)