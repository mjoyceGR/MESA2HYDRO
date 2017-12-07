#!/usr/bin/python3

import subprocess
import argparse
import os

def run(function):
    def wrapper(work_path):
        commands = ["cd {}".format(work_path)]
        commands.extend(function(work_path))
        for command in commands:
            try:
                subprocess.call(command, shell=True)
            except:
                print("{} command failed".format(command))
    return wrapper


@run
def compile_MESA(work_path):
    return ["echo 'compiling'"]


@run
def execute_MESA(work_path):
    return ["echo 'executing'"]

this_dir = os.path.dirname(os.path.realpath(__file__))

parser = argparse.ArgumentParser(description='MESA run tool')

parser.add_argument('-c', '--compile', action='store_true', help='Compiles MESA')
parser.add_argument('-e', '--execute', action='store_true', help='Executes MESA')
parser.add_argument('-b', '--both', action='store_true', help='Compiles then executes MESA')
parser.add_argument('-p', '--path', default=this_dir, help='Path MESA is stored')

args = parser.parse_args()

if args.compile:
    compile_MESA(args.path)
elif args.execute:
    execute_MESA(args.path)
elif args.both:
    compile_MESA(args.path)
    execute_MESA(args.path)
else:
    print("Please specify to compile or execute or both for MESA system")

