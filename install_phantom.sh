#!/bin/bash


pushd lib

# -m is name the module pygfunc
# generates the pygfunc.pyf file
f2py3 pygfunc.f90 -m pygfunc -h pygfunc.pyf

# generates .so file
f2py3 -c pygfunc.pyf write_data_phantom.f90 pygfunc.f90 -m pygfunc

popd
