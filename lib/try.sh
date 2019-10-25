
gfortran -c write_data_phantom.f90 -std=f2003


## -m is name the module pygfunc
f2py3 pygfunc.f90 -m pygfunc -h pygfunc.pyf
f2py3 -c pygfunc.pyf --opt='-std=f2003' write_data_phantom.f90 pygfunc.f90 -m pygfunc
