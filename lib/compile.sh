# ! /bin/bash

rm pygfunc.c

gfortran -fPIC -c pygfunc.f90  phantread_module.f90

python setup.py build_ext --inplace

python ../work/run.py ../work/mainsequence.cfg