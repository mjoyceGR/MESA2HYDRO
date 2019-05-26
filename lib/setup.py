from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

######################################################

print "addition to setup.py 5/23/19-- TO BE INTEGRATED: contact Lianne about combining this with setup.py"
## command: python setup.py build_ext --inplace

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran phantread_module.f90 -c -o phantread_module.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)
shared_obj_comp = 'gfortran pygfunc.f90 -c -o pygfunc.o -O3 -fPIC'
print shared_obj_comp
system(shared_obj_comp)

ext_modules = [Extension(# module name:
                         'pygfunc',
                         # source file:
                         ['pygfunc.pyx'],
                         libraries=['gfortran'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3'],
                         # other files to link to
                         extra_objects=['phantread_module.o', 'pygfunc.o'])]

setup(name = 'pygfunc',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = cythonize(ext_modules))
