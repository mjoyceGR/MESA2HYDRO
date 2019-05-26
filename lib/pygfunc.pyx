from numpy import linspace, empty
from numpy cimport ndarray as ar
from cpython cimport array
import array

cdef extern from "pygfunc.h":
    void c_gfunc(int ngas, int* mgas, double* x, double* y, double* z, double* h, double* u, double* msink);


def to_cdef(int ngas, int[::1] mgas, double[::1] x, double[::1] y, double[::1] z, double[::1] h, double[::1] u, double[::1] central_point_mass):
    ## this function's job is to pass things to c_gfunc()
    #cdef:
    #    int c_ngas = ngas
    #    array.array c_mgas = array.array('i', mgas)

    # figure out basic C syntax or whatever this is
    print("calling")
    c_gfunc(ngas, &mgas[0], &x[0], &y[0], &z[0], &h[0], &u[0], &central_point_mass[0])
    print("called")
