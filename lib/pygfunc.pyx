from numpy import linspace, empty
from cpython cimport array
import array

cdef extern from "pygfunc.h":
    void c_gfunc(int* ngas, double* mgas, double* x, double* y, double* z, double* h, double* u, double* msink, double* hsoft_sink);


def to_cdef(int ngas, double[::1] mgas,\
 			double[::1] x, double[::1] y, double[::1] z,\
 			double[::1] h, double[::1] u, double central_point_mass, double hsoft_sink):
    c_gfunc(&ngas, &mgas[0], &x[0], &y[0], &z[0], &h[0], &u[0], &central_point_mass, &hsoft_sink)
    print("called, loc (7)")
