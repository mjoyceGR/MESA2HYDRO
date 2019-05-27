from numpy import linspace, empty
from numpy cimport ndarray as ar
from cpython cimport array
import array

cdef extern from "pygfunc.h":
    void c_gfunc(int* ngas, double* mgas, double* x, double* y, double* z, double* h, double* u, double* msink, double* hsoft_sink);
    #void c_gfunc(int* ngas, double* msink); #, double* mgas, double* x, double* y, double* z, double* h, double* u, double* msink);


print("\n\nloc 1 pyx, calling...\n")

#def to_cdef(int ngas ,double central_point_mass):
def to_cdef(int ngas, double[::1] mgas,\
 			double[::1] x, double[::1] y, double[::1] z,\
 			double[::1] h, double[::1] u, double central_point_mass, double hsoft_sink):

# double[::1] x, double[::1] y, double[::1] z, double[::1] h, double[::1] u, double central_point_mass):
    
    ## this function's job is to pass things to c_gfunc()
    #cdef:
    #    int c_ngas = ngas
    #    array.array c_mgas = array.array('i', mgas)

    # figure out basic C syntax or whatever this is
    #print("calling, loc 4")
    #print "(loc 4), type(central_point_mass)", type(central_point_mass)
    ### these used to be all zeros

    # zero, 1, or ngas??
    #c_gfunc(&ngas,  &mgas[0], &x[0], &y[0], &z[0], &h[0], &u[0], &central_point_mass)
    c_gfunc(&ngas, &mgas[0], &x[0], &y[0], &z[0], &h[0], &u[0], &central_point_mass, &hsoft_sink)
    print("called, loc (7)")
