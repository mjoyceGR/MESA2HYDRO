from numpy import linspace, empty
from numpy cimport ndarray as ar

cdef extern from "pygfunc.h":
    void c_gfunc(int* ngas,int* mgas, double* x, double* y, double* z, double*h, double* u, double* msink);
    ## arguments: (double* a, int* n, int* m, double* a, double* b, double* c)


def to_cdef(ngas, mgas, x, y, z, h, u, central_point_mass, *args, **kwargs):
	## this function's job is to pass things to c_gfunc()

	#double x, double a=-10.0, double b=10.0,
    #cdef:
    #     ar[double] ax = linspace(a, b, n)
    #     ar[double,ndim=2] c = empty((n, n), order='F')
 

    ## cast my python types as c types ?
    cdef:
    	int ngas = ngas
    	int mgas = mgas
		ar[double] xcoords = x
		ar[double] ycoords = y
		ar[double] zcoords = z
		ar[double] hsml = h
		ar[double] uvals = u
		double msink = central_point_mass

 	# pointers?
 	# figure out basic C syntax or whatever this is
    c_gfunc(&ngas, &mgas,\
    		<double*> xcoords, <double*> ycoords, <double*> zcoords,\
    		<double*> hsml, <double*> uvals, &msink)
    #c_gfunc(&x, &n, &n, <double*> ax.data, <double*> ax.data, <double*> c.data)
    
    return 