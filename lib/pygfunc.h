/* ngas is an int, the rest are arrays, msink is a float

	is double* the correct prefix for arrays of doubles AND for stand-alone doubles?

extern void c_gfunc(int ngas, double* msink);

/* extremely careful with asterisks (pointers) here !!! */
extern void c_gfunc(int* ngas, double* mgas, double* x, double* y, double* z, double* h, double* u, double* msink,double* hsink_soft);
/*
, double* x, double* y, double* z, double* h, double* u, double* msink);

 
whatever the arguments in pygfunc.f90 are, in here
ex) extern void c_gfunc(double* x, int* n, int* m, double* a, double* b, double* c);
*/
