/*extern void c_gfunc(int* n);
*/
extern void c_gfunc(int ngas, int* mgas, double* x, double* y, double* z, double* h, double* u, double* msink);

/* whatever the arguments in pygfunc.f90 are, in here
ex) extern void c_gfunc(double* x, int* n, int* m, double* a, double* b, double* c);
*/
