/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */

#ifndef SPU_IMPLICIT_H_INCLUDED
#define SPU_IMPLICIT_H_INCLUDED

typedef struct Sparse {
    int     m;
    int     n;
    int    *ia;
    int    *ja;
    double *sa;
} sparse_t;



double interpolate(int n, double h, double x0, double *tab, double x);
double differentiate(int n, double h, double x0, double *tab, double x);
void 
spu_implicit_assemble(grid_t *g, double *s0, double *s, double *mob, double *dmob,
                      double *dflux, double *gflux, double *src, double dt, sparse_t *S, 
                      double *b);

double 
spu_implicit(grid_t *g, double *s0, double *s, double *mob, double *dmob,
              double *dflux, double *gflux, double *src, double dt, 
              void (*linear_solver)(int, int*, int*, double*, double*, double*));

void 
compute_mobilities(int n, double *s, double *mob, double *dmob, 
                   int ntab, double h, double x0, double *tab);
#endif /* SPU_IMPLICIT_H_INCLUDED */

