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



void
spu_implicit_assemble(struct UnstructuredGrid *g, double *s0, double *s, double *mob, double *dmob,
                      double *dflux, double *gflux, double *src, double dt, sparse_t *S,
                      double *b);

double
spu_implicit(struct UnstructuredGrid *g, double *s0, double *s, double h, double x0, int ntab, double *tab,
              double *dflux, double *gflux, double *src, double dt,
             void (*linear_solver)(int, int*, int*, double *, double *, double *));

#endif /* SPU_IMPLICIT_H_INCLUDED */

