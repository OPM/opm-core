/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */

#ifndef SPU_IMPLICIT_H_INCLUDED
#define SPU_IMPLICIT_H_INCLUDED
double
spu_implicit(grid_t *g, double *s0, double *s, double *mob, double *dmob,
             double *dflux, double *gflux, double *src, double dt);
double 
spu_implicit2(grid_t *g, double *s0, double *s, double *mob, double *dmob,
              double *dflux, double *gflux, double *src, double dt);

void 
compute_mobilities(int n, double *s, double *mob, double *dmob, 
                   int ntab, double h, double x0, double *tab);
#endif /* SPU_IMPLICIT_H_INCLUDED */

