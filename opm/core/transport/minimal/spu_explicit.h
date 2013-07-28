/*
 * Copyright 2010 (c) SINTEF ICT, Applied Mathematics.
 * Jostein R. Natvig <Jostein.R.Natvig at sintef.no>
 */

#ifndef SPU_EXPLICIT_H_INCLUDED
#define SPU_EXPLICIT_H_INCLUDED
void
spu_explicit(struct UnstructuredGrid *g,
             double *s0,
             double *s,
             double *mob,
             double *dflux,
             double *gflux,
             double *src,
             double dt);

#endif /* SPU_EXPLICIT_H_INCLUDED */

