/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#ifndef IFSH_MS_H_INCLUDED
#define IFSH_MS_H_INCLUDED

#include <stddef.h>

#include "grid.h"
#include "coarse_sys.h"

struct CSRMatrix;
struct ifsh_ms_impl;

struct ifsh_ms_data {
    /* Linear system */
    struct CSRMatrix    *A;     /* Coefficient matrix */
    double              *b;     /* System RHS */
    double              *x;     /* Solution */

    /* Private implementational details. */
    struct ifsh_ms_impl *pimpl;
};


struct ifsh_ms_data *
ifsh_ms_construct(grid_t       *G,
                  const int    *p,
                  const double *perm,
                  const double *src,
                  const double *totmob,
                  LocalSolver   linsolve);

void
ifsh_ms_destroy(struct ifsh_ms_data *h);

void
ifsh_ms_assemble(const double        *src,
                 const double        *totmob,
                 struct ifsh_ms_data *h);

void
ifsh_ms_press_flux(grid_t *G, struct ifsh_ms_data *h,
                   double *cpress, double *fflux);

#endif /* IFSH_MS_H_INCLUDED */
