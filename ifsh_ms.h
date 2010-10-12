/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IFSH_MS_H_INCLUDED
#define IFSH_MS_H_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

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


#ifdef __cplusplus
}
#endif

#endif /* IFSH_MS_H_INCLUDED */
