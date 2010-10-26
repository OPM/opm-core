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

#ifndef OPM_COARSE_SYS_HEADER_INCLUDED
#define OPM_COARSE_SYS_HEADER_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

#include "grid.h"

/* ---------------------------------------------------------------------- */

struct coarse_sys {
    int *dof2conn;           /* Map dof->connection (coarse interface) */
    int *blkdof_pos;         /* Start pointers to each block's dofs */
    int *basis_pos;          /* Start pointers to each block's bf's */
    int *cell_ip_pos;        /* Start pointers to each block's IP */

    int    *blkdof;          /* Each block's dofs */
    double *basis;           /* All basis functions */
    double *cell_ip;         /* Fine-scale IP contributions */
    double *Binv;            /* Coarse-scale inverse IP per block */
};


/* ---------------------------------------------------------------------- */

struct coarse_topology;
struct CSRMatrix;

typedef void (*LocalSolver)(struct CSRMatrix *A,
                            double           *b,
                            double           *x);

struct coarse_sys *
coarse_sys_construct(grid_t *g, const int   *p,
                     struct coarse_topology *ct,
                     const double           *perm,
                     const double           *src,
                     const double           *totmob,
                     LocalSolver             linsolve);

void
coarse_sys_destroy(struct coarse_sys *sys);

void
coarse_sys_compute_cell_ip(int                nc,
                           int                max_nconn,
                           int                nb,
                           const int         *pconn,
                           const double      *Binv,
                           const int         *b2c_pos,
                           const int         *b2c,
                           struct coarse_sys *sys);

void
coarse_sys_compute_Binv(int                nb,
                        int                max_bcells,
                        const double      *totmob,
                        const int         *b2c_pos,
                        const int         *b2c,
                        struct coarse_sys *sys,
                        double            *work);


#ifdef __cplusplus
}
#endif

#endif  /* OPM_COARSE_SYS_HEADER_INCLUDED */
