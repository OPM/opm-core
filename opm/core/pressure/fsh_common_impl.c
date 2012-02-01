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

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/grid.h>
#include <opm/core/well.h>
#include <opm/core/pressure/flow_bc.h>

#include <opm/core/pressure/fsh.h>
#include <opm/core/pressure/fsh_common_impl.h>

#include <opm/core/pressure/mimetic/hybsys.h>
#include <opm/core/pressure/mimetic/hybsys_global.h>

#if defined MAX
#undef MAX
#endif
#define MAX(a,b) (((a) > (b)) ? (a) : (b))



/* ---------------------------------------------------------------------- */
/* Release dynamic memory resources associated to internal data of a
 * particular (incompressible) flow solver instance. */
/* ---------------------------------------------------------------------- */
static void
fsh_destroy_impl(struct fsh_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        hybsys_well_free(pimpl->wsys );
        hybsys_free     (pimpl->sys  );
        free            (pimpl->ddata);
        free            (pimpl->idata);
    }

    free(pimpl);
}


/* ---------------------------------------------------------------------- */
/* Eliminate 'npp' prescribed (Dirichlet) conditions (locdof,dofval)
 * from global system of linear equations.  Move known values to RHS
 * whilst zeroing coefficient matrix contributions.  Local system of
 * dimension 'n'. */
/* ---------------------------------------------------------------------- */
static void
fsh_explicit_elimination(int n, int npp,
                         int *locdof, double *dofval,
                         double *S, double *r)
/* ---------------------------------------------------------------------- */
{
    int i, p;

    for (i = 0; i < npp; i++) {
        for (p = (locdof[i] + 1) % n; p != locdof[i]; p = (p + 1) % n) {
            /* Subtract known values from RHS. */
            r[p] -= S[p + locdof[i]*n] * dofval[i];

            /* Eliminate DOF (zero row/col "locdof[i]"; leave diagonal). */
            S[p + locdof[i]*n] = 0.0;
            S[locdof[i] + p*n] = 0.0;
        }

        /* Produce trivial equation S(i,i)*x(i) = S(i,i)*x0(i). */
        r[p] = S[p * (n + 1)] * dofval[i];
    }
}


/* ---------------------------------------------------------------------- */
/* Release memory resources for ifsh data-handle 'h' */
/* ---------------------------------------------------------------------- */
void
fsh_destroy(struct fsh_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        fsh_destroy_impl(h->pimpl);
        csrmatrix_delete(h->A    );
    }

    free(h);
}


/* ---------------------------------------------------------------------- */
struct fsh_impl *
fsh_impl_allocate_basic(size_t idata_sz, size_t ddata_sz)
/* ---------------------------------------------------------------------- */
{
    struct fsh_impl *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->idata = malloc(idata_sz * sizeof *new->idata);
        new->ddata = malloc(ddata_sz * sizeof *new->ddata);

        new->sys   = NULL;
        new->wsys  = NULL;

        if ((new->idata == NULL) || (new->ddata == NULL)) {
            fsh_destroy_impl(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
/* Determine nnz (=sum(diff(facePos)^2)) and max(diff(facePos) for grid */
/* ---------------------------------------------------------------------- */
void
fsh_count_grid_dof(struct UnstructuredGrid *G, int *max_ngdof, size_t *sum_ngdof2)
/* ---------------------------------------------------------------------- */
{
    int c, n;

    *max_ngdof  = INT_MIN;
    *sum_ngdof2 = 0;

    for (c = 0; c < G->number_of_cells; c++) {
        n = G->cell_facepos[c + 1] - G->cell_facepos[c];

        *max_ngdof   = MAX(*max_ngdof, n);
        *sum_ngdof2 += n * n;
    }
}


/* ---------------------------------------------------------------------- */
/* Impose boundary conditions on local contribution to global system. */
/* ---------------------------------------------------------------------- */
int
fsh_impose_bc(int nconn, int *conn, flowbc_t *bc,
              struct fsh_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    int i, npp, f;

    npp = 0;
    for (i = 0; i < nconn; i++) {
        f = conn[i];

        if (bc->type[f] == PRESSURE) {
            pimpl->work [npp] = bc->bcval[f];
            pimpl->iwork[npp] = i;

            npp += 1;
        } else if (bc->type[f] == FLUX) {
            pimpl->sys->r[i] -= bc->bcval[f];
        }
    }

    if (npp > 0) {
        fsh_explicit_elimination(nconn, npp,
                                 pimpl->iwork,
                                 pimpl->work,
                                 pimpl->sys->S,
                                 pimpl->sys->r);
    }

    return npp;
}


/* ---------------------------------------------------------------------- */
void
fsh_define_impl_arrays(size_t           nc,
                       size_t           nnu,
                       size_t           nhf,
                       size_t           max_ncf,
                       well_t          *W,
                       struct fsh_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    pimpl->cflux    = pimpl->ddata  + 2 * nnu;
    pimpl->work     = pimpl->cflux  + nhf;

    pimpl->gdof_pos = pimpl->idata;
    pimpl->gdof     = pimpl->gdof_pos + (nc + 1);
    pimpl->iwork    = pimpl->gdof     + nhf;

    if (W != NULL) {
        pimpl->cwell_pos = pimpl->iwork     + max_ncf;
        pimpl->cwells    = pimpl->cwell_pos + nc + 1;

        pimpl->WI  = pimpl->work + max_ncf;
        pimpl->wdp = pimpl->WI + W->well_connpos[ W->number_of_wells ];
    }
}


/* ---------------------------------------------------------------------- */
void
fsh_define_cell_wells(size_t nc, well_t *W, struct fsh_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < nc + 1; i++) {
        pimpl->cwell_pos[i] = 0;
    }

    derive_cell_wells(nc, W, pimpl->cwell_pos, pimpl->cwells);
}


/* ---------------------------------------------------------------------- */
void
fsh_define_linsys_arrays(struct fsh_data *h)
/* ---------------------------------------------------------------------- */
{
    h->b = h->pimpl->ddata;
    h->x = h->b + h->A->m;
}


/* ---------------------------------------------------------------------- */
void
fsh_compute_table_sz(struct UnstructuredGrid *G, well_t *W, int max_ngconn,
                     size_t *nnu, size_t *idata_sz, size_t *ddata_sz)
/* ---------------------------------------------------------------------- */
{
    int nc, ngconn_tot;

    *nnu = G->number_of_faces;

    nc         = G->number_of_cells;
    ngconn_tot = G->cell_facepos[nc];

    *idata_sz  = nc + 1;        /* gdof_pos */
    *idata_sz += ngconn_tot;    /* gdof */
    *idata_sz += max_ngconn;    /* iwork */

    *ddata_sz  = 2 * (*nnu);    /* rhs + soln */
    *ddata_sz += ngconn_tot;    /* cflux */
    *ddata_sz += max_ngconn;    /* work */

    if (W != NULL) {
        *nnu += W->number_of_wells;

        /* cwell_pos */
        *idata_sz += nc + 1;

        /* cwells */
        *idata_sz += 2 * W->well_connpos[ W->number_of_wells ];

        /* rhs + soln */
        *ddata_sz += 2 * W->number_of_wells;

        /* WI, wdp */
        *ddata_sz += 2 * W->well_connpos[ W->number_of_wells ];
    }
}
