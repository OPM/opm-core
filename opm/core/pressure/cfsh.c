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

#include <opm/core/pressure/fsh.h>
#include <opm/core/pressure/fsh_common_impl.h>
#include <opm/core/pressure/mimetic/hybsys.h>
#include <opm/core/pressure/mimetic/hybsys_global.h>




/* ---------------------------------------------------------------------- */
static int
cfsh_assemble_grid(flowbc_t        *bc,
                   const double    *Binv,
                   const double    *gpress,
                   const double    *src,
                   struct fsh_data *h)
/* ---------------------------------------------------------------------- */
{
    int     c, n, nc, p1, p2;
    int     npp;
    int    *pgconn, *gconn;

    nc     = h->pimpl->nc;
    pgconn = h->pimpl->gdof_pos;
    gconn  = h->pimpl->gdof;

    p1 = p2 = npp = 0;
    for (c = 0; c < nc; c++) {
        n = pgconn[c + 1] - pgconn[c];

        hybsys_cellcontrib_unsymm(c, n, p1, p2, gpress, src, Binv,
                                  h->pimpl->sys);

        npp += fsh_impose_bc(n, gconn + p1, bc, h->pimpl);

        hybsys_global_assemble_cell(n, gconn + p1,
                                    h->pimpl->sys->S,
                                    h->pimpl->sys->r, h->A, h->b);

        p1 += n;
        p2 += n * n;
    }

    return npp;
}


/* ======================================================================
 * Public routines follow.
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
/* Allocate and define supporting structures for assembling the global
 * system of linear equations to couple the grid (reservoir)
 * connections represented by 'G' and, if present (i.e., non-NULL),
 * the well connections represented by 'W'. */
/* ---------------------------------------------------------------------- */
struct fsh_data *
cfsh_construct(struct UnstructuredGrid *G, well_t *W)
/* ---------------------------------------------------------------------- */
{
    int              nc, ngconn_tot;
    size_t           idata_sz, ddata_sz, nnu;
    struct fsh_data *new;

    assert (G != NULL);

    /* Allocate master structure, define system matrix sparsity */
    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->A     = hybsys_define_globconn(G, W);
        new->pimpl = NULL;

        if (new->A == NULL) {
            fsh_destroy(new);
            new = NULL;
        }
    }


    /* Allocate implementation structure */
    if (new != NULL) {
        fsh_count_grid_dof(G, &new->max_ngconn, &new->sum_ngconn2);

        fsh_compute_table_sz(G, W, new->max_ngconn,
                             &nnu, &idata_sz, &ddata_sz);

        new->pimpl = fsh_impl_allocate_basic(idata_sz, ddata_sz);

        if (new->pimpl == NULL) {
            fsh_destroy(new);
            new = NULL;
        }
    }


    /* Allocate Schur complement contributions.  Unsymmetric system. */
    if (new != NULL) {
        nc         = G->number_of_cells;
        ngconn_tot = G->cell_facepos[nc];

        fsh_define_linsys_arrays(new);
        fsh_define_impl_arrays(nc, nnu, ngconn_tot, new->max_ngconn,
                               W, new->pimpl);

        new->pimpl->sys = hybsys_allocate_unsymm(new->max_ngconn,
                                                 nc, ngconn_tot);

        if (W != NULL) {
            fsh_define_cell_wells(nc, W, new->pimpl);

            new->pimpl->wsys =
                hybsys_well_allocate_unsymm(new->max_ngconn, nc,
                                            new->pimpl->cwell_pos);
        }

        if ((new->pimpl->sys == NULL) ||
            ((W != NULL) && (new->pimpl->wsys == NULL))) {
            /* Failed to allocate ->sys or ->wsys (if W != NULL) */
            fsh_destroy(new);
            new = NULL;
        }
    }


    if (new != NULL) {
        /* All allocations succeded.  Fill metadata and return. */
        new->pimpl->nc = nc;
        new->pimpl->nf = G->number_of_faces;
        new->pimpl->nw = (W != NULL) ? W->number_of_wells : 0;

        memcpy(new->pimpl->gdof_pos,
               G->cell_facepos     ,
               (nc + 1)   * sizeof *new->pimpl->gdof_pos);

        memcpy(new->pimpl->gdof    ,
               G->cell_faces       ,
               ngconn_tot * sizeof *new->pimpl->gdof);

        hybsys_init(new->max_ngconn, new->pimpl->sys);
    }

    return new;
}


/* ---------------------------------------------------------------------- */
/* Assemble global system of linear equations
 *
 *     fsh->A * fsh->x = fsh->b
 */
/* ---------------------------------------------------------------------- */
void
cfsh_assemble(flowbc_t        *bc,
              const double    *src,
              const double    *Binv,
              const double    *Biv,
              const double    *P,
              const double    *gpress,
              well_control_t  *wctrl,
              const double    *WI,
              const double    *BivW,
              const double    *wdp,
              struct fsh_data *h)
/* ---------------------------------------------------------------------- */
{
    int npp;                /* Number of prescribed pressure values */

    /* Suppress warnings about unused parameters. */
    (void) wctrl;  (void) WI;  (void) BivW;  (void) wdp;

    hybsys_schur_comp_unsymm(h->pimpl->nc,
                             h->pimpl->gdof_pos,
                             Binv, Biv, P, h->pimpl->sys);

    npp = cfsh_assemble_grid(bc, Binv, gpress, src, h);

    if (npp == 0) {
        h->A->sa[0] *= 2;        /* Remove zero eigenvalue */
    }
}
