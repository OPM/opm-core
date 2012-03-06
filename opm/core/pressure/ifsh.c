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
static void
ifsh_set_effective_well_params(const double    *WI,
                               const double    *wdp,
                               struct fsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int     c, nc, i, perf;
    int    *cwpos, *cwells;
    double *wsys_WI, *wsys_wdp;

    nc       = ifsh->pimpl->nc;
    cwpos    = ifsh->pimpl->cwell_pos;
    cwells   = ifsh->pimpl->cwells;
    wsys_WI  = ifsh->pimpl->WI;
    wsys_wdp = ifsh->pimpl->wdp;

    for (c = i = 0; c < nc; c++) {
        for (; i < cwpos[c + 1]; i++) {
            perf = cwells[2*i + 1];

            wsys_WI [i] = WI [perf];
            wsys_wdp[i] = wdp[perf];
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
ifsh_assemble_grid(struct FlowBoundaryConditions *bc,
                   const double    *Binv,
                   const double    *gpress,
                   const double    *src,
                   struct fsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int     c, n, nc, p1, p2;
    int     npp;
    int    *pgconn, *gconn;

    nc     = ifsh->pimpl->nc;
    pgconn = ifsh->pimpl->gdof_pos;
    gconn  = ifsh->pimpl->gdof;

    p1 = p2 = npp = 0;
    for (c = 0; c < nc; c++) {
        n = pgconn[c + 1] - pgconn[c];

        hybsys_cellcontrib_symm(c, n, p1, p2, gpress, src, Binv,
                                ifsh->pimpl->sys);

        npp += fsh_impose_bc(n, gconn + p1, bc, ifsh->pimpl);

        hybsys_global_assemble_cell(n, gconn + p1,
                                    ifsh->pimpl->sys->S,
                                    ifsh->pimpl->sys->r,
                                    ifsh->A, ifsh->b);

        p1 += n;
        p2 += n * n;
    }

    return npp;
}


/* ---------------------------------------------------------------------- */
static void
ifsh_impose_well_control(int              c,
                         struct FlowBoundaryConditions *bc,
                         well_control_t  *wctrl,
                         struct fsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int  ngconn, nwconn, i, j, w1, w2, wg, f;
    int *pgconn, *gconn, *pwconn, *wconn;

    double bhp;
    double *r, *r2w, *w2w;

    /* Enforce symmetric system */
    assert (ifsh->pimpl->wsys->r2w == ifsh->pimpl->wsys->w2r);

    pgconn = ifsh->pimpl->gdof_pos;
    pwconn = ifsh->pimpl->cwell_pos;

    gconn  = ifsh->pimpl->gdof   +   pgconn[c];
    wconn  = ifsh->pimpl->cwells + 2*pwconn[c];

    ngconn = pgconn[c + 1] - pgconn[c];
    nwconn = pwconn[c + 1] - pwconn[c];

    r2w = ifsh->pimpl->wsys->r2w;
    w2w = ifsh->pimpl->wsys->w2w;
    r   = ifsh->pimpl->wsys->r  ;

    /* Adapt local system to prescribed boundary pressures (r->w) */
    for (i = 0; i < ngconn; i++) {
        f = gconn[i];
        j = ifsh->pimpl->bdry_condition[ f ];

        if (j != -1) {
            for (w1 = 0; w1 < nwconn; w1++) {
                /* Eliminate prescribed (boundary) pressure value */
                r  [ngconn + w1]   -= r2w[i + w1*ngconn] * bc->value[j];
                r2w[i + w1*ngconn]  = 0.0;
            }

            r[i] = 0.0;         /* RHS value handled in *reservoir* asm */
        }
    }

    /* Adapt local system to prescribed well (bottom-hole) pressures;
     * w->r and w->w. */
    for (w1 = 0; w1 < nwconn; w1++) {
        wg = wconn[2*w1 + 0];

        if (wctrl->ctrl[wg] == BHP) {
            bhp = wctrl->target[wg];

            /* Well->reservoir */
            for (i = 0; i < ngconn; i++) {
#ifndef NDEBUG
                j = ifsh->pimpl->bdry_condition[ gconn[i] ];

                assert ((j == -1)                    ||
                        (bc->type[j] != BC_PRESSURE) ||
                        !(fabs(r2w[i + w1*ngconn]) > 0.0));
#endif

                r  [i]             -= r2w[i + w1*ngconn] * bhp;
                r2w[i + w1*ngconn]  = 0.0;
            }

            /* Well->well */
            for (w2 = (w1 + 1) % nwconn; w2 != w1; w2 = (w2 + 1) % nwconn) {
                r  [ngconn + w2]    -= w2w[w2 + w1*nwconn] * bhp;
                w2w[w2 + w1*ngconn]  = 0.0;
                w2w[w1 + w2*ngconn]  = 0.0;
            }

            /* Assemble final well equation of the form S*p_bh = S*p_bh^0 */
            assert (fabs(w2w[w1 * (nwconn + 1)]) > 0.0);

            r[ngconn + w1] = w2w[w1 * (nwconn + 1)] * bhp;
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
ifsh_assemble_well(struct FlowBoundaryConditions *bc,
                   well_control_t  *wctrl,
                   struct fsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int npp;
    int ngconn, nwconn, c, nc, w;

    int *pgconn, *gconn, *pwconn, *wconn;

    nc = ifsh->pimpl->nc;

    pgconn = ifsh->pimpl->gdof_pos;
    gconn  = ifsh->pimpl->gdof;
    pwconn = ifsh->pimpl->cwell_pos;
    wconn  = ifsh->pimpl->cwells;

    for (c = 0; c < nc; c++) {
        ngconn = pgconn[c + 1] - pgconn[c];
        nwconn = pwconn[c + 1] - pwconn[c];

        if (nwconn > 0) {
            hybsys_well_cellcontrib_symm(c, ngconn, pgconn[c],
                                         pwconn,
                                         ifsh->pimpl->WI,
                                         ifsh->pimpl->wdp,
                                         ifsh->pimpl->sys,
                                         ifsh->pimpl->wsys);

            ifsh_impose_well_control(c, bc, wctrl, ifsh);

            hybsys_global_assemble_well_sym(ifsh->pimpl->nf,
                                            ngconn, gconn +   pgconn[c],
                                            nwconn, wconn + 2*pwconn[c] + 0,
                                            ifsh->pimpl->wsys->r2w,
                                            ifsh->pimpl->wsys->w2w,
                                            ifsh->pimpl->wsys->r,
                                            ifsh->A, ifsh->b);
        }
    }

    npp = 0;
    for (w = 0; w < ifsh->pimpl->nw; w++) {
        if (wctrl->ctrl[w] == BHP) {
            npp += 1;
        } else if (wctrl->ctrl[w] == RATE) {
            /* Impose total rate constraint.
             *
             * Note sign resulting from ->target[w] denoting
             * *injection* flux. */
            ifsh->b[ifsh->pimpl->nf + w] -= - wctrl->target[w];
        }
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
ifsh_construct(struct UnstructuredGrid *G, well_t *W)
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


    /* Allocate Schur complement contributions.  Symmetric system. */
    if (new != NULL) {
        nc         = G->number_of_cells;
        ngconn_tot = G->cell_facepos[nc];

        fsh_define_linsys_arrays(new);
        fsh_define_impl_arrays(nc, G->number_of_faces, nnu,
                               ngconn_tot, new->max_ngconn, W, new->pimpl);

        new->pimpl->sys = hybsys_allocate_symm(new->max_ngconn,
                                               nc, ngconn_tot);

        if (W != NULL) {
            fsh_define_cell_wells(nc, W, new->pimpl);

            new->pimpl->wsys = hybsys_well_allocate_symm(new->max_ngconn, nc,
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
 *     ifsh->A * ifsh->x = ifsh->b
 *
 * from effective local inner product matrices Binv, effective gravity
 * pressure gpress, boundary conditions bc, and source terms src. */
/* ---------------------------------------------------------------------- */
void
ifsh_assemble(struct FlowBoundaryConditions *bc,
              const double     *src,
              const double     *Binv,
              const double     *gpress,
              well_control_t   *wctrl,
              const double     *WI,
              const double     *wdp,
              struct fsh_data  *ifsh)
/* ---------------------------------------------------------------------- */
{
    int npp;                /* Number of prescribed pressure values */

    fsh_map_bdry_condition(bc, ifsh->pimpl);

    hybsys_schur_comp_symm(ifsh->pimpl->nc,
                           ifsh->pimpl->gdof_pos,
                           Binv, ifsh->pimpl->sys);

    if (ifsh->pimpl->nw > 0) {
        ifsh_set_effective_well_params(WI, wdp, ifsh);

        hybsys_well_schur_comp_symm(ifsh->pimpl->nc,
                                    ifsh->pimpl->cwell_pos,
                                    ifsh->pimpl->WI,
                                    ifsh->pimpl->sys,
                                    ifsh->pimpl->wsys);
    }

    npp = ifsh_assemble_grid(bc, Binv, gpress, src, ifsh);

    if (ifsh->pimpl->nw > 0) {
        npp += ifsh_assemble_well(bc, wctrl, ifsh);
    }

    if (npp == 0) {
        ifsh->A->sa[0] *= 2;        /* Remove zero eigenvalue */
    }
}
