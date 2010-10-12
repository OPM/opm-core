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

#include "ifsh.h"
#include "hybsys.h"
#include "hybsys_global.h"

#if defined MAX
#undef MAX
#endif
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

struct ifsh_impl {
    int nc, nf, nw;             /* Number of cells, faces, wells */

    double        *Binv;        /* Effective inverse IP */
    double        *gpress;      /* Effective gravity pressure */
    double        *cflux;       /* Cell (half-contact) fluxes */

    double        *WI;          /* Effective inverse IP, wells */
    double        *wdp;         /* Effective grav. press, wells */

    int           *pgconn;      /* Start pointers, grid connections */
    int           *gconn;       /* Grid connections */

    int           *cwpos;       /* Start pointers, c->w mapping */
    int           *cwells;      /* c->w mapping */

    struct hybsys      *sys;    /* Hybrid cell contribs */
    struct hybsys_well *wsys;   /* Hybrid cell contribs from wells */

    double *work;               /* Scratch array, floating point */
    int    *iwork;              /* Scratch array, integers */

    /* Linear storage goes here... */
    int    *idata;              /* Actual storage array, integers */
    double *ddata;              /* Actual storage array, floating point */
};


/* ---------------------------------------------------------------------- */
/* Determine nnz (=sum(diff(facePos)^2)) and max(diff(facePos) for grid */
/* ---------------------------------------------------------------------- */
static void
count_grid_connections(grid_t *G, int *max_ngconn, size_t *sum_ngconn2)
/* ---------------------------------------------------------------------- */
{
    int c, n;

    *max_ngconn  = INT_MIN;
    *sum_ngconn2 = 0;

    for (c = 0; c < G->number_of_cells; c++) {
        n = G->cell_facepos[c + 1] - G->cell_facepos[c];

        *max_ngconn   = MAX(*max_ngconn, n);
        *sum_ngconn2 += n * n;
    }
}


/* ---------------------------------------------------------------------- */
static void
ifsh_compute_table_sz(grid_t *G, well_t *W,
                      int max_ngconn, size_t sum_ngconn2,
                      size_t *nnu, size_t *idata_sz, size_t *ddata_sz)
/* ---------------------------------------------------------------------- */
{
    int nc, ngconn_tot;

    *nnu = G->number_of_faces;

    nc         = G->number_of_cells;
    ngconn_tot = G->cell_facepos[nc];

    *idata_sz  = max_ngconn;    /* iwork */

    *ddata_sz  = 2 * (*nnu);    /* rhs + soln */
    *ddata_sz += sum_ngconn2;   /* Binv */
    *ddata_sz += ngconn_tot;    /* gpress */
    *ddata_sz += ngconn_tot;    /* cflux */
    *ddata_sz += max_ngconn;    /* work */

    if (W != NULL) {
        *nnu += W->number_of_wells;

        /* cwpos */
        *idata_sz += nc + 1;

        /* cwells */
        *idata_sz += 2 * W->well_connpos[ W->number_of_wells ];

        /* rhs + soln */
        *ddata_sz += 2 * W->number_of_wells;

        /* WI, wdp */
        *ddata_sz += 2 * W->well_connpos[ W->number_of_wells ];
    }
}


/* ---------------------------------------------------------------------- */
/* Allocate and define supporting structures for assembling the global
 * system of linear equations to couple the grid (reservoir)
 * connections represented by 'G' and, if present (i.e., non-NULL),
 * the well connections represented by 'W'. */
/* ---------------------------------------------------------------------- */
struct ifsh_data *
ifsh_construct(grid_t *G, well_t *W)
/* ---------------------------------------------------------------------- */
{
    int               nc, ngconn_tot;
    size_t            idata_sz, ddata_sz, nnu;
    struct ifsh_data *new;

    assert (G != NULL);

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->A     = hybsys_define_globconn(G, W);
        new->pimpl = NULL;

        if (new->A == NULL) {
            ifsh_destroy(new);
            new = NULL;
        }
    }

    if (new != NULL) {
        new->pimpl = malloc(1 * sizeof *new->pimpl);

        if (new->pimpl == NULL) {
            ifsh_destroy(new);
            new = NULL;
        } else {
            new->pimpl->sys   = NULL;     new->pimpl->wsys  = NULL;
            new->pimpl->idata = NULL;     new->pimpl->ddata = NULL;
        }
    }

    if (new != NULL) {
        count_grid_connections(G, &new->max_ngconn, &new->sum_ngconn2);

        ifsh_compute_table_sz(G, W, new->max_ngconn, new->sum_ngconn2,
                              &nnu, &idata_sz, &ddata_sz);

        nc         = G->number_of_cells;
        ngconn_tot = G->cell_facepos[nc];

        new->pimpl->idata = malloc(idata_sz * sizeof *new->pimpl->idata);
        new->pimpl->ddata = malloc(ddata_sz * sizeof *new->pimpl->ddata);
        new->pimpl->sys   = hybsys_allocate_symm(new->max_ngconn,
                                                 nc, ngconn_tot);
        if ((new->pimpl->idata == NULL) ||
            (new->pimpl->ddata == NULL) || (new->pimpl->sys == NULL)) {
            ifsh_destroy(new);
            new = NULL;
        } else {
            new->b             = new->pimpl->ddata;
            new->x             = new->b             + nnu;
            new->pimpl->Binv   = new->x             + nnu;
            new->pimpl->gpress = new->pimpl->Binv   + new->sum_ngconn2;
            new->pimpl->cflux  = new->pimpl->gpress + ngconn_tot;
            new->pimpl->work   = new->pimpl->cflux  + ngconn_tot;

            new->pimpl->iwork  = new->pimpl->idata;

            hybsys_init(new->max_ngconn, new->pimpl->sys);
        }
    }

    if ((new != NULL) && (W != NULL)) {
        new->pimpl->cwpos  = new->pimpl->iwork + new->max_ngconn;
        new->pimpl->cwells = new->pimpl->cwpos + nc + 1;

        memset(new->pimpl->cwpos, 0, (nc + 1) * sizeof *new->pimpl->cwpos);

        derive_cell_wells(nc, W, new->pimpl->cwpos, new->pimpl->cwells);

        new->pimpl->wsys = hybsys_well_allocate_symm(new->max_ngconn, nc,
                                                     new->pimpl->cwpos);

        if (new->pimpl->wsys == NULL) {
            ifsh_destroy(new);
            new = NULL;
        } else {
            new->pimpl->WI  = new->pimpl->work + new->max_ngconn;
            new->pimpl->wdp = new->pimpl->WI + W->well_connpos[ W->number_of_wells ];
        }
    }

    if (new != NULL) {
        new->pimpl->nc = nc;
        new->pimpl->nf = G->number_of_faces;
        new->pimpl->nw = (W != NULL) ? W->number_of_wells : 0;

        /* Contentious.  Imposes severe restrictions on lifetime(G). */
        new->pimpl->pgconn = G->cell_facepos;
        new->pimpl->gconn  = G->cell_faces;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
/* Release dynamic memory resources associated to internal data of a
 * particular (incompressible) flow solver instance. */
/* ---------------------------------------------------------------------- */
static void
ifsh_destroy_impl(struct ifsh_impl *pimpl)
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
/* Release memory resources for ifsh data-handle 'h' */
/* ---------------------------------------------------------------------- */
void
ifsh_destroy(struct ifsh_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        ifsh_destroy_impl(h->pimpl);
        csrmatrix_delete (h->A    );
    }

    free(h);
}


/* ---------------------------------------------------------------------- */
/* Eliminate 'npp' prescribed (Dirichlet) conditions (locdof,dofval)
 * from global system of linear equations.  Move known values to RHS
 * whilst zeroing coefficient matrix contributions.  Local system of
 * dimension 'n'. */
/* ---------------------------------------------------------------------- */
static void
ifsh_explicit_elimination(int n, int npp,
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
/* Impose boundary conditions on local contribution to global system. */
/* ---------------------------------------------------------------------- */
static int
ifsh_impose_bc(int nconn, int *conn, flowbc_t *bc,
               struct ifsh_impl *pimpl)
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
        ifsh_explicit_elimination(nconn, npp,
                                  pimpl->iwork,
                                  pimpl->work,
                                  pimpl->sys->S,
                                  pimpl->sys->r);
    }

    return npp;
}


/* ---------------------------------------------------------------------- */
static void
ifsh_set_effective_params(const double     *Binv,
                          const double     *gpress,
                          const double     *totmob,
                          const double     *omega,
                          struct ifsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int     c, n, nc, p1, p2, i;
    int    *pgconn, *gconn;
    double *BI, *gp;

    nc     = ifsh->pimpl->nc;
    pgconn = ifsh->pimpl->pgconn;
    gconn  = ifsh->pimpl->gconn;

    BI     = ifsh->pimpl->Binv;
    gp     = ifsh->pimpl->gpress;

    p1 = p2 = 0;
    for (c = 0; c < nc; c++) {
        n = pgconn[c + 1] - pgconn[c];

        for (i = 0; i < n    ; i++) {
            gp[p1 + i] = omega[c] * gpress[p1 + i];
        }

        for (i = 0; i < n * n; i++) {
            BI[p2 + i] = totmob[c] * Binv[p2 + i];
        }

        p1 += n;
        p2 += n * n;
    }
}



/* ---------------------------------------------------------------------- */
static int
ifsh_assemble_grid(flowbc_t         *bc,
                   const double     *src,
                   struct ifsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int     c, n, nc, p1, p2;
    int     npp;
    int    *pgconn, *gconn;
    double *BI, *gp;

    nc     = ifsh->pimpl->nc;
    pgconn = ifsh->pimpl->pgconn;
    gconn  = ifsh->pimpl->gconn;

    BI     = ifsh->pimpl->Binv;
    gp     = ifsh->pimpl->gpress;

    p1 = p2 = npp = 0;
    for (c = 0; c < nc; c++) {
        n = pgconn[c + 1] - pgconn[c];

        hybsys_cellcontrib_symm(c, n, p1, p2, gp, src, BI,
                                ifsh->pimpl->sys);

        npp += ifsh_impose_bc(n, gconn + p1, bc, ifsh->pimpl);

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
ifsh_set_effective_well_params(const double     *WI,
                               const double     *wdp,
                               const double     *totmob,
                               const double     *omega,
                               struct ifsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int     c, nc, i, perf;
    int    *cwpos, *cwells;
    double *wsys_WI, *wsys_wdp;

    nc       = ifsh->pimpl->nc;
    cwpos    = ifsh->pimpl->cwpos;
    cwells   = ifsh->pimpl->cwells;
    wsys_WI  = ifsh->pimpl->WI;
    wsys_wdp = ifsh->pimpl->wdp;

    for (c = i = 0; c < nc; c++) {
        for (; i < cwpos[c + 1]; i++) {
            perf = cwells[2*i + 1];

            wsys_WI [i] = totmob[c] * WI [perf];
            wsys_wdp[i] = omega [c] * wdp[perf];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
ifsh_impose_well_control(int               c,
                         flowbc_t         *bc,
                         well_control_t   *wctrl,
                         struct ifsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int  ngconn, nwconn, i, w1, w2, wg, f;
    int *pgconn, *gconn, *pwconn, *wconn;

    double bhp;
    double *r, *r2w, *w2w;

    /* Enforce symmetric system */
    assert (ifsh->pimpl->wsys->r2w == ifsh->pimpl->wsys->w2r);

    pgconn = ifsh->pimpl->pgconn;
    pwconn = ifsh->pimpl->cwpos;

    gconn  = ifsh->pimpl->gconn  +   pgconn[c];
    wconn  = ifsh->pimpl->cwells + 2*pwconn[c];

    ngconn = pgconn[c + 1] - pgconn[c];
    nwconn = pwconn[c + 1] - pwconn[c];

    r2w = ifsh->pimpl->wsys->r2w;
    w2w = ifsh->pimpl->wsys->w2w;
    r   = ifsh->pimpl->wsys->r  ;

    /* Adapt local system to prescribed boundary pressures (r->w) */
    for (i = 0; i < ngconn; i++) {
        f = gconn[i];

        if (bc->type[f] == PRESSURE) {
            for (w1 = 0; w1 < nwconn; w1++) {
                /* Eliminate prescribed (boundary) pressure value */
                r  [ngconn + w1]   -= r2w[i + w1*ngconn] * bc->bcval[f];
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
                assert ((bc->type[gconn[i]] != PRESSURE) ||
                        !(fabs(r2w[i + w1*ngconn]) > 0.0));

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
ifsh_assemble_well(flowbc_t         *bc,
                   well_control_t   *wctrl,
                   struct ifsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int npp;
    int ngconn, nwconn, c, nc, w;

    int *pgconn, *gconn, *pwconn, *wconn;

    nc = ifsh->pimpl->nc;

    pgconn = ifsh->pimpl->pgconn;
    gconn  = ifsh->pimpl->gconn;
    pwconn = ifsh->pimpl->cwpos;
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


/* ---------------------------------------------------------------------- */
/* Assemble global system of linear equations
 *
 *     ifsh->A * ifhs->x = ifhs->b
 *
 * from local inner product matrices Binv, gravity pressure gpress,
 * boundary conditions bc, source terms src and fluid properties
 * totmob and omega. */
/* ---------------------------------------------------------------------- */
void
ifsh_assemble(flowbc_t         *bc,
              double           *src,
              double           *Binv,
              double           *gpress,
              well_control_t   *wctrl,
              double           *WI,
              double           *wdp,
              double           *totmob, /* \sum_i \lambda_i */
              double           *omega,  /* \sum_i \rho_i f_i */
              struct ifsh_data *ifsh)
/* ---------------------------------------------------------------------- */
{
    int npp;                /* Number of prescribed pressure values */

    ifsh_set_effective_params(Binv, gpress, totmob, omega, ifsh);

    hybsys_schur_comp_symm(ifsh->pimpl->nc,
                           ifsh->pimpl->pgconn,
                           ifsh->pimpl->Binv,
                           ifsh->pimpl->sys);

    if (ifsh->pimpl->nw > 0) {
        ifsh_set_effective_well_params(WI, wdp, totmob, omega, ifsh);

        hybsys_well_schur_comp_symm(ifsh->pimpl->nc,
                                    ifsh->pimpl->cwpos,
                                    ifsh->pimpl->WI,
                                    ifsh->pimpl->sys,
                                    ifsh->pimpl->wsys);
    }

    npp = ifsh_assemble_grid(bc, src, ifsh);

    if (ifsh->pimpl->nw > 0) {
        npp += ifsh_assemble_well(bc, wctrl, ifsh);
    }

    if (npp == 0) {
        ifsh->A->sa[0] *= 2;        /* Remove zero eigenvalue */
    }
}


/* ---------------------------------------------------------------------- */
/* Compute cell pressures (cpress) and interface fluxes (fflux) from
 * current solution of system of linear equations, h->x.  Back
 * substitution process, projected half-contact fluxes. */
/* ---------------------------------------------------------------------- */
void
ifsh_press_flux(grid_t *G, struct ifsh_data *h,
                double *cpress, double *fflux,
                double *wpress, double *wflux)
/* ---------------------------------------------------------------------- */
{
    int c, f, i;
    double s;

    hybsys_compute_press_flux(G->number_of_cells,
                              G->cell_facepos,
                              G->cell_faces,
                              h->pimpl->gpress,
                              h->pimpl->Binv,
                              h->pimpl->sys,
                              h->x, cpress, h->pimpl->cflux,
                              h->pimpl->work);

    if (h->pimpl->nw > 0) {
        assert ((wpress != NULL) && (wflux != NULL));
        hybsys_compute_press_flux_well(G->number_of_cells, G->cell_facepos,
                                       G->number_of_faces, h->pimpl->nw,
                                       h->pimpl->cwpos, h->pimpl->cwells,
                                       h->pimpl->Binv, h->pimpl->WI,
                                       h->pimpl->wdp, h->pimpl->sys,
                                       h->pimpl->wsys, h->x, cpress,
                                       h->pimpl->cflux, wpress, wflux,
                                       h->pimpl->work);
    }

    for (f = 0; f < G->number_of_faces; f++) { fflux[f] = 0.0; }

    i = 0;
    for (c = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];
            s = 2.0*(G->face_cells[2*f + 0] == c) - 1.0;

            fflux[f] += s * h->pimpl->cflux[i];
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        i = (G->face_cells[2*f + 0] >= 0) +
            (G->face_cells[2*f + 1] >= 0);

        fflux[f] /= i;
    }
}
