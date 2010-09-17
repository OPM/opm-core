#include <assert.h>
#include <limits.h>
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

    int           *pgconn;
    int           *gconn;

    int           *cwpos;
    int           *cwells;

    struct hybsys *sys;

    double *work;
    int    *iwork;

    /* Linear storage goes here... */
    int    *idata;
    double *ddata;
};


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
        new->pimpl = malloc(1 * sizeof *new->pimpl);

        if ((new->A == NULL) || (new->pimpl == NULL)) {
            ifsh_destroy(new);
            new = NULL;
        }
    }

    if (new != NULL) {
        count_grid_connections(G, &new->max_ngconn, &new->sum_ngconn2);

        idata_sz = new->max_ngconn;     /* iwork */

        nc         = G->number_of_cells;
        ngconn_tot = G->cell_facepos[nc];

        nnu = G->number_of_faces;
        if (W != NULL) {
            nnu += W->number_of_wells;

            /* cwpos and cwells */
            idata_sz  = nc + 1;
            idata_sz += 2 * W->well_connpos[ W->number_of_wells ];
        }

        ddata_sz  = 2 * nnu;             /* rhs + soln */
        ddata_sz += new->sum_ngconn2;    /* Binv */
        ddata_sz += ngconn_tot;          /* gpress */
        ddata_sz += ngconn_tot;          /* cflux */
        ddata_sz += new->max_ngconn;     /* work */

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
    }

    if (new != NULL) {
        new->pimpl->nc = nc;
        new->pimpl->nf = G->number_of_faces;
        new->pimpl->nw = (W != NULL) ? W->number_of_wells : 0;

        new->pimpl->pgconn = G->cell_facepos;
        new->pimpl->gconn  = G->cell_faces;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
ifsh_destroy_impl(struct ifsh_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        hybsys_free(pimpl->sys  );
        free       (pimpl->ddata);
        free       (pimpl->idata);
    }

    free(pimpl);
}


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
    int     c, n, nc, p1, p2, i;
    int     npp;                /* Number of prescribed pressure values */
    int    *pgconn, *gconn;
    double *BI, *gp;

    nc     = ifsh->pimpl->nc;
    BI     = ifsh->pimpl->Binv;
    gp     = ifsh->pimpl->gpress;
    pgconn = ifsh->pimpl->pgconn;
    gconn  = ifsh->pimpl->gconn;

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

    hybsys_schur_comp_symm(nc, pgconn, BI, ifsh->pimpl->sys);

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

    if (npp == 0) {
        ifsh->A->sa[0] *= 2;        /* Remove zero eigenvalue */
    }
}


/* ---------------------------------------------------------------------- */
void
ifsh_press_flux(grid_t *G, struct ifsh_data *h, double *src,
                double *cpress, double *fflux)
/* ---------------------------------------------------------------------- */
{
    int c, f, i;
    double s;

    hybsys_compute_press_flux(G->number_of_cells,
                              G->cell_facepos,
                              G->cell_faces,
                              h->pimpl->gpress,
                              src, h->pimpl->Binv,
                              h->pimpl->sys,
                              h->x, cpress, h->pimpl->cflux,
                              h->pimpl->work);

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
