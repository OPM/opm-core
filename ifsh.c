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

    double        *S;           /* Schur complement contrib per cell */
    double        *gpress;
    double        *cflux;       /* Cell (half-contact) fluxes */

    int           *pgconn;
    int           *gconn;

    int           *cwpos;
    int           *cwells;

    struct hybsys *sys;

    /* Linear storage goes here... */
    int     idata_default;
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

        idata_sz = 0;

        nc         = G->number_of_cells;
        ngconn_tot = G->cell_facepos[nc];

        nnu = G->number_of_faces;
        if (W != NULL) {
            nnu += W->number_of_wells;

            idata_sz  = nc + 1;
            idata_sz += 2 * W->well_connpos[ W->number_of_wells ];
        }

        ddata_sz  = 2 * nnu;             /* rhs + soln */
        ddata_sz += new->sum_ngconn2;    /* Binv */
        ddata_sz += ngconn_tot;          /* gpress */
        ddata_sz += ngconn_tot;          /* cflux */

        if (idata_sz == 0) {
            /* malloc(0) may or may not return NULL.  Handle. */
            new->pimpl->idata = &new->pimpl->idata_default;
        } else {
            new->pimpl->idata = malloc(idata_sz * sizeof *new->pimpl->idata);
        }

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
            new->pimpl->S      = new->x             + nnu;
            new->pimpl->gpress = new->pimpl->S      + new->sum_ngconn2;
            new->pimpl->cflux  = new->pimpl->gpress + ngconn_tot;
        }
    }

    if ((new != NULL) &&
        (new->pimpl->idata != &new->pimpl->idata_default)) {
        new->pimpl->cwpos  = new->pimpl->idata;
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

        if (pimpl->idata != &pimpl->idata_default) {
            free   (pimpl->idata);
        }
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
    int    *pgconn;
    double *S, *gp;

    nc     = ifsh->pimpl->nc;
    S      = ifsh->pimpl->S;
    gp     = ifsh->pimpl->gpress;
    pgconn = ifsh->pimpl->pgconn;

    p2 = 0;
    for (c = 0; c < ifsh->pimpl->nc; c++) {
        n = pgconn[c + 1] - pgconn[c];

        for (i = 0; i < n    ; i++) {
           gp[p1 + i] = omega[c] * gpress[p1 + i];
        }

        for (i = 0; i < n * n; i++) {
            S [p2 + i] = totmob[c] * Binv[p2 + i];
        }

        p1 += n    ;
        p2 += n * n;
    }
}
