#include <assert.h>
#include <stddef.h>
#include <stdlib.h>

#include "hash_set.h"
#include "hybsys_global.h"


#if defined MAX
#undef MAX
#endif
#define MAX(a,b) (((a) > (b)) ? (a) : (b))


/* ---------------------------------------------------------------------- */
static void
delete_cell_wells(int *cwpos, int *cwells)
/* ---------------------------------------------------------------------- */
{
    free(cwells);
    free(cwpos);
}


/* ---------------------------------------------------------------------- */
static int
build_cell_wells(size_t nc, well_t *W, int **cwpos, int **cwells)
/* ---------------------------------------------------------------------- */
{
    int  ret, w;
    int *p, *cw, *c, *connpos;

    size_t i;

    ret = 0;

    p  = calloc(nc + 1, sizeof *p);
    cw = NULL;

    if (p != NULL) {
        connpos = W->well_connpos;

        c = W->well_cells;
        for (w = 0; w < W->number_of_wells; w++) {
            for (; c != W->well_cells + connpos[w + 1]; c++) {
                p[*c + 1] += 1;
            }
        }

        for (i = 1; i <= nc; i++) {
            p[0] += p[i];
            p[i]  = p[0] - p[i];
        }

        cw = malloc(p[0] * sizeof *cw);

        if (cw != NULL) {
            ret  = p[0];        /* Success */

            p[0] = 0;
            c    = W->well_cells;
            for (w = 0; w < W->number_of_wells; w++) {
                for (; c != W->well_cells + connpos[w + 1]; c++) {
                    cw[ p[*c + 1] ++ ] = w;
                }
            }
        } else {
            free(p);  p = NULL;
        }
    }

    *cwpos  = p;
    *cwells = cw;

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
deallocate_well_dofset(size_t nw, struct hash_set **wia)
/* ---------------------------------------------------------------------- */
{
    size_t w;

    if (wia != NULL) {
        for (w = 0; w < nw; w++) {
            hash_set_deallocate(wia[w]);
        }
    }

    free(wia);
}


/* ---------------------------------------------------------------------- */
static struct hash_set **
allocate_well_dofset(grid_t *G, well_t *W)
/* ---------------------------------------------------------------------- */
{
    int w, i, c, ok;
    int set_size;
    struct hash_set **wia;

    wia = malloc(W->number_of_wells * sizeof *wia);

    if (wia != NULL) {
        ok = 1;
        for (w = 0; ok && (w < W->number_of_wells); w++) {
            set_size = 1;       /* Approximate expected set size */

            for (i = W->well_connpos[w]; i < W->well_connpos[w + 1]; i++) {
                c = W->well_cells[i];

                set_size += G->cell_facepos[c + 1] - G->cell_facepos[c];
            }

            wia[w] = hash_set_allocate(set_size);

            ok = wia[w] != NULL;
        }

        if (!ok) {
            deallocate_well_dofset(w, wia);
            wia = NULL;
        }
    }

    return wia;
}


/* ---------------------------------------------------------------------- */
static void
count_conn_per_row_grid(grid_t *G, struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    int    c, nc, *ia, *ja;
    size_t m;

    nc = G->number_of_cells;
    ia = G->cell_facepos;
    ja = G->cell_faces;

    A->ia[0] = 0;
    for (m = 0; m < A->m; m++) {
        A->ia[m + 1] = 1;
    }

    for (c = 0; c < nc; c++) {
        for (; ja != G->cell_faces + ia[c + 1]; ja++) {
            A->ia[ *ja + 1 ] += ia[c + 1] - ia[c] - 1;
        }
    }
}


/* ---------------------------------------------------------------------- */
static int
count_conn_per_row_well(grid_t *G, well_t *W,
                        int              *cwpos,
                        int              *cwells,
                        struct hash_set  **wia,
                        struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    int c, w, nc, nf, nwdof, ok, i, dof, *ia, *ja, *wja;

    size_t m;

    nc  = G->number_of_cells;
    nf  = G->number_of_faces;
    ia  = G->cell_facepos;
    wja = cwells;

    ok = 1;
    for (c = 0; c < nc; c++) {
        for (; ok && (wja != cwells + cwpos[c + 1]); wja++) {
            for (       ja  = G->cell_faces + ia[c + 0];
                 ok && (ja != G->cell_faces + ia[c + 1]); ja++) {
                dof = *ja;
                ok  = hash_set_insert(dof, wia[*wja]) == dof;
            }

            for (i = cwpos[c]; ok && (i < cwpos[c + 1]); i++) {
                dof = nf + cwells[i];
                ok  = hash_set_insert(dof, wia[*wja]) == dof;
            }
        }
    }

    if (ok) {
        for (w = 0; w < W->number_of_wells; i++) {
            nwdof = 0;

            for (m = 0; m < wia[w]->m; m++) {
                if ((dof = wia[w]->s[m]) >= 0) {
                    A->ia[ dof + 1 ] += 1; /* face to well */
                    nwdof            += 1; /* well to face */
                }
            }

            A->ia[ nf + i + 1 ] = nwdof;
        }
    }

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
fill_self_connections(struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    size_t r;

    for (r = 0; r < A->m; r++) {
        A->ja[ A->ia[r + 1] ++ ] = r;
    }

    A->n = A->m;
}


/* ---------------------------------------------------------------------- */
static void
fill_grid_connections(grid_t *G, struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    int c, i, j, n;
    int dof1, dof2;

    int *ia, *ja;

    ia = G->cell_facepos;
    ja = G->cell_faces;

    for (c = 0; c < G->number_of_cells; c++) {
        n = ia[c + 1] - ia[c];

        for (i = 0; i < n; i++) {
            dof1 = ja[ ia[c] + i ];

            if (dof1 >= 0) {
                A->n = MAX(A->n, (size_t) dof1);
            }

            for (j = (i + 1) % n; j != i; j = (j + 1) % n) {
                dof2 = ja[ ia[c] + j ];

                A->ja[ A->ia[dof1 + 1] ++ ] = dof2;
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
fill_well_connections(int nf, int nw,
                      struct hash_set **wia,
                      struct CSRMatrix *A)
/* ---------------------------------------------------------------------- */
{
    int    w, dof;
    size_t i;

    for (w = 0; w < nw; w++) {
        for (i = 0; i < wia[w]->m; i++) {
            dof = wia[w]->s[i];

            if (dof >= 0) {
                A->n = MAX(A->n, (size_t)dof);

                if (dof < nf) {     /* Connect face to well */
                    A->ja[ A->ia[ dof + 1 ] ++ ] = nf + w;
                }

                if (dof != nf + w) {/* Connect well to dof (avoid self) */
                    A->ja[ A->ia[ nf + w + 1 ] ++ ] = dof;
                }
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
struct CSRMatrix *
hybsys_define_globconn(grid_t *G, well_t *W)
/* ---------------------------------------------------------------------- */
{
    int nw, ok;
    int *cwell_pos = NULL, *cwells = NULL;

    struct hash_set  **wia;
    struct CSRMatrix *A;

    assert (G != NULL);

    ok = 1;
    nw = 0;
    wia = NULL;
    if (W != NULL) {
        ok = build_cell_wells(G->number_of_cells, W, &cwell_pos, &cwells);
        nw = ok && W->number_of_wells;

        if (nw > 0) {
            wia = allocate_well_dofset(G, W);
            if (wia == NULL) { nw = 0; }
        }
    }

    A = csrmatrix_new_count_nnz(G->number_of_faces + nw);

    if (A != NULL) {
        count_conn_per_row_grid(G, A);

        if (nw > 0) {
            assert (nw == W->number_of_wells);
            ok = count_conn_per_row_well(G, W, cwell_pos,
                                         cwells, wia, A) > 0;
        } else {
            ok = 1;
        }

        ok = ok && (csrmatrix_new_elms_pushback(A) > 0);

        if (ok) {
            fill_self_connections(A);
            fill_grid_connections(G, A);
            fill_well_connections(G->number_of_faces, nw, wia, A);

            csrmatrix_sortrows(A);
        } else {
            csrmatrix_delete(A);
            A = NULL;
        }
    }

    deallocate_well_dofset(nw, wia);
    delete_cell_wells(cwell_pos, cwells);

    return A;
}
