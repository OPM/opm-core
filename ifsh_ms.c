/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "coarse_conn.h"
#include "coarse_sys.h"
#include "hybsys.h"
#include "hybsys_global.h"
#include "ifsh_ms.h"
#include "partition.h"
#include "sparse_sys.h"


#if defined(MAX)
#undef MAX
#endif
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))


struct ifsh_ms_impl {
    int    max_bcells;        /* Maximum number of cells per block */
    size_t ntotdof;           /* Total number of degrees of freedom */
    size_t max_nblkdof;       /* Max_i ndof(i) */
    size_t sum_nblkdof2;      /* Sum_i ndof(i)^2 */

    double *gpress;           /* Coarse gravity contributions */
    double *bsrc;             /* Coarse-scale source terms */

    double *bpress;           /* Block pressure */
    double *hflux;            /* Coarse-scale half-contact fluxes */
    double *flux;             /* Coarse-scale interface fluxes */

    double *work;             /* Assembly and back subst. work array */

    double *fs_hflux;         /* Fine-scale half-contact fluxes */

    int    *p;                /* Copy of partition vector */
    int    *pb2c, *b2c;       /* Bloc->cell mapping (inverse partition) */

    struct hybsys          *hsys;
    struct coarse_topology *ct;
    struct coarse_sys      *sys;

    /* Linear storage */
    double *ddata;
    int    *idata;
};


/* ---------------------------------------------------------------------- */
static void
ifsh_ms_impl_destroy(struct ifsh_ms_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    if (pimpl != NULL) {
        free                   (pimpl->idata);
        free                   (pimpl->ddata);
        hybsys_free            (pimpl->hsys );
        coarse_sys_destroy     (pimpl->sys  );
        coarse_topology_destroy(pimpl->ct   );
    }

    free(pimpl);
}


/* ---------------------------------------------------------------------- */
static int
max_block_cells(size_t nc, size_t nb, const int *p)
/* ---------------------------------------------------------------------- */
{
    int    ret, *cnt;
    size_t c;

    ret = -1;

    cnt = calloc(nb, sizeof *cnt);

    if (cnt != NULL) {
        ret = 0;

        for (c = 0; c < nc; c++) {
            cnt[p[c]] += 1;
            ret = MAX(ret, cnt[p[c]]);
        }
    }

    free(cnt);

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
count_blkdof(struct ifsh_ms_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    int b, ndof, p;

    pimpl->ntotdof      = 0;
    pimpl->max_nblkdof  = 0;
    pimpl->sum_nblkdof2 = 0;

    for (b = p = 0; b < pimpl->ct->nblocks; b++) {
        ndof = pimpl->sys->blkdof_pos[b + 1] - pimpl->sys->blkdof_pos[b];

        for (; p < pimpl->sys->blkdof_pos[b + 1]; p++) {
            pimpl->ntotdof = MAX(pimpl->ntotdof,
                                 (size_t) pimpl->sys->blkdof[p]);
        }

        pimpl->max_nblkdof   = MAX(pimpl->max_nblkdof, (size_t) ndof);
        pimpl->sum_nblkdof2 += ndof * ndof;
    }

    pimpl->ntotdof += 1;        /* Account for zero-based numbering */
}


/* ---------------------------------------------------------------------- */
static int
ifsh_ms_vectors_construct(grid_t *G, struct ifsh_ms_impl *pimpl)
/* ---------------------------------------------------------------------- */
{
    size_t nb, n_fs_gconn;
    size_t nblkdof_tot, max_pairs, work_sz;
    size_t alloc_sz;

    nb          = pimpl->ct->nblocks;
    nblkdof_tot = pimpl->sys->blkdof_pos[nb];
    max_pairs   = pimpl->max_nblkdof * (pimpl->max_nblkdof + 1) / 2;
    n_fs_gconn  = G->cell_facepos[ G->number_of_cells ];

    work_sz     = pimpl->max_bcells + max_pairs;

    alloc_sz    = pimpl->ntotdof;     /* b */
    alloc_sz   += pimpl->ntotdof;     /* x */

    alloc_sz   += nblkdof_tot;        /* gpress */
    alloc_sz   += nb;                 /* bsrc */
    alloc_sz   += nb;                 /* bpress */
    alloc_sz   += nblkdof_tot;        /* hflux */
    alloc_sz   += pimpl->ntotdof;     /* flux */

    alloc_sz   += work_sz;            /* work */

    alloc_sz   += n_fs_gconn;         /* fs_hflux */

    pimpl->ddata = malloc(alloc_sz * sizeof *pimpl->ddata);

    alloc_sz  = G->number_of_cells; /* p */
    alloc_sz += nb + 1;             /* pb2c */
    alloc_sz += G->number_of_cells; /* b2c */

    pimpl->idata = malloc(alloc_sz * sizeof *pimpl->idata);

    return (pimpl->ddata != NULL) && (pimpl->idata != NULL);
}


/* ---------------------------------------------------------------------- */
static void
set_impl_pointers(grid_t *G, struct ifsh_ms_data *h)
/* ---------------------------------------------------------------------- */
{
    size_t nb, n_fs_gconn;
    size_t nblkdof_tot, max_pairs, work_sz;

    nb          = h->pimpl->ct->nblocks;
    nblkdof_tot = h->pimpl->sys->blkdof_pos[nb];
    max_pairs   = h->pimpl->max_nblkdof * (h->pimpl->max_nblkdof + 1) / 2;
    n_fs_gconn  = G->cell_facepos[ G->number_of_cells ];

    work_sz     = h->pimpl->max_bcells + max_pairs;

    /* Linear system */
    h->b               = h->pimpl->ddata;
    h->x               = h->b             + h->pimpl->ntotdof;

    /* Coarse-scale back-substitution results */
    h->pimpl->gpress   = h->x             + h->pimpl->ntotdof;
    h->pimpl->bsrc     = h->pimpl->gpress + nblkdof_tot;
    h->pimpl->bpress   = h->pimpl->bsrc   + nb;
    h->pimpl->hflux    = h->pimpl->bpress + nb;
    h->pimpl->flux     = h->pimpl->hflux  + nblkdof_tot;

    /* Back-substitution work array */
    h->pimpl->work     = h->pimpl->flux   + h->pimpl->ntotdof;

    /* Fine-scale hflux accumulation array */
    h->pimpl->fs_hflux = h->pimpl->work   + work_sz;

    /* Partition vector */
    h->pimpl->p = h->pimpl->idata;

    /* Block->cell mapping (inverse partition vector) */
    h->pimpl->pb2c = h->pimpl->p    + G->number_of_cells;
    h->pimpl->b2c  = h->pimpl->pb2c + nb + 1;
}


/* ---------------------------------------------------------------------- */
static struct CSRMatrix *
ifsh_ms_matrix_construct(size_t m, size_t nnz, size_t nb,
                         const int *pdof, const int *dof)
/* ---------------------------------------------------------------------- */
{
    int               p, n, i1, dof1, i2, dof2;
    size_t            i, b;
    struct CSRMatrix *new;

    new = csrmatrix_new_known_nnz(m, nnz);

    if (new != NULL) {
        for (i = 0; i < m; i++) { new->ia[ i + 1 ] = 1; } /* Count self */

        /* Count others */
        for (b = 0, p = 0; b < nb; b++) {
            n = pdof[b + 1] - pdof[b];

            for (; p < pdof[b + 1]; p++) {
                new->ia[ dof[p] + 1 ] += n - 1;
            }
        }

        /* Set start pointers */
        new->ia[0] = 0;
        for (i = 1; i <= m; i++) {
            new->ia[0] += new->ia[i];
            new->ia[i]  = new->ia[0] - new->ia[i];
        }

        /* Insert self */
        for (i = 0; i < m; i++) {
            new->ja[ new->ia[ i + 1 ] ++ ] = (MAT_SIZE_T) i;
        }

        /* Insert others */
        for (b = 0; b < nb; b++) {
            p = pdof[b];
            n = pdof[b + 1] - p;

            for (i1 = 0; i1 < n; i1++) {
                dof1 = dof[p + i1];

                for (i2 = (i1 + 1) % n; i2 != i1; i2 = (i2 + 1) % n) {
                    dof2 = dof[p + i2];

                    new->ja[ new->ia[ dof1 + 1 ] ++ ] = dof2;
                }
            }
        }

        assert (new->ia[m] <= (MAT_SIZE_T) nnz);

        new->ia[0] = 0;
        csrmatrix_sortrows(new);
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static struct ifsh_ms_impl *
ifsh_ms_impl_construct(grid_t       *G     ,
                       const int    *p     ,
                       const double *perm  ,
                       const double *src   ,
                       const double *totmob,
                       LocalSolver   linsolve)
/* ---------------------------------------------------------------------- */
{
    int max_nconn, nb, nconn_tot;
    int expected_nconn, alloc_ok;

    struct ifsh_ms_impl *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->hsys  = NULL;      /* Crucial NULL-initialisation */
        new->sys   = NULL;
        new->ddata = NULL;
        new->idata = NULL;

        expected_nconn = 256;   /* CHANGE THIS! */
        alloc_ok       = 0;

        new->ct = coarse_topology_create(G->number_of_cells,
                                         G->number_of_faces,
                                         expected_nconn, p,
                                         G->face_cells);


        if (new->ct != NULL) {
            new->max_bcells = max_block_cells(G->number_of_cells,
                                              new->ct->nblocks, p);

            new->sys = coarse_sys_construct(G, p, new->ct, perm,
                                            src, totmob, linsolve);
        }

        if ((new->sys != NULL) && (new->max_bcells > 0)) {
            count_blkdof(new);

            max_nconn = new->max_nblkdof;
            nb        = new->ct->nblocks;
            nconn_tot = new->sys->blkdof_pos[ nb ];

            new->hsys = hybsys_allocate_symm(max_nconn, nb, nconn_tot);
        }

        if (new->hsys != NULL) {
            alloc_ok = ifsh_ms_vectors_construct(G, new);
        }

        if (! alloc_ok) {
            ifsh_ms_impl_destroy(new);
            new = NULL;
        } else {
            hybsys_init(max_nconn, new->hsys);
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
vector_zero(size_t m, double *v)
/* ---------------------------------------------------------------------- */
{
    size_t i;

    for (i = 0; i < m; i++) { v[i] = 0.0; }
}


/* ---------------------------------------------------------------------- */
static void
average_flux(size_t nf, const int *N, double *flux)          /* Poor name */
/* ---------------------------------------------------------------------- */
{
    size_t f;

    for (f = 0; f < nf; f++) {
        flux[f] /= ((N[2*f + 0] >= 0) + (N[2*f + 1] >= 0));
    }
}


/* ====================================================================== *
 * Public routines below.                                                 *
 * ====================================================================== */


/* ---------------------------------------------------------------------- */
struct ifsh_ms_data *
ifsh_ms_construct(grid_t       *G     ,
                  const int    *p     ,
                  const double *perm  ,
                  const double *src   ,
                  const double *totmob,
                  LocalSolver   linsolve)
/* ---------------------------------------------------------------------- */
{
    struct ifsh_ms_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        new->pimpl = ifsh_ms_impl_construct(G, p, perm, src,
                                            totmob, linsolve);
        new->A     = NULL;

        if (new->pimpl != NULL) {
            new->A = ifsh_ms_matrix_construct(new->pimpl->ntotdof,
                                              new->pimpl->sum_nblkdof2,
                                              new->pimpl->ct->nblocks,
                                              new->pimpl->sys->blkdof_pos,
                                              new->pimpl->sys->blkdof);
        }

        if (new->A == NULL) {
            ifsh_ms_destroy(new);
            new = NULL;
        } else {
            set_impl_pointers(G, new);

            memcpy(new->pimpl->p, p, G->number_of_cells * sizeof *p);
            memset(new->pimpl->pb2c, 0,
                   new->pimpl->ct->nblocks * sizeof *new->pimpl->pb2c);

            partition_invert(G->number_of_cells, new->pimpl->p,
                             new->pimpl->pb2c, new->pimpl->b2c);
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
ifsh_ms_destroy(struct ifsh_ms_data *h)
/* ---------------------------------------------------------------------- */
{
    if (h != NULL) {
        ifsh_ms_impl_destroy(h->pimpl);
        csrmatrix_delete    (h->A    );
    }

    free(h);
}


/* ---------------------------------------------------------------------- */
void
ifsh_ms_assemble(const double        *src   ,
                 const double        *totmob,
                 struct ifsh_ms_data *h)
/* ---------------------------------------------------------------------- */
{
    int b, i, nconn, p1, p2;

    int *pb2c, *b2c;
    double *bsrc;

    pb2c = h->pimpl->pb2c;
    b2c  = h->pimpl->b2c;
    bsrc = h->pimpl->bsrc;

    for (b = i = 0; b < h->pimpl->ct->nblocks; b++) {
        bsrc[b] = 0.0;

        for (; i < pb2c[b + 1]; i++) { bsrc[b] += src[ b2c[i] ]; }
    }

    coarse_sys_compute_Binv(h->pimpl->ct->nblocks,
                            h->pimpl->max_bcells, totmob,
                            pb2c, b2c, h->pimpl->sys, h->pimpl->work);

    hybsys_schur_comp_symm(h->pimpl->ct->nblocks,
                           h->pimpl->sys->blkdof_pos,
                           h->pimpl->sys->Binv, h->pimpl->hsys);

    csrmatrix_zero(         h->A);
    vector_zero   (h->A->m, h->b);

    /* Exclude gravity */
    vector_zero(h->pimpl->sys->blkdof_pos[ h->pimpl->ct->nblocks ],
                h->pimpl->gpress);

    p2 = 0;
    for (b = 0; b < h->pimpl->ct->nblocks; b++) {
        p1    = h->pimpl->sys->blkdof_pos[b + 0]     ;
        nconn = h->pimpl->sys->blkdof_pos[b + 1] - p1;

        hybsys_cellcontrib_symm(b, nconn, p1, p2, h->pimpl->gpress,
                                bsrc, h->pimpl->sys->Binv,
                                h->pimpl->hsys);

        hybsys_global_assemble_cell(nconn, h->pimpl->sys->blkdof + p1,
                                    h->pimpl->hsys->S, h->pimpl->hsys->r,
                                    h->A, h->b);

        p2 += nconn * nconn;
    }

    /* Remove zero eigenvalue */
    h->A->sa[0] *= 2;
}


/* ---------------------------------------------------------------------- */
void
ifsh_ms_press_flux(grid_t *G, struct ifsh_ms_data *h,
                   double *cpress, double *fflux)
/* ---------------------------------------------------------------------- */
{
    int    b, f, i, j, n, dof, *c;
    double s;

    MAT_SIZE_T nrows, ncols, lda, incx, incy;
    double a1, a2;

    hybsys_compute_press_flux(h->pimpl->ct->nblocks,
                              h->pimpl->sys->blkdof_pos,
                              h->pimpl->sys->blkdof,
                              h->pimpl->gpress, h->pimpl->sys->Binv,
                              h->pimpl->hsys, h->x, h->pimpl->bpress,
                              h->pimpl->hflux, h->pimpl->work);

    vector_zero(h->pimpl->ct->nfaces, h->pimpl->flux);

    for (b = i = 0; b < h->pimpl->ct->nblocks; b++) {
        for (; i < h->pimpl->sys->blkdof_pos[b + 1]; i++) {
            dof = h->pimpl->sys->blkdof[ i ];

            f = h->pimpl->sys->dof2conn[ dof ];
            s = 2.0*(h->pimpl->ct->neighbours[2*f + 0] == b) - 1.0;

            assert (f < h->pimpl->ct->nfaces);

            h->pimpl->flux[ f ] += s * h->pimpl->hflux[ i ];
        }
    }

    average_flux(h->pimpl->ct->nfaces,
                 h->pimpl->ct->neighbours, h->pimpl->flux);

    vector_zero(G->number_of_faces, fflux);

    incx = incy = 1;
    a1   = 1.0;
    a2   = 0.0;

    for (b = i = 0; b < h->pimpl->ct->nblocks; b++) {
        /* Derive coarse-scale (projected) hc fluxes for block */
        for (; i < h->pimpl->sys->blkdof_pos[b + 1]; i++) {
            dof = h->pimpl->sys->blkdof[ i ];

            f = h->pimpl->sys->dof2conn[ dof ];
            s = 2.0*(h->pimpl->ct->neighbours[2*f + 0] == b) - 1.0;

            h->pimpl->hflux[ i ] = s * h->pimpl->flux[ f ];
        }

        /* Construct fs hc fluxes in block (\Psi_b * v_c) */
        ncols  = i - h->pimpl->sys->blkdof_pos[b]; /* ndof in block */
        nrows  = h->pimpl->sys->basis_pos[b + 1] -
                 h->pimpl->sys->basis_pos[b + 0];  /* NUMEL(Psi(:)) */
        nrows /= ncols;

        lda = nrows;
        dgemv_("No Transpose", &nrows, &ncols, &a1,
               h->pimpl->sys->basis + h->pimpl->sys->basis_pos[b],
               &lda, h->pimpl->hflux + h->pimpl->sys->blkdof_pos[b],
               &incx, &a2, h->pimpl->fs_hflux, &incy);

        /* Derive cell pressure (constant per block) and accumulate fs
         * interface fluxes (internal interfaces visited twice). */
        n = 0;
        for (c  = h->pimpl->b2c + h->pimpl->pb2c[b + 0];
             c != h->pimpl->b2c + h->pimpl->pb2c[b + 1]; c++) {

            cpress[*c] = h->pimpl->bpress[b];

            for (j = G->cell_facepos[*c + 0];
                 j < G->cell_facepos[*c + 1]; j++, n++) {
                f = G->cell_faces[j];
                s = 2.0*(G->face_cells[2*f + 0] == *c) - 1.0;

                fflux[ f ] += s * h->pimpl->fs_hflux[ n ];
            }
        }
    }

    average_flux(G->number_of_faces, G->face_cells, fflux);
}
