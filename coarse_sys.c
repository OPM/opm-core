/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathemathics
 */
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


#ifndef DEBUG_OUTPUT
#define DEBUG_OUTPUT 0
#endif

#if DEBUG_OUTPUT
#include <stdio.h>
#endif


#include "blas_lapack.h"
#include "coarse_conn.h"
#include "coarse_sys.h"
#include "hybsys.h"
#include "hybsys_global.h"
#include "mimetic.h"
#include "partition.h"
#include "sparse_sys.h"


#if defined(MAX)
#undef MAX
#endif
#define MAX(a,b)  (((a) > (b)) ? (a) : (b))


struct coarse_sys_meta {
    size_t max_ngconn;
    size_t sum_ngconn2;

    size_t max_blk_cells;
    size_t max_blk_nhf;
    size_t max_blk_nintf;
    size_t max_blk_sum_nhf2;
    size_t max_cf_nf;

    int *blk_nhf;
    int *blk_nintf;

    int *loc_fno;
    int *ncf;
    int *pconn2;

    int *pb2c, *b2c;
    
    int *data;
};


struct bf_asm_data {
    struct hybsys *fsys;

    struct CSRMatrix *A;
    double           *b;
    double           *x;

    double           *v;
    double           *p;

    double *store;
};


/* ---------------------------------------------------------------------- */
static void
coarse_sys_meta_destroy(struct coarse_sys_meta *m)
/* ---------------------------------------------------------------------- */
{
    if (m != NULL) {
        free(m->data);
    }

    free(m);
}


/* ---------------------------------------------------------------------- */
struct coarse_sys_meta *
coarse_sys_meta_allocate(size_t nblocks, size_t nfaces, size_t nc)
/* ---------------------------------------------------------------------- */
{
    size_t                  alloc_sz;
    struct coarse_sys_meta *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        alloc_sz  = nblocks;     /* blk_nhf */
        alloc_sz += nblocks;     /* blk_nintf */
        alloc_sz += nc;          /* ncf */
        alloc_sz += nc + 1;      /* pconn2 */
        alloc_sz += nfaces;      /* loc_fno */
        alloc_sz += nblocks + 1; /* pb2c */
        alloc_sz += nc;          /* b2c */

        new->data = calloc(alloc_sz, sizeof *new->data);

        if (new->data == NULL) {
            coarse_sys_meta_destroy(new);

            new = NULL;
        } else {
            new->blk_nhf   = new->data;
            new->blk_nintf = new->blk_nhf   + nblocks;
            new->ncf       = new->blk_nintf + nblocks;
            new->pconn2    = new->ncf       + nc;
            new->loc_fno   = new->pconn2    + nc + 1;

            new->pb2c      = new->loc_fno   + nfaces;
            new->b2c       = new->pb2c      + nblocks + 1;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
coarse_sys_meta_fill(int nc, const int *pgconn,
                     size_t nneigh, const int *neigh,
                     const int *p,
                     struct coarse_topology *ct,
                     struct coarse_sys_meta *m)
/* ---------------------------------------------------------------------- */
{
    int    c1, b1, c2, b2, i, n;
    size_t f, blk_sum_nhf2;

    m->max_blk_nhf = m->max_blk_nintf = 0;

    for (f = 0; f < nneigh; f++) {
        c1 = neigh[2*f + 0];   b1 = (c1 >= 0) ? p[c1] : -1;
        c2 = neigh[2*f + 1];   b2 = (c2 >= 0) ? p[c2] : -1;

        assert ((b1 >= 0) || (b2 >= 0));

        if (b1 == b2) {
            m->blk_nintf[b1] += 1;
            m->max_blk_nintf  = MAX(m->max_blk_nintf,
                                    (size_t) m->blk_nintf[b1]);
        }

        if (b1 >= 0) {
            m->blk_nhf[b1] += 1;
            m->max_blk_nhf  = MAX(m->max_blk_nhf,
                                  (size_t) m->blk_nhf[b1]);
        }

        if (b2 >= 0) {
            m->blk_nhf[b2] += 1;
            m->max_blk_nhf  = MAX(m->max_blk_nhf,
                                  (size_t) m->blk_nhf[b2]);
        }
    }

    memset(m->loc_fno, -1, nneigh * sizeof *m->loc_fno);

    m->max_cf_nf = 0;

    for (f = 0; f < (size_t) ct->nfaces; f++) {
        m->max_cf_nf = MAX(m->max_cf_nf,
                           (size_t) (ct->subfacepos[f + 1] -
                                     ct->subfacepos[f]));
    }

    m->max_ngconn = m->sum_ngconn2 = 0;
    for (c1 = 0; c1 < nc; c1++) {
        n = pgconn[c1 + 1] - pgconn[c1];

        m->max_ngconn   = MAX(m->max_ngconn, (size_t) n);
        m->sum_ngconn2 += n * n;

        m->ncf[c1]        = n;
        m->pconn2[c1 + 1] = m->pconn2[c1] + (n * n);
    }

    partition_invert(nc, p, m->pb2c, m->b2c);

    m->max_blk_cells = 0;
    for (b1 = 0; b1 < ct->nblocks; b1++) {
        m->max_blk_cells = MAX(m->max_blk_cells,
                               (size_t) (m->pb2c[b1 + 1] -
                                         m->pb2c[b1]));
    }

    m->max_blk_sum_nhf2 = 0;
    for (b1 = 0; b1 < ct->nblocks; b1++) {
        blk_sum_nhf2 = 0;

        for (i = m->pb2c[b1]; i < m->pb2c[b1 + 1]; i++) {
            c1 = m->b2c[i];
            blk_sum_nhf2 += m->pconn2[c1 + 1] - m->pconn2[c1];
        }

        m->max_blk_sum_nhf2 = MAX(m->max_blk_sum_nhf2, blk_sum_nhf2);
    }
}


/* ---------------------------------------------------------------------- */
static struct coarse_sys_meta *
coarse_sys_meta_construct(size_t nc, const int *pgconn,
                          size_t nneigh, const int *neigh,
                          const int *p, struct coarse_topology *ct)
/* ---------------------------------------------------------------------- */
{
    struct coarse_sys_meta *m;

    m = coarse_sys_meta_allocate(ct->nblocks, ct->nfaces, nc);

    if (m != NULL) {
        coarse_sys_meta_fill(nc, pgconn, nneigh, neigh, p, ct, m);
    }

    return m;
}


/* max(diff(p(1:n))) */
/* ---------------------------------------------------------------------- */
static int
max_diff(int n, int *p)
/* ---------------------------------------------------------------------- */
{
    int i, d, ret;

    assert ((n > 0) && (p != NULL));

    ret = p[1] - p[0];               assert (ret >= 0);
    for (i = 1; i < n; i++) {
        d   = p[i + 1] - p[i];
        ret = (d > ret) ? d : ret;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
void
coarse_sys_destroy(struct coarse_sys *sys)
/* ---------------------------------------------------------------------- */
{
    if (sys != NULL) {
        free(sys->Binv);
        free(sys->cell_ip);
        free(sys->basis);
        free(sys->blkdof);

        free(sys->ip_pos);
        free(sys->bf_pos);
        free(sys->dof_pos);
    }

    free(sys);
}


/* ---------------------------------------------------------------------- */
static double *
compute_fs_ip(grid_t *g, double *perm, const struct coarse_sys_meta *m)
/* ---------------------------------------------------------------------- */
{
    double *Binv;

    Binv = malloc(m->sum_ngconn2 * sizeof *Binv);

    if (Binv != NULL) {
        mim_ip_simple_all(g->number_of_cells, g->dimensions, m->max_ngconn,
                          m->ncf, g->cell_facepos, g->cell_faces,
                          g->face_cells, g->face_centroids, g->face_normals,
                          g->face_areas, g->cell_centroids, g->cell_volumes,
                          perm, Binv);
    }
    
    return Binv;
}


/* ---------------------------------------------------------------------- */
static double *
perm_weighting(size_t nc, size_t nd,
               const double *perm,
               const double *cvol)
/* ---------------------------------------------------------------------- */
{
    size_t c, d, off;
    double t;
    double *w;

    w = malloc(nc * sizeof *w);

    if (w != NULL) {
        for (c = off = 0; c < nc; c++, off += nd*nd) {
            t = 0.0;

            for (d = 0; d < nd; d++) {
                t += perm[off + d*(nd + 1)];
            }

            w[c] = t * cvol[c];
        }
    }

    return w;
}


/* ---------------------------------------------------------------------- */
static int
enforce_explicit_source(size_t nc, size_t nb, const int *p,
                        const double *src, struct coarse_sys_meta *m,
                        double *w)
/* ---------------------------------------------------------------------- */
{
    int     i, ret;
    size_t  c, b;
    int    *has_src;

    ret     = 0;
    has_src = calloc(nb, sizeof *has_src);

    if (has_src != NULL) {
        for (c = 0; c < nc; c++) {
            has_src[p[c]] += fabs(src[c]) > 0.0;
        }

        for (b = 0; b < nb; b++) {
            if (has_src[b]) {
                for (i = m->pb2c[b]; i < m->pb2c[b + 1]; i++) {
                    w[m->b2c[i]] = 0.0;
                }
            }
        }

        for (c = 0; c < nc; c++) {
            if (fabs(src[c]) > 0.0) {
                w[c] = src[c];
            }
        }

        ret = 1;
    }

    free(has_src);

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
normalize_weighting(size_t nc, size_t nb, const int *p, double *w)
/* ---------------------------------------------------------------------- */
{
    int     ret;
    size_t  c, b;
    double *bw;

    ret = 0;
    bw  = malloc(nb * sizeof *bw);

    if (bw != NULL) {
        for (b = 0; b < nb; b++) { bw[ b  ]  =  0.0; }
        for (c = 0; c < nc; c++) { bw[p[c]] += w[c]; }

        for (c = 0; c < nc; c++) {
            assert (fabs(bw[p[c]]) > 0.0);

            w[c] /= bw[p[c]];
        }

        ret = 1;
    }

    free(bw);

    return ret;
}


/* ---------------------------------------------------------------------- */
static double *
coarse_weight(grid_t *g, size_t nb,
              const int              *p,
              struct coarse_sys_meta *m,
              const double           *perm, const double *src)
/* ---------------------------------------------------------------------- */
{
    int     ok;
    double *w;

    ok = 0;
    w  = perm_weighting(nb, g->dimensions, perm, g->cell_volumes);

    if (w != NULL) {
        ok = enforce_explicit_source(g->number_of_cells,
                                     nb, p, src, m, w);
    }

    if (ok) {
        ok = normalize_weighting(g->number_of_cells, nb, p, w);
    }

    if (!ok) { free(w);  w = NULL; }

    return w;
}


/* ---------------------------------------------------------------------- */
static void
bf_asm_data_deallocate(struct bf_asm_data *data)
/* ---------------------------------------------------------------------- */
{
    if (data != NULL) {
        free            (data->store);
        csrmatrix_delete(data->A);
        hybsys_free     (data->fsys);
    }

    free(data);
}


/* ---------------------------------------------------------------------- */
static struct bf_asm_data *
bf_asm_data_allocate(struct coarse_topology *ct,
                     struct coarse_sys_meta *m)
/* ---------------------------------------------------------------------- */
{
    size_t              max_nhf, max_cells, max_faces;
    struct bf_asm_data *new;

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        max_nhf   = 2 * m->max_blk_nhf;
        max_cells = 2 * m->max_blk_cells;
        max_faces = 2 * m->max_blk_nintf + m->max_cf_nf;
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static struct coarse_sys *
coarse_sys_allocate(struct coarse_topology *ct,
                    struct coarse_sys_meta *m)
/* ---------------------------------------------------------------------- */
{
}


/* ---------------------------------------------------------------------- */
static struct coarse_sys *
coarse_sys_define(grid_t *g, const double *Binv, const double *w,
                  struct coarse_topology *ct, struct coarse_sys_meta *m)
/* ---------------------------------------------------------------------- */
{
    struct coarse_sys *sys;

    sys = coarse_sys_allocate(ct, m);

    if (sys != NULL) {
    }

    return sys;
}


/* ---------------------------------------------------------------------- */
struct coarse_sys *
coarse_sys_construct(grid_t *g, const int   *p,
                     struct coarse_topology *ct,
                     const double           *perm,
                     const double           *src,
                     const double           *totmob)
/* ---------------------------------------------------------------------- */
{
    double                 *Binv, *w;
    struct coarse_sys_meta *m;
    struct coarse_sys      *sys;

    sys = NULL;  Binv = NULL;  w = NULL;

    m = coarse_sys_meta_construct(g->number_of_cells, g->cell_facepos,
                                  g->number_of_faces, g->face_cells,
                                  p, ct);

    if (m != NULL) {
        Binv = compute_fs_ip(g, perm, m);
        w    = coarse_weight(g, ct->nblocks, p, m, perm, src);
    }

#if 0
    if ((Binv != NULL) && (w != NULL)) {
        sys = coarse_sys_define();
    }
#endif

    free(w);
    free(Binv);
    coarse_sys_meta_destroy(m);

    return sys;
}


/* Compute \Psi'_i * B * \Psi_j for all basis function pairs (i,j) for
 * all cells.  Inverts inv(B) (i.e., Binv) in each cell.  Iterates
 * over blocks (CSR representation b2c_pos, b2c).  Result store in
 * sys->cell_ip, a packed representation of the IP pairs (one col per
 * cell per block).
 *
 * Allocates work arrays and may fail.  Does currently not report failure.*/
/* ---------------------------------------------------------------------- */
void
coarse_sys_compute_cell_ip(int                nc,
                           int                max_nconn,
                           int                nb,
                           const int         *pconn,
                           const double      *Binv,
                           const int         *b2c_pos,
                           const int         *b2c,
                           struct coarse_sys *sys)
/* ---------------------------------------------------------------------- */
{
    int i, i1, i2, b, c, n, bf, *pconn2;
    int max_nbf, nbf, loc_nc;

    size_t p, nbf_pairs, bf_off, bf_sz;

    MAT_SIZE_T mm, nn, kk, nrhs, ld1, ld2, ld3, info;

    double a1, a2;
    double *work, *BI, *Psi, *IP;

#if DEBUG_OUTPUT
    FILE *fp;
#endif


    max_nbf = max_diff(nb, sys->dof_pos);

    pconn2  = malloc((nc + 1) * sizeof *pconn2);
    work    = malloc(((max_nconn * max_nconn) + /* BI */
                      (max_nconn * max_nbf  ) + /* Psi */
                      (max_nbf   * max_nbf  ))  /* IP */
                     * sizeof *work);

    if ((pconn2 != NULL) && (work != NULL)) {
        BI  = work + 0                      ;
        Psi = BI   + (max_nconn * max_nconn);
        IP  = Psi  + (max_nconn * max_nbf  );

        pconn2[0] = 0;

        for (i = 1; i <= nc; i++) {
            n         = pconn[i] - pconn[i - 1];
            pconn2[i] = pconn2[i - 1] + (n * n);
        }

#if DEBUG_OUTPUT
        fp = fopen("debug_out.m", "wt");
#endif

        for (b = 0; b < nb; b++) {
            loc_nc = b2c_pos[b + 1] - b2c_pos[b];
            bf_off = 0;
            nbf    = sys->dof_pos[b + 1] - sys->dof_pos[b];

            assert ((sys->bf_pos[b + 1] -
                     sys->bf_pos[b + 0]) % nbf == 0);

            bf_sz  = (sys->bf_pos[b + 1] - sys->bf_pos[b]) / nbf;

            nbf_pairs = nbf * (nbf + 1) / 2;

            for (i = 0; i < loc_nc; i++) {
                c = b2c[b2c_pos[b] + i];
                n = pconn[c + 1] - pconn[c];

                /* Linearise (collect) BF values for cell */
                p = sys->bf_pos[b] + bf_off;
                for (bf = 0; bf < nbf; bf++, p += bf_sz) {
                    memcpy(Psi + bf*n, &sys->basis[p], n * sizeof *Psi);
                }
#if DEBUG_OUTPUT
                fprintf(fp, "Psi{%d,%d} = [\n", b+1, i+1);
                for (i2 = 0; i2 < nbf; i2++) {
                    for (i1 = 0; i1 < n; i1++) {
                        fprintf(fp, "%22.15e ", Psi[i1 + i2*n]);
                    }
                    fprintf(fp, ";\n");
                }
                fprintf(fp, "].';\n");
#endif

                /* Extract cell's inv(B) values... */
                memcpy(BI, &Binv[pconn2[c]], n * n * sizeof *BI);
#if DEBUG_OUTPUT
                fprintf(fp, "BI{%d,%d} = [\n", b+1, i+1);
                for (i2 = 0; i2 < n; i2++) {
                    for (i1 = 0; i1 < n; i1++) {
                        fprintf(fp, "%22.15e ", BI[i1 + i2*n]);
                    }
                    fprintf(fp, ";\n");
                }
                fprintf(fp, "].';\n");
#endif

                /* ...and (Cholesky) factor it... */
                nn  = n;
                ld1 = n;
                dpotrf_("Upper Triangular", &nn, BI, &ld1, &info);

                /* ...and solve BI*X = Psi (=> Psi (=X) <- B*Psi) */
                nrhs = nbf;
                ld2  = n;
                dpotrs_("Upper Triangular", &nn, &nrhs, BI, &ld1,
                        Psi, &ld2, &info);
#if DEBUG_OUTPUT
                fprintf(fp, "BPsi{%d,%d} = [\n", b+1, i+1);
                for (i2 = 0; i2 < nbf; i2++) {
                    for (i1 = 0; i1 < n; i1++) {
                        fprintf(fp, "%22.15e ", Psi[i1 + i2*n]);
                    }
                    fprintf(fp, ";\n");
                }
                fprintf(fp, "].';\n");
#endif

                /* Finally, compute IP = Psi'*X = Psi'*B*Psi... */
                mm  = nn = nbf;
                kk  = n;
                ld1 = bf_sz;   ld2 = n;    ld3 = nbf;
                a1  = 1.0;     a2  = 0.0;
                dgemm_("Transpose", "No Transpose", &mm, &nn, &kk,
                       &a1, &sys->basis[sys->bf_pos[b] + bf_off], &ld1,
                       Psi, &ld2, &a2, IP, &ld3);
#if DEBUG_OUTPUT
                fprintf(fp, "PsiTBPsi{%d,%d} = [\n", b+1, i+1);
                for (i2 = 0; i2 < nbf; i2++) {
                    for (i1 = 0; i1 < nbf; i1++) {
                        fprintf(fp, "%22.15e ", IP[i1 + i2*nbf]);
                    }
                    fprintf(fp, ";\n");
                }
                fprintf(fp, "].';\n");
#endif

                /* ...and fill results into ip-contrib for this cell... */
                p = sys->ip_pos[b] + i*nbf_pairs;
                for (i2 = 0; i2 < nbf; i2++) {
                    for (i1 = 0; i1 <= i2; i1++, p++) {
                        sys->cell_ip[p] = IP[i1 + i2*nbf];
                    }
                }

                /* ...and prepare for next cell. */
                bf_off += n;
            }
        }
#if DEBUG_OUTPUT
        fclose(fp);
#endif
    }

    free(work);  free(pconn2);
}


/* Compute inv(B) on coarse scale from fine-scale contributions.
 * Specifically, this function computes the inverse of
 * B=sum(1/lambda_c * B_c, c\in Blk_j) for all blocks, 'j'.  The
 * fine-scale B contributions are computed in
 * coarse_sys_compute_cell_ip() above.
 *
 * Does not fail. */
/* ---------------------------------------------------------------------- */
void
coarse_sys_compute_Binv(int                nb,
                        int                max_bcells,
                        const double      *totmob,
                        const int         *b2c_pos,
                        const int         *b2c,
                        struct coarse_sys *sys,
                        double            *work)
/* ---------------------------------------------------------------------- */
{
    int    b, i, i1, i2, nbf_pairs, loc_nc, nbf;

    double a1, a2;
    double *Lti, *B;

    size_t p1, p2;

    MAT_SIZE_T mm, nn, ld, incx, incy, info;

#if DEBUG_OUTPUT
    FILE *fp;
#endif

    Lti = work + 0;
    B   = work + max_bcells;

    incx = incy = 1;
    p2   = 0;
    for (b = 0; b < nb; b++) {
        loc_nc = b2c_pos[b + 1] - b2c_pos[b];

        for (i = 0; i < loc_nc; i++) {
            Lti[i] = 1.0 / totmob[b2c[b2c_pos[b] + i]];
        }

        /* Form coarse inner-product matrix for block 'b' as (inverse)
         * mobility weighted sum of cell contributions */
        nbf       = sys->dof_pos[b + 1] - sys->dof_pos[b];
        nbf_pairs = nbf * (nbf + 1) / 2;

        mm = ld = nbf_pairs;
        nn = loc_nc;
        a1 = 1.0;
        a2 = 0.0;
        dgemv_("No Transpose", &mm, &nn, &a1,
               &sys->cell_ip[sys->ip_pos[b]], &ld,
               Lti, &incx, &a2, B, &incy);

        /* Factor (packed) SPD inner-product matrix... */
        mm = nbf;
        dpptrf_("Upper Triangular", &mm, B, &info);
        if (info == 0) {
            /* ...and invert it... */
            dpptri_("Upper Triangular", &mm, B, &info);
        } else {
#if DEBUG_OUTPUT
            fp = fopen("debug_out.m", "at");
            mm = ld = nbf_pairs;
            nn = loc_nc;
            a1 = 1.0;
            a2 = 0.0;
            dgemv_("No Transpose", &mm, &nn, &a1,
                   &sys->cell_ip[sys->ip_pos[b]], &ld,
                   Lti, &incx, &a2, B, &incy);
            fprintf(fp, "IP{%d} = [\n", b + 1);
            for (i2 = 0; i2 < nn; i2++) {
                for (i1 = 0; i1 < mm; i1++) {
                    fprintf(fp, "%18.12e ",
                            sys->cell_ip[sys->ip_pos[b] + i1 + i2*mm]);
                }
                fprintf(fp, ";\n");
            }
            fprintf(fp, "].';\n");

            fprintf(fp, "B{%d} = [ %% info = %lu\n", b + 1,
                    (unsigned long) info);
            for (i1 = 0; i1 < nbf_pairs; i1++)
                fprintf(fp, "%22.15e ", B[i1]);
            fprintf(fp, "].';\n");
            fclose(fp);
#endif
        }

        /* ...and write result to permanent storage suitable for
         * reduction functions hybsys_schur_comp*() */
        p1 = 0;
        for (i2 = 0; i2 < nbf; i2++) {
            for (i1 = 0; i1 <= i2; i1++, p1++) {
                sys->Binv[p2 + i1 + i2*nbf] = B[p1]; /* col i2 */
                sys->Binv[p2 + i2 + i1*nbf] = B[p1]; /* row i2 */
            }
        }

        p2 += nbf * nbf;
    }
}
