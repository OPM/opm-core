#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "trans_tpfa.h"

/* ---------------------------------------------------------------------- */
/* htrans <- sum(C(:,i) .* K(cellNo,:) .* N(:,j), 2) ./ sum(C.*C, 2) */
/* ---------------------------------------------------------------------- */
void
tpfa_htrans_compute(grid_t *G, const double *perm, double *htrans)
/* ---------------------------------------------------------------------- */
{
    int    c, d, f, i, j;
    double s, dist, denom;

    double Kn[3];
    double *cc, *fc, *n;
    const double *K;

    MAT_SIZE_T nrows, ncols, ldA, incx, incy;
    double a1, a2;

    d = G->dimensions;

    nrows = ncols    = ldA = d;
    incx  = incy     = 1      ;
    a1    = 1.0;  a2 = 0.0    ;

    for (c = i = 0; c < G->number_of_cells; c++) {
        K  = perm + (c * d * d);
        cc = G->cell_centroids + (c * d);

        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];
            s = 2.0*(G->face_cells[2*f + 0] == c) - 1.0;

            n  = G->face_normals   + (f * d);
            fc = G->face_centroids + (f * d);

            dgemv_("No Transpose", &nrows, &ncols,
                   &a1, K, &ldA, n, &incx, &a2, &Kn[0], &incy);

            htrans[i] = denom = 0.0;
            for (j = 0; j < d; j++) {
                dist = fc[j] - cc[j];

                htrans[i] += s * dist * Kn[j];
                denom     +=     dist * dist;
            }

            assert (denom > 0);
            htrans[i] /= denom;
            htrans[i]  = fabs(htrans[i]);
        }
    }
}


/* ---------------------------------------------------------------------- */
void
tpfa_trans_compute(grid_t *G, const double *htrans, double *trans)
/* ---------------------------------------------------------------------- */
{
    int c, i, f;

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 0.0;
    }

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            trans[f] += 1.0 / htrans[i];
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 1.0 / trans[f];
    }
}


/* ---------------------------------------------------------------------- */
void
tpfa_eff_trans_compute(grid_t       *G,
                       const double *totmob,
                       const double *htrans,
                       double       *trans)
/* ---------------------------------------------------------------------- */
{
    int c, i, f;

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 0.0;
    }

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            trans[f] += 1.0 / (totmob[c] * htrans[i]);
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 1.0 / trans[f];
    }
}


/* ---------------------------------------------------------------------- */
void
small_matvec(size_t n, int sz, const double *A, const double *X, double *Y)
/* ---------------------------------------------------------------------- */
{
    size_t i, p1, p2;

    MAT_SIZE_T nrows, ncols, ld, incx, incy;
    double     a1, a2;

    nrows = ncols = ld = sz;
    incx  = incy  = 1;

    a1 = 1.0;
    a2 = 0.0;
    for (i = p1 = p2 = 0; i < n; i++) {
        dgemv_("No Transpose", &nrows, &ncols,
               &a1, A + p2, &ld, X + p1, &incx,
               &a2,              Y + p1, &incy);

        p1 += sz;
        p2 += sz * sz;
    }
}


/* ---------------------------------------------------------------------- */
static void
compr_htran_mult_core(grid_t       *G,
                      int           np,
                      const double *Ac,
                      const double *xf,
                      double       *ht_mult,
                      double       *luAc,
                      double       *v,
                      MAT_SIZE_T   *ipiv)
/* ---------------------------------------------------------------------- */
{
    int c, i, f, p, np2;

    size_t p2;

    MAT_SIZE_T nrows, ncols, ldA, ldX, nrhs, info;

    np2   = np * np;
    nrows = ncols = ldA = ldX = np;
    info  = 0;

    for (c = 0, p2; c < G->number_of_cells; c++, p2 += np2) {
        /* Factor Ac */
        memcpy(luAc, Ac + p2, np2 * sizeof *luAc);
        dgetrf_(&nrows, &ncols, luAc, &ldA, ipiv, &info);

        /* Define right-hand sides for local tran-mult systems */
        for (i = G->cell_facepos[c + 0], nrhs = 0;
             i < G->cell_facepos[c + 1]; i++, nrhs++) {
            f = G->cell_faces[i];

            for (p = 0; p < np; p++) {
                v[nrhs*np + p] = xf[f*np + p];
            }
        }

        /* Solve local tran-mult systems */
        dgetrs_("No Transpose", &nrows, &nrhs,
                luAc, &ldA, ipiv, v, &ldX, &info);

        /* Compute local tran-multipliers by summing over phases */
        for (i = G->cell_facepos[c + 0], nrhs = 0;
             i < G->cell_facepos[c + 1]; i++, nrhs++) {
            ht_mult[i] = 0.0;
            for (p = 0; p < np; p++) {
                ht_mult[i] += v[nrhs*np + p];
            }
        }
    }
}


/* ---------------------------------------------------------------------- */
/* ht_mult <- Ac \ xf
 *
 * np is number of phases.
 * Ac is R/B per cell.
 * xf is (R/B)*{\lambda_\alpha}_f, pre-computed in small_matvec().
 *
 * Result, ht_mult, is a scalar per half-face.
 *
 * Algorithm:
 *   for each cell,
 *      Ac <- LU(Ac)
 *      for each face(cell),
 *         Solve Ac*v = xf(f)
 *         ht_mult[c,f] = sum(v) over phases */
/* ---------------------------------------------------------------------- */
void
tpfa_compr_htran_mult(grid_t       *G ,
                      int           np,
                      size_t        max_ngconn,
                      const double *Ac,
                      const double *xf,
                      double       *ht_mult)
/* ---------------------------------------------------------------------- */
{
    double     *luAc, *v;
    MAT_SIZE_T *ipiv;

    luAc = malloc(np * np         * sizeof *luAc);
    v    = malloc(np * max_ngconn * sizeof *v);
    ipiv = malloc(np              * sizeof *ipiv);

    if ((luAc != NULL) && (v != NULL) && (ipiv != NULL)) {
        compr_htran_mult_core(G, np, Ac, xf, ht_mult,
                              luAc, v, ipiv);
    }

    free(ipiv);
    free(v);
    free(luAc);
}
