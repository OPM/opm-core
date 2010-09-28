/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "blas_lapack.h"
#include "hybsys.h"


#if defined(MAX)
#undef MAX
#endif

#define MAX(a,b) (((a) > (b)) ? (a) : (b))


/* ---------------------------------------------------------------------- */
struct hybsys *
hybsys_allocate_symm(int max_nconn, int nc, int nconn_tot)
/* ---------------------------------------------------------------------- */
{
    struct hybsys *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->one = malloc(max_nconn             * sizeof *new->one);
        new->r   = malloc(max_nconn             * sizeof *new->r  );
        new->S   = malloc(max_nconn * max_nconn * sizeof *new->S  );
        new->L   = malloc(nc                    * sizeof *new->L  );
        new->q   = malloc(nc                    * sizeof *new->q  );
        new->F1  = malloc(nconn_tot             * sizeof *new->F1 );

        if ((new->one == NULL) || (new->r == NULL) || (new->S  == NULL) ||
            (new->L   == NULL) || (new->q == NULL) || (new->F1 == NULL)) {
            hybsys_free(new);

            new = NULL;
        } else {
            new->F2 = new->F1;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
struct hybsys *
hybsys_allocate_unsymm(int max_nconn, int nc, int nconn_tot)
/* ---------------------------------------------------------------------- */
{
    struct hybsys *new;

    new = hybsys_allocate_symm(max_nconn, nc, nconn_tot);

    if (new != NULL) {
        new->F2 = malloc(nconn_tot * sizeof *new->F2);

        if (new->F2 == NULL) {
            hybsys_free(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
static void
hybsys_well_count_conn(int nc, const int *cwpos,
                       int *max_nw, size_t *sum_nwc)
/* ---------------------------------------------------------------------- */
{
    int c, nw;

    *max_nw  = 0;
    *sum_nwc = 0;

    for (c = 0; c < nc; c++) {
        nw = cwpos[c + 1] - cwpos[c];

        assert (nw >= 0);

        *max_nw   = MAX(*max_nw, nw);
        *sum_nwc += nw;
    }
}


/* ---------------------------------------------------------------------- */
struct hybsys_well *
hybsys_well_allocate_symm(int max_nconn, int nc, int *cwpos)
/* ---------------------------------------------------------------------- */
{
    int                 max_nw;
    size_t              sum_nwc, alloc_sz;

    struct hybsys_well *new;

    assert (cwpos[nc] > cwpos[0]); /* Else no wells. */

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        hybsys_well_count_conn(nc, cwpos, &max_nw, &sum_nwc);

        alloc_sz  = sum_nwc;            /* F1 */
        alloc_sz += max_nconn + max_nw; /* r */
        alloc_sz += max_nw * max_nconn; /* w2r */
        alloc_sz += max_nw * max_nw;    /* w2w */

        new->data = malloc(alloc_sz * sizeof *new->data);

        if (new->data != NULL) {
            new->F1  = new->data;
            new->F2  = new->F1;

            new->r   = new->F2  + sum_nwc;

            new->w2r = new->r   + max_nconn + max_nw;
            new->r2w = new->w2r;

            new->w2w = new->r2w + (max_nw * max_nconn);
        } else {
            hybsys_well_free(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
struct hybsys_well *
hybsys_well_allocate_unsymm(int max_nconn, int nc, int *cwpos)
/* ---------------------------------------------------------------------- */
{
    int                 max_nw;
    size_t              sum_nwc, alloc_sz;

    struct hybsys_well *new;

    assert (cwpos[nc] > cwpos[0]); /* Else no wells. */

    new = malloc(1 * sizeof *new);

    if (new != NULL) {
        hybsys_well_count_conn(nc, cwpos, &max_nw, &sum_nwc);

        alloc_sz  = 2 * sum_nwc;            /* F1, F2 */
        alloc_sz +=     max_nconn + max_nw; /* r */
        alloc_sz += 2 * max_nw * max_nconn; /* w2r, r2w */
        alloc_sz +=     max_nw * max_nw;    /* w2w */

        new->data = malloc(alloc_sz * sizeof *new->data);

        if (new->data != NULL) {
            new->F1  = new->data;
            new->F2  = new->F1  + sum_nwc;

            new->r   = new->F2  + sum_nwc;

            new->w2r = new->r   + max_nconn + max_nw;
            new->r2w = new->w2r + (max_nw * max_nconn);

            new->w2w = new->r2w + (max_nw * max_nconn);
        } else {
            hybsys_well_free(new);
            new = NULL;
        }
    }

    return new;
}


/* ---------------------------------------------------------------------- */
void
hybsys_free(struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    if (sys != NULL) {
        if (sys->F2 != sys->F1) { free(sys->F2); } /* unsymmetric system */

        free(sys->F1 );
        free(sys->q  );
        free(sys->L  );
        free(sys->S  );
        free(sys->r  );
        free(sys->one);
    }

    free(sys);
}


/* ---------------------------------------------------------------------- */
void
hybsys_well_free(struct hybsys_well *wsys)
/* ---------------------------------------------------------------------- */
{
    if (wsys != NULL) {
        free(wsys->data);
    }

    free(wsys);
}


/* ---------------------------------------------------------------------- */
void
hybsys_init(int max_nconn, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    int i;

    for (i = 0; i < max_nconn; i++) {
        sys->one[i] = 1.0;
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_schur_comp_symm(int nc, const int *pconn,
                       const double *Binv, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    int    c, p1, p2, nconn;
    double a1, a2;

    MAT_SIZE_T incx, incy;
    MAT_SIZE_T nrows, ncols, lda;

    incx = incy = 1;
    p1 = p2 = 0;

    for (c = 0; c < nc; c++) {
        p1    = pconn[c + 0];
        nconn = pconn[c + 1] - pconn[c];
        nrows = ncols = lda = nconn;

        /* F <- C' * inv(B) == (inv(B) * ones(n,1))' in single cell */
        a1 = 1.0;  a2 = 0.0;
        dgemv_("No Transpose"   , &nrows, &ncols,
               &a1, &Binv[p2]   , &lda, sys->one, &incx,
               &a2, &sys->F1[p1],                 &incy);

        /* L <- C' * inv(B) * C == SUM(F) == ones(n,1)' * F */
        sys->L[c] = ddot_(&nrows, sys->one, &incx, &sys->F1[p1], &incy);

        p2 += nconn * nconn;
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_schur_comp_unsymm(int nc, const int *pconn,
                         const double *Binv, const double *BIV,
                         const double *P, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    int    c, p1, p2, nconn;
    double a1, a2;

    MAT_SIZE_T incx, incy;
    MAT_SIZE_T nrows, ncols, lda;

    assert ((sys->F2 != sys->F1) &&
            (sys->F2 != NULL));

    incx = incy = 1;
    p2 = 0;

    for (c = 0; c < nc; c++) {
        p1    = pconn[c + 0];
        nconn = pconn[c + 1] - pconn[c];

        nrows = ncols = lda = nconn;

        /* F1 <- C' * inv(B) */
        a1 = 1.0;  a2 = 0.0;
        dgemv_("No Transpose"   , &nrows, &ncols,
               &a1, &Binv[p2]   , &lda, sys->one, &incx,
               &a2, &sys->F1[p1],                 &incy);

        /* F2 <- (C - V)' * inv(B) == F1 - V'*inv(B) */
        a1 = -1.0;
        memcpy(&sys->F2[p1], &sys->F1[p1], nconn * sizeof sys->F2[p1]);
        daxpy_(&nrows, &a1, &BIV[p1], &incx, &sys->F2[p1], &incy);

        /* L <- (C - V)' * inv(B) * C - P */
        sys->L[c]  = ddot_(&nrows, sys->one, &incx, &sys->F1[p1], &incy);
        sys->L[c] -= ddot_(&nrows, sys->one, &incx, &BIV[p1]    , &incy);
        sys->L[c] -= P[c];

        p2 += nconn * nconn;
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_schur_comp_gen(int nc, const int *pconn,
                      const double *Binv, const double *C2,
                      const double *P, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    int    c, p1, p2, nconn;
    double a1, a2;

    MAT_SIZE_T incx, incy;
    MAT_SIZE_T nrows, ncols, lda;

    assert ((sys->F2 != sys->F1) &&
            (sys->F2 != NULL));

    incx = incy = 1;
    p2 = 0;

    for (c = 0; c < nc; c++) {
        p1    = pconn[c + 0];
        nconn = pconn[c + 1] - pconn[c];

        nrows = ncols = lda = nconn;

        /* F1 <- C' * inv(B) */
        a1 = 1.0;  a2 = 0.0;
        dgemv_("No Transpose"   , &nrows, &ncols,
               &a1, &Binv[p2]   , &lda, sys->one, &incx,
               &a2, &sys->F1[p1],                 &incy);

        /* F2 <- C2' * inv(B) */
        dgemv_("No Transpose"   , &nrows, &ncols,
               &a1, &Binv[p2]   , &lda, &C2[p1], &incx,
               &a2, &sys->F2[p1],                &incy);

        /* L <- C2' * inv(B) * C - P == F2'*ones(n,1) - P */
        sys->L[c]  = ddot_(&nrows, sys->one, &incx, &sys->F2[p1], &incy);
        sys->L[c] -= P[c];

        p2 += nconn * nconn;
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_well_schur_comp_symm(int nc, const int *cwpos,
                            double             *WI,
                            struct hybsys      *sys,
                            struct hybsys_well *wsys)
/* ---------------------------------------------------------------------- */
{
    int c, i;

    for (c = i = 0; c < nc; c++) {
        for (; i < cwpos[c + 1]; i++) {
            wsys->F1[i]  = WI[i];
            sys->L  [c] += WI[i];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
hybsys_cellmat_symm_core(int nconn, const double *Binv, double L,
                         const double *F, double *S)
/* ---------------------------------------------------------------------- */
{
    int        i, j;
    MAT_SIZE_T n, k, ldA, ldC;
    double     a1, a2;

    /* S <- D' * inv(B) * D == inv(B) in single cell */
    memcpy(S, Binv, nconn * nconn * sizeof *S);

    /* S <- S - F'*inv(L)*F */
    n  = ldA = ldC = nconn;
    k  = 1;
    a1 = -1.0 / L;
    a2 = 1.0;
    dsyrk_("Upper Triangular", "No Transpose", &n, &k,
           &a1, F, &ldA, &a2, S, &ldC);

    /* Account for DSYRK only updating the upper triangular part of S */
    for (j = 0; j < nconn; j++) {
        for (i = j + 1; i < nconn; i++) {
            S[i + j*nconn] = S[j + i*nconn];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
hybsys_cellmat_unsymm_core(int nconn, const double *Binv, double L,
                           const double *F1, const double *F2,
                           double *S)
/* ---------------------------------------------------------------------- */
{
    MAT_SIZE_T m, n, k, ldF1, ldF2, ldS;
    double     a1, a2;

    /* S <- D' * inv(B) * D == inv(B) in single cell */
    memcpy(S, Binv, nconn * nconn * sizeof *S);

    /* S <- S - F1'*inv(L)*F2 */
    a1 = -1.0 / L;
    a2 =  1.0;

    m = n = nconn;
    k = 1;
    ldF1 = ldF2 = 1;
    ldS  = nconn;

    dgemm_("Transpose", "No Transpose", &m, &n, &k,
           &a1, F1, &ldF1, F2, &ldF2, &a2, S, &ldS);
}


/* ---------------------------------------------------------------------- */
static double
hybsys_cellrhs_core(int nconn, const double *gpress, double src,
                    const double *Binv, double L, const double *F1,
                    const double *F2, double *R)
/* ---------------------------------------------------------------------- */
{
    MAT_SIZE_T n, k, ldA, incx, incy;
    double     a1, a2;

    /* r <- inv(B)*gpress + F1'*inv(L)*(src - F2*gpress)
     *   == inv(B)*gpress + F1'*inv(L)*(src - C2'*inv(B)*gpress) */
    k  = 1;
    a1 = 1.0;  a2 = 0.0;
    incx = incy = 1;

    n = k = ldA = nconn;

    dgemv_("No Transpose", &n, &k,
           &a1, Binv, &ldA, gpress, &incx,
           &a2, R   ,               &incy);

    src -= ddot_(&n, F2, &incx, gpress, &incy);

    a1 = src / L;
    daxpy_(&n, &a1, F1, &incx, R, &incy);

    return src;
}


/* ---------------------------------------------------------------------- */
void
hybsys_cellcontrib_symm(int c, int nconn, int p1, int p2,
                        const double *gpress, const double *src,
                        const double *Binv, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    hybsys_cellmat_symm_core(nconn, &Binv[p2],
                             sys->L[c], &sys->F1[p1],
                             sys->S);

    sys->q[c] = hybsys_cellrhs_core(nconn, &gpress[p1], src[c], &Binv[p2],
                                    sys->L[c], &sys->F1[p1], &sys->F1[p1],
                                    sys->r);
}


/* ---------------------------------------------------------------------- */
void
hybsys_cellcontrib_unsymm(int c, int nconn, int p1, int p2,
                          const double *gpress, const double *src,
                          const double *Binv, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    assert ((sys->F2 != sys->F1) &&
            (sys->F2 != NULL));

    hybsys_cellmat_unsymm_core(nconn, &Binv[p2],
                               sys->L[c], &sys->F1[p1], &sys->F2[p1],
                               sys->S);

    sys->q[c] = hybsys_cellrhs_core(nconn, &gpress[p1], src[c], &Binv[p2],
                                    sys->L[c], &sys->F1[p1], &sys->F2[p1],
                                    sys->r);
}


/* ---------------------------------------------------------------------- */
void
hybsys_well_cellcontrib_symm(int c, int ngconn, int p1,
                             const int *cwpos,
                             const double *WI, const double *wdp,
                             struct hybsys *sys, struct hybsys_well *wsys)
/* ---------------------------------------------------------------------- */
{
    int        i, w, nw, wp1;
    MAT_SIZE_T mm, nn, kk, ld1, ld2, ld3, incx, incy;
    double     a1, a2, q;

    nw  = cwpos[c + 1] - cwpos[c];
    wp1 = cwpos[c];

    /* -------------------------------------------------------------- */
    /* w2r = - F1(r)'*F2(w)/L, r2w = w2r' */
    mm = ngconn;   ld1 = 1;
    nn = nw;       ld2 = 1;
    kk = 1;        ld3 = ngconn;

    a1 = -1.0 / sys->L[c];
    a2 = 0.0;

    dgemm_("Transpose", "No Transpose", &mm, &nn, &kk,
           &a1, &sys->F1[p1], &ld1, &wsys->F2[wp1], &ld2,
           &a2, wsys->w2r, &ld3);

    /* -------------------------------------------------------------- */
    /* w2w = BI - F1(w)'*F2(w)/L */
    mm = nw;   ld1 = 1;
    nn = nw;   ld2 = 1;
    kk = 1;    ld3 = nw;

    a1 = -1.0 / sys->L[c];
    a2 = 0.0;

    dgemm_("Transpose", "No Transpose", &mm, &nn, &kk,
           &a1, &wsys->F1[wp1], &ld1, &wsys->F2[wp1], &ld2,
           &a2, wsys->w2w, &ld3);

    for (w = 0; w < nw; w++) {
        wsys->w2w[w * (nw + 1)] += WI[wp1 + w];
    }

    /* -------------------------------------------------------------- */
    /* Global RHS contributions */
    mm   = nw;
    incx = incy = 1;
    q    = ddot_(&mm, &wsys->F2[wp1], &incx, &wdp[wp1], &incy);

    a1 = -q / sys->L[c];
    for (i = 0; i < ngconn; i++) {
        wsys->r[i] = a1 * sys->F1[p1 + i];
    }

    sys->q[c] -= q;
    a1 = sys->q[c] / sys->L[c];
    for (w = 0; w < nw; w++) {
        wsys->r[ngconn + w] = a1*wsys->F1[wp1 + w] +
                              WI[wp1 + w] * wdp[wp1 + w];
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_compute_press_flux(int nc, const int *pconn, const int *conn,
                          const double *gpress,
                          const double *Binv, const struct hybsys *sys,
                          const double *pi, double *press, double *flux,
                          double *work)
/* ---------------------------------------------------------------------- */
{
    int    c, i, nconn, p1, p2;
    double a1, a2;

    MAT_SIZE_T incx, incy, nrows, ncols, lda;

    incx = incy = 1;

    p2 = 0;
    a1 = 1.0;
    a2 = 0.0;
    for (c = 0; c < nc; c++) {
        p1    = pconn[c + 0];
        nconn = pconn[c + 1] - p1;

        /* Serialise interface pressures for cell */
        for (i = 0; i < nconn; i++) {
            /* work[i] = pi[conn[p1 + i]] - gpress[p1 + i]; */
            work[i] = pi[conn[p1 + i]];
        }

        nrows = ncols = lda = nconn;

        /* Solve Lp = g - F2*f + F2*pi (for cell pressure) */
        press[c]  = sys->q[c]; /* src[c]; */
        press[c] += ddot_(&nrows, &sys->F2[p1], &incx, work, &incy);
        press[c] /= sys->L[c];

        /* Form rhs of system B*v = f + C*p - D*pi */
        for (i = 0; i < nconn; i++) {
            work[i] = gpress[p1 + i] + press[c] - work[i];
        }

        /* Solve resulting system (-> half face fluxes) */
        dgemv_("No Transpose", &nrows, &ncols,
               &a1, &Binv[p2], &lda, work, &incx,
               &a2, &flux[p1],             &incy);

        p2 += nconn * nconn;
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_compute_press_flux_well(int nc, const int *pgconn, int nf,
                               int nw, const int *pwconn, const int *wconn,
                               const double *Binv,
                               const double *WI,
                               const double *wdp,
                               const struct hybsys      *sys,
                               const struct hybsys_well *wsys,
                               const double             *pi,
                               double *cpress, double *cflux,
                               double *wpress, double *wflux,
                               double *work)
/* ---------------------------------------------------------------------- */
{
    int    c, w, wg, perf;
    int    ngconn, nwconn;
    size_t gp1, gp2, wp1;

    MAT_SIZE_T mm, nn, incx, incy, ld;

    double dcp, one;

    gp2 = 0;
    for (c = 0; c < nc; c++) {
        ngconn = pgconn[c + 1] - pgconn[c];
        nwconn = pwconn[c + 1] - pwconn[c];

        if (nwconn > 0) {
            dcp = 0.0;
            
            gp1 = pgconn[c];
            wp1 = pwconn[c];

            for (w = 0; w < nwconn; w++) {
                wg      = wconn[2*(wp1 + w) + 0];
                work[w] = pi[nf + wg];
            }

            mm   = nwconn;  incx = incy = 1;
            dcp  = ddot_(&mm, &wsys->F2[wp1], &incx, work, &incy);
            dcp /= sys->L[c];

            cpress[c] += dcp;

            mm  = nn = ld = ngconn;
            one = 1.0;
            dgemv_("No Transpose", &mm, &nn,
                   &dcp, &Binv[gp2], &ld, sys->one  , &incx,
                   &one,                 &cflux[gp1], &incy);

            for (w = 0; w < nwconn; w++) {
                perf         = wconn[2*(wp1 + w) + 1];

                wflux[perf]  = wdp[wp1 + w] + cpress[c] - work[w];
                wflux[perf] *= - WI [wp1 + w]; /* Sign => positive inj. */
            }
        }

        gp2 += ngconn + ngconn;
    }

    /* Assign well BHP from linsolve output */
    memcpy(wpress, pi + nf, nw * sizeof *wpress);
}


#if 0
/*
 * Routines to assemble global matrix
 *
 */

/* ---------------------------------------------------------------------- */
static MAT_SIZE_T *
hybsys_build_ia(int nc, int nf,
                const int *pconn, const int *conn)
/* ---------------------------------------------------------------------- */
{
   MAT_SIZE_T    *ia = malloc((nf+1) * sizeof *ia);

   int i;
   for(i=0; i<nf+1; ++i)
   {
      ia[i] = 0;
   }

   /*
    *   Compute rowsizes
    */
   int c, pos = 0;
   for(c=0; c<nc; ++c)
   {
      int n = pconn[c + 1] - pconn[c];
      for (i=pos; i<pos+n; ++i)
      {
         mxAssert(conn[i]<nf, "conn out of bounds");
         ia[1+conn[i]] += n - 1;
      }
      pos += n;
   }

   /*
    *   cumulative sum...
    */
   for(i=1; i<nf+1; ++i)
   {
      ia[i] = ia[i-1] + ia[i]+1;
   }
   return ia;
}

/* ---------------------------------------------------------------------- */
static MAT_SIZE_T*
hybsys_build_ja(int nc, int nf,
                const int *pconn, const int *conn,
                MAT_SIZE_T *ia, int *work)
/* ---------------------------------------------------------------------- */
{
   MAT_SIZE_T *ja = malloc(ia[nf] * sizeof *ja);
   int  i,j;

   /*
    *   For each row, diagonal entries are positioned first.
    */
   for(i=0; i<nf; ++i)
   {
      work[i] = 1;
      ja[ia[i]] = i;
   }

   int c, pos = 0;
   for(c=0; c<nc; ++c)
   {
      int n = pconn[c + 1] - pconn[c];
      for (i=pos; i<pos+n; ++i)
      {
         int fi = conn[i];
         /* mxAssert(fi<nf, "fi out of bounds"); */

         for (j=pos; j<i; ++j)
         {
            int fj = conn[j];
            /* mxAssert(fj<nf, "fj out of bounds"); */

            /*
             *   No conditionals since off-diagonals entries are
             *   visited only once.
             */
            ja[ia[fi] + work[fi]++] = fj;
            ja[ia[fj] + work[fj]++] = fi;
         }
      }
      pos += n;
   }
   return ja;
}
/* ---------------------------------------------------------------------- */
static void
hybsys_build_sa_and_b(int nc, int nf,
                      const int *pconn, const int *conn, MAT_SIZE_T *ia,
                      const double *S, const double *R,
                      int *work, double **sa, double **b)
/* ---------------------------------------------------------------------- */
{
   *sa = malloc(ia[nf] * sizeof **sa);
   *b  = malloc(nf * sizeof **b);
   int     i,j;

   /*
    *   Clear diagonal and work array
    */
   for(i=0; i<nf; ++i)
   {
      work[i]      = 1;
      (*sa)[ia[i]] = 0;
      (*b) [i]     = 0;
   }

   const double *s = S;
   const double *r = R;

   int c, pos = 0;
   for(c=0; c<nc; ++c)
   {
      int n = pconn[c + 1] - pconn[c];
      for (i=pos; i<pos+n; ++i)
      {
         int fi = conn[i];
         int ii = i-pos;

         for (j=pos; j<i; ++j)
         {
            int fj = conn[j];
            int jj = j-pos;

            /*
             *   We can use assignment since off-diagonal entries are
             *   visited only once.
             */
            (*sa)[ia[fi] + work[fi]++] = s[ii + jj*n];
            (*sa)[ia[fj] + work[fj]++] = s[jj + ii*n];
         }

         /*
          *   Diagonal entries are visited more than once.
          */
         (*sa)[ia[fi]] += s[ii + ii*n];
         (*b) [fi]     += r[ii];
      }

      s   += n*n;
      r   += n;
      pos += n;
   }
}


/* ---------------------------------------------------------------------- */
static void
hybsys_build_matrix_structure(int nc, int nf,
                              const int *pconn, const int *conn,
                              MAT_SIZE_T **ia, MAT_SIZE_T **ja)
/* ---------------------------------------------------------------------- */
{
   int *work = malloc(nf * sizeof *work);

   *ia       = hybsys_build_ia(nc, nf, pconn, conn);
   *ja       = hybsys_build_ja(nc, nf, pconn, conn, *ia, work);

   free(work);
}


/* ---------------------------------------------------------------------- */
static void
hybsys_assemble_global_system(int nc, int nf,
                              const int *pconn, const int *conn,
                              const double *S, const double *R,
                              double **sa, double **b, MAT_SIZE_T *ia)
/* ---------------------------------------------------------------------- */
{
   int *work  = malloc(nf * sizeof *work);

   hybsys_build_sa_and_b(nc, nf, pconn, conn, ia, S, R, work, sa, b);

   free(work);
}


/* ---------------------------------------------------------------------- */
void
hybsys_assemble(int nc, int nf,
                const int *pconn, const int *conn,
                const double *S, const double *R,
                struct Sparse*A, double **b)
/* ---------------------------------------------------------------------- */
{
   A->m    =  A->n   = nf;
   hybsys_build_matrix_structure(nc, nf, pconn, conn, &A->ia, &A->ja);
   hybsys_assemble_global_system(nc, nf, pconn, conn, S, R, &A->sa, b, A->ia);
}
#endif
