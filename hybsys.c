#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


#if defined(COMPILING_FOR_MATLAB) && COMPILING_FOR_MATLAB
#include <mex.h>
#define MAT_SIZE_T mwSignedIndex
#endif

#ifndef MAT_SIZE_T
#define MAT_SIZE_T int
#endif


#include "blas_lapack.h"
#include "hybsys.h"

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
        new->F1  = malloc(nconn_tot             * sizeof *new->F1 );

        if ((new->one == NULL) || (new->r  == NULL) || (new->S == NULL) ||
            (new->L   == NULL) || (new->F1 == NULL)) {
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
void
hybsys_free(struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    if (sys != NULL) {
        if (sys->F2 != sys->F1) { free(sys->F2); } /* unsymmetric system */

        free(sys->F1 );
        free(sys->L  );
        free(sys->S  );
        free(sys->r  );
        free(sys->one);
    }

    free(sys);
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
    int    c, i, p1, p2, nconn;
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
    int    c, i, p1, p2, nconn;
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
    int    c, i, p1, p2, nconn;
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
static void
hybsys_cellmat_symm_core(int nconn, const double *Binv, double L,
                         const double *F, double *S)
/* ---------------------------------------------------------------------- */
{
    int        i, j;
    MAT_SIZE_T n, k, ldA, ldC, incx, incy;
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
static void
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

    hybsys_cellrhs_core(nconn, &gpress[p1], src[c], &Binv[p2],
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

    hybsys_cellrhs_core(nconn, &gpress[p1], src[c], &Binv[p2],
                        sys->L[c], &sys->F1[p1], &sys->F2[p1],
                        sys->r);
}


/* ---------------------------------------------------------------------- */
void
hybsys_compute_press_flux(int nc, const int *pconn, const int *conn,
                          const double *gpress, const double *src,
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
            work[i] = pi[conn[p1 + i]] - gpress[p1 + i];
        }

        nrows = ncols = lda = nconn;

        /* Solve Lp = g - F2*f + F2*pi (for cell pressure) */
        press[c]  = src[c];
        press[c] += ddot_(&nrows, &sys->F2[p1], &incx, work, &incy);
        press[c] /= sys->L[c];

        /* Form rhs of system B*v = f + C*p - D*pi */
        for (i = 0; i < nconn; i++) {
            work[i] = press[c] - work[i];
        }

        /* Solve resulting system (-> half face fluxes) */
        dgemv_("No Transpose", &nrows, &ncols,
               &a1, &Binv[p2], &lda, work, &incx,
               &a2, &flux[p1],             &incy);

        p2 += nconn * nconn;
    }
}


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
