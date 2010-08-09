#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>


#if defined(COMPILING_FOR_MATLAB) && COMPILING_FOR_MATLAB
#include "mex.h"
#include "matrix.h"
#define MAT_SIZE_T mwSignedIndex
#endif

#ifndef MAT_SIZE_T
#define MAT_SIZE_T int
#endif


#include "blas_lapack.h"
#include "hybsys.h"

/* ---------------------------------------------------------------------- */
struct hybsys *
hybsys_allocate(int max_ncf, int nc, int ncf_tot)
/* ---------------------------------------------------------------------- */
{
    struct hybsys *new;

    new = malloc(1 * sizeof *new);
    if (new != NULL) {
        new->one = malloc(max_ncf           * sizeof *new->one);
        new->S   = malloc(max_ncf * max_ncf * sizeof *new->S  );
        new->L   = malloc(nc                * sizeof *new->L  );
        new->F   = malloc(ncf_tot           * sizeof *new->F  );
        new->r   = malloc(ncf_tot           * sizeof *new->r  );

        if ((new->one == NULL) || (new->S == NULL) ||
            (new->L   == NULL) || (new->F == NULL) || (new->r == NULL)) {
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
        free(sys->r  );
        free(sys->F  );
        free(sys->L  );
        free(sys->S  );
        free(sys->one);
    }

    free(sys);
}


/* ---------------------------------------------------------------------- */
void
hybsys_init(int max_ncf, int ncf_tot, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    int i;

    for (i = 0; i < max_ncf; i++) {
        sys->one[i] = 1.0;
    }

    for (i = 0; i < ncf_tot; i++) {
        sys->r[i] = 0.0;
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_compute_components(int nc, const int *nconn,
                          const double *gflux, const double *src,
                          const double *Binv, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    int    c, i, p1, p2;
    double csrc, a1, a2;

    MAT_SIZE_T incx, incy;
    MAT_SIZE_T nrows, ncols, lda;

    incx = incy = 1;
    p1 = p2 = 0;

    for (c = 0; c < nc; c++) {
        nrows = ncols = lda = nconn[c];

        /* F <- C' * inv(B) == (inv(B) * ones(n,1))' in single cell */
        a1 = 1.0;  a2 = 0.0;
        dgemv_("No Transpose", &nrows, &ncols,
               &a1, &Binv[p2]  , &lda, sys->one, &incx,
               &a2, &sys->F[p1],                 &incy);

        /* L <- C' * inv(B) * C == SUM(F) == ones(n,1)' * F */
        sys->L[c] = ddot_(&nrows, sys->one, &incx, &sys->F[p1], &incy);

        /* csrc <- g[c] - C'*v_g */
        csrc = src[c] - ddot_(&nrows, sys->one, &incx, &gflux[p1], &incy);

        /* r <- v_g */
        for (i = 0; i < nconn[c]; i++) {
            sys->r[p1 + i] = gflux[p1 + i];
        }

        /* r <- r + F'*inv(L)*(g - C'*v_g) == r + (csrc / L[c])*F */
        a1 = csrc / sys->L[c];
        daxpy_(&nrows, &a1, &sys->F[p1], &incx, &sys->r[p1], &incy);

        p1 += nconn[c];
        p2 += nconn[c] * nconn[c];
    }
}


/* ---------------------------------------------------------------------- */
void
hybsys_compute_cellmatrix_core(int nconn, const double *Binv,
                               double L, const double *F, double *S)
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
void
hybsys_compute_cellmatrix(int c, int nconn, int p1, int p2,
                          const double *Binv, struct hybsys *sys)
/* ---------------------------------------------------------------------- */
{
    hybsys_compute_cellmatrix_core(nconn, &Binv[p2], sys->L[c],
                                   &sys->F[p1], sys->S);
}


/* ---------------------------------------------------------------------- */
void
hybsys_compute_press_flux(int nc, const int *nconn, const int *conn,
                          const double *gflux, const double *src,
                          const double *Binv, const struct hybsys *sys,
                          const double *pi, double *press, double *flux,
                          double *work, const int lwork)
/* ---------------------------------------------------------------------- */
{
    int c, i, p1, p2;
    double a1, a2;

    MAT_SIZE_T incx, incy, nrows, ncols, lda;

    incx = incy = 1;

    p1 = p2 = 0;
    a1 = a2 = 1.0;
    for (c = 0; c < nc; c++) {
        assert (lwork >= nconn[c]);

        /* Serialise interface pressures for cell */
        for (i = 0; i < nconn[c]; i++) {
            work[i] = pi[conn[p1 + i]];
        }

        nrows = ncols = lda = nconn[c];

        /* Solve Lp = g - C'*v_g + F'*pi (for cell pressure) */
        press[c]  = src[c];
        press[c] -= ddot_(&nrows, sys->one   , &incx, &gflux[p1], &incy);
        press[c] += ddot_(&nrows, &sys->F[p1], &incx, work      , &incy);
        press[c] /= sys->L[c];

        /* Form rhs of system B*v = B*v_g + C*p - D*pi */
        for (i = 0; i < nconn[c]; i++) {
            flux[p1 + i] = gflux[p1 + i];
            work[i]      = press[c] - work[i];
        }

        /* Solve resulting system (-> half face fluxes) */
        dgemv_("No Transpose", &nrows, &ncols,
               &a1, &Binv[p2], &lda, work, &incx,
               &a2, &flux[p1],             &incy);

        p1 += nconn[c];
        p2 += nconn[c] * nconn[c];
    }
}


/*
 * Routines to assemble global matrix
 *
 */

/* ---------------------------------------------------------------------- */
static MAT_SIZE_T *
hybsys_build_ia(int nc, int nf, int *nconn, int *conn)
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
      int n = nconn[c];
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
hybsys_build_ja(int nc, int nf, int *nconn, int *conn, 
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
      int n = nconn[c];
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
hybsys_build_sa_and_b(int nc, int nf, int *nconn, int *conn, MAT_SIZE_T *ia, 
                      double *S, double *R, int *work, double **sa, double **b)
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
   
   double *s = S;
   double *r = R;

   int c, pos = 0; 
   for(c=0; c<nc; ++c)
   {
      int n = nconn[c];
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
void hybsys_build_matrix_structure(int nc, int nf, int *nconn, int *conn, 
                                   MAT_SIZE_T **ia, MAT_SIZE_T **ja)
/* ---------------------------------------------------------------------- */
{
   int *work = malloc(nf * sizeof *work);

   *ia       = hybsys_build_ia(nc, nf, nconn, conn);
   *ja       = hybsys_build_ja(nc, nf, nconn, conn, *ia, work);

   free(work);
}
/* ---------------------------------------------------------------------- */
void hybsys_assemble_global_system(int nc, int nf, int *nconn, int *conn, 
                                   double *S, double *R, 
                                   double **sa, double **b, MAT_SIZE_T *ia)
/* ---------------------------------------------------------------------- */
{
   int *work  = malloc(nf * sizeof *work);

   hybsys_build_sa_and_b(nc, nf, nconn, conn, ia, S, R, work, sa, b); 

   free(work);
}

/* ---------------------------------------------------------------------- */
void
hybsys_assemble(int nc, int nf, int *nconn, int *conn, double *S, double *R, struct Sparse*A, double **b)
/* ---------------------------------------------------------------------- */
{
   A->m    =  A->n   = nf;
   hybsys_build_matrix_structure(nc, nf, nconn, conn, &A->ia, &A->ja);
   hybsys_assemble_global_system(nc, nf, nconn, conn, S, R, &A->sa, b, A->ia);
}
