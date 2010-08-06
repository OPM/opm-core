#include <string.h>

#include "mex.h"
#include "matrix.h"

#include "hybsys.h"

#define MAX(a,b) (((a) > (b)) ? (a) : (b))

/* ---------------------------------------------------------------------- */
static int
verify_args(int nlhs, int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok;

    ok =       (nlhs == 4) && (nrhs == 3);
    ok = ok && mxIsDouble(prhs[0]);
    ok = ok && (mxIsDouble(prhs[1]) || mxIsInt32(prhs[1]));
    ok = ok && (mxGetNumberOfElements(prhs[0]) >
                mxGetNumberOfElements(prhs[1]));

    return ok;
}


/* ---------------------------------------------------------------------- */
static void
count_cf(const mxArray *nconn, int *nc, int *max_ncf, int *ncf_tot)
/* ---------------------------------------------------------------------- */
{
    int    c, *pi;
    double    *pd;

    *nc = mxGetNumberOfElements(nconn);

    *max_ncf = *ncf_tot = 0;

    if (mxIsDouble(nconn)) {
        pd = mxGetPr(nconn);

        for (c = 0; c < *nc; c++) {
            *max_ncf  = MAX(*max_ncf, pd[c]);
            *ncf_tot += pd[c];
        }
    } else {
        pi = mxGetData(nconn);

        for (c = 0; c < *nc; c++) {
            *max_ncf  = MAX(*max_ncf, pi[c]);
            *ncf_tot += pi[c];
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
deallocate_aux_arrays(int *nconn, double *src, double *gflux)
/* ---------------------------------------------------------------------- */
{
    /* Apparently mxFree() makes no guarantee regarding NULL arguments.
     * Institute a belt-and-suspenders approach to releasing resources.
     */
    if (gflux != NULL) { mxFree(gflux); }
    if (src   != NULL) { mxFree(src);   }
    if (nconn != NULL) { mxFree(nconn); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int nc, int ncf_tot,
                    int **nconn, double **src, double **gflux)
/* ---------------------------------------------------------------------- */
{
    int    ret, *n;
    double *s, *g;

    n = mxMalloc(nc      * sizeof *n);
    s = mxMalloc(nc      * sizeof *s);
    g = mxMalloc(ncf_tot * sizeof *g);

    if ((n == NULL) || (s == NULL) || (g == NULL)) {
        deallocate_aux_arrays(n, s, g);

        *nconn = NULL;
        *src   = NULL;
        *gflux = NULL;

        ret = 0;
    } else {
        *nconn = n;
        *src   = s;
        *gflux = g;

        ret = 1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static void
get_nconn(const mxArray *M_nconn, int *nconn)
/* ---------------------------------------------------------------------- */
{
    size_t c, nc;

    int    *pi;
    double *pd;

    nc = mxGetNumberOfElements(M_nconn);

    if (mxIsDouble(M_nconn)) {
        pd = mxGetPr(M_nconn);

        for (c = 0; c < nc; c++) { nconn[c] = pd[c]; }
    } else {
        pi = mxGetData(M_nconn);

        for (c = 0; c < nc; c++) { nconn[c] = pi[c]; };
    }
}

/* ---------------------------------------------------------------------- */
static void
get_conn(const mxArray *M_conn, int *conn)
/* ---------------------------------------------------------------------- */
{
    size_t nel, i;

    int    *pi;
    double *pd;

    nel = mxGetNumberOfElements(M_conn);

    if (mxIsDouble(M_conn)) {
        pd = mxGetPr(M_conn);

        for (i = 0; i < nel; i++) { conn[i] = pd[i] - 1; }
    } else {
        pi = mxGetData(M_conn);

        for (i = 0; i < nel; i++) { conn[i] = pi[i] - 1; }
    }
}

/* ---------------------------------------------------------------------- */
static int 
get_number_of_faces(int nc, int *nconn, int *conn)
/* ---------------------------------------------------------------------- */
{
   int N = 0;
   int i;

   for (i=0; i<nc; ++i)
   {
      N += nconn[i];
   }

   int nf=0; 
   for(i=0; i<N; ++i)
   {
      nf = MAX(nf, conn[i]+1);
   }
   
   return nf;
}


/*
 * [S, r, F, L] = mex_schur_comp_symm(BI, nconn)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, nc, ncf_tot, max_ncf, sum_ncf2, p1, p2, c, i, *nconn;
    double *Binv, *ptr, *src, *gflux;
    struct hybsys *sys;

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        count_cf(prhs[1], &nc, &max_ncf, &ncf_tot);

        allocate_aux_arrays(nc, ncf_tot, &nconn, &src, &gflux);

        sum_ncf2 = mxGetNumberOfElements(prhs[0]);
        plhs[0]  = mxCreateDoubleMatrix(sum_ncf2, 1, mxREAL);
        plhs[1]  = mxCreateDoubleMatrix(ncf_tot,  1, mxREAL);
        plhs[2]  = mxCreateDoubleMatrix(ncf_tot,  1, mxREAL);
        plhs[3]  = mxCreateDoubleMatrix(nc,       1, mxREAL);

        sys = hybsys_allocate(max_ncf, nc, ncf_tot);
        hybsys_init(max_ncf, ncf_tot, sys);

        for (i = 0; i < nc; i++)      { src[i]   = 0.0; } /* No sources */
        for (i = 0; i < ncf_tot; i++) { gflux[i] = 0.0; } /* No gravity */

        Binv = mxGetPr(prhs[0]);
        get_nconn(prhs[1], nconn);

        hybsys_compute_components(nc, nconn, gflux, src, Binv, sys);

        ptr = mxGetPr(plhs[0]);
        p1 = p2 = 0;
        for (c = 0; c < nc; c++) {
            hybsys_compute_cellmatrix(c, nconn[c], p1, p2, Binv, sys);

            memcpy(ptr + p2, sys->S, nconn[c] * nconn[c] * sizeof *ptr);

            p1 += nconn[c];
            p2 += nconn[c] * nconn[c];
        }
        
        int *conn= mxMalloc(mxGetNumberOfElements(prhs[2])*sizeof *conn);
        get_conn(prhs[2], conn);

        int  nf          = get_number_of_faces(nc, nconn, conn);           
        struct Sparse *A = hybsys_assemble(nc, nf, nconn, conn, ptr, sys->r);
        free(A->ia);  free(A->ja); free(A->sa); free(A);
        
        ptr = mxGetPr(plhs[1]);
        memcpy(ptr, sys->r, ncf_tot * sizeof *ptr);

        ptr = mxGetPr(plhs[2]);
        memcpy(ptr, sys->F, ncf_tot * sizeof *ptr);

        ptr = mxGetPr(plhs[3]);
        memcpy(ptr, sys->L, nc      * sizeof *ptr);

        hybsys_free(sys);
        deallocate_aux_arrays(nconn, src, gflux);
        mxFree(conn);
    }
}
