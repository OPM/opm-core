#include <string.h>

#include <mex.h>
#define MAT_SIZE_T mwSignedIndex

#include "call_umfpack.h"
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
static int
count_cellconn(int nc, const int *pconn)
/* ---------------------------------------------------------------------- */
{
    int c, nconn, max_nconn;

    max_nconn = 0;

    for (c = 0; c < nc; c++) {
        nconn = pconn[c + 1] - pconn[c];

        max_nconn = MAX(max_nconn, nconn);
    }

    return max_nconn;
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
allocate_aux_arrays(int nc, int nconn_tot,
                    double **src, double **gflux)
/* ---------------------------------------------------------------------- */
{
    int    ret;
    double *s, *g;

    s = mxMalloc(nc        * sizeof *s);
    g = mxMalloc(nconn_tot * sizeof *g);

    if ((s == NULL) || (g == NULL)) {
        deallocate_aux_arrays(NULL, s, g);

        *src   = NULL;
        *gflux = NULL;

        ret = 0;
    } else {
        *src   = s;
        *gflux = g;

        ret = 1;
    }

    return ret;
}


/* ---------------------------------------------------------------------- */
static int
get_pconn(const mxArray *M_pconn, int **pconn)
/* ---------------------------------------------------------------------- */
{
    int ret;

    size_t e, ne;

    int    *pi;
    double *pd;

    ne = mxGetNumberOfElements(M_pconn);

    *pconn = mxMalloc(ne * sizeof **pconn);

    if (*pconn != NULL) {
        if (mxIsDouble(M_pconn)) {
            pd = mxGetPr(M_pconn);

            for (e = 0; e < ne; e++) { (*pconn)[e] = pd[e] - 1; }
        } else {
            pi = mxGetData(M_pconn);

            for (e = 0; e < ne; e++) { (*pconn)[e] = pi[e] - 1; };
        }

        ret = ne - 1;
    } else {
        ret = -1;
    }

    return ret;
}


#if defined(ASSEMBLE_AND_SOLVE_UMFPACK) && ASSEMBLE_AND_SOLVE_UMFPACK
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
   int N, i, nf;

   N = 0;
   for (i=0; i<nc; ++i)
   {
      N += nconn[i];
   }

   nf = 0;
   for(i=0; i<N; ++i)
   {
      nf = MAX(nf, conn[i]+1);
   }

   return nf;
}
#endif


/*
 * [S, r, F, L] = mex_schur_comp_symm(BI, connPos, conns)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, nc, nconn_tot, max_nconn, sum_nconn2;
    int p2, c, i, nconn, *pconn;
    double *Binv, *ptr, *src, *gflux;
    struct hybsys *sys;

#if defined(ASSEMBLE_AND_SOLVE_UMFPACK) && ASSEMBLE_AND_SOLVE_UMFPACK
    int *conn, nf;
    double *b, *x;
    struct Sparse A;
#endif

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        nc = get_pconn(prhs[1], &pconn);

        nconn_tot = pconn[nc];
        max_nconn = count_cellconn(nc, pconn);

        allocate_aux_arrays(nc, nconn_tot, &src, &gflux);

        sum_nconn2 = mxGetNumberOfElements(prhs[0]);
        plhs[0]    = mxCreateDoubleMatrix(sum_nconn2, 1, mxREAL);
        plhs[1]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL);
        plhs[2]    = mxCreateDoubleMatrix(nconn_tot,  1, mxREAL);
        plhs[3]    = mxCreateDoubleMatrix(nc,         1, mxREAL);

        sys = hybsys_allocate(max_nconn, nc, nconn_tot);
        hybsys_init(max_nconn, nconn_tot, sys);

        for (i = 0; i < nc; i++)        { src[i]   = 0.0; } /* No sources */
        for (i = 0; i < nconn_tot; i++) { gflux[i] = 0.0; } /* No gravity */
#if 0
        src[0] = 1;
        src[nc-1]=-1;
#endif
        Binv = mxGetPr(prhs[0]);

        hybsys_compute_components(nc, pconn, gflux, src, Binv, sys);

        ptr = mxGetPr(plhs[0]);
        p2 = 0;
        for (c = 0; c < nc; c++) {
            nconn = pconn[c + 1] - pconn[c];

            hybsys_compute_cellmatrix(c, nconn, pconn[c], p2, Binv, sys);

            memcpy(ptr + p2, sys->S, nconn * nconn * sizeof *ptr);

            p2 += nconn * nconn;
        }

#if defined(ASSEMBLE_AND_SOLVE_UMFPACK) && ASSEMBLE_AND_SOLVE_UMFPACK
        conn = mxMalloc(mxGetNumberOfElements(prhs[2]) * sizeof *conn);
        get_conn(prhs[2], conn);

        nf          = get_number_of_faces(nc, nconn, conn);
        hybsys_assemble(nc, nf, nconn, conn, ptr, sys->r, &A, &b);
        x = malloc(A.n * sizeof *x);

        callMWUMFPACK(A.n, A.ia, A.ja, A.sa, b, x);

        for (i = 0; i < nf; ++i) mexPrintf("x[%d]=%f\n", i, x[i]);
        free(A.ia); free(A.ja); free(A.sa); free(x);  free(b);

        mxFree(conn);
#endif

        ptr = mxGetPr(plhs[1]);
        memcpy(ptr, sys->r, nconn_tot * sizeof *ptr);

        ptr = mxGetPr(plhs[2]);
        memcpy(ptr, sys->F, nconn_tot * sizeof *ptr);

        ptr = mxGetPr(plhs[3]);
        memcpy(ptr, sys->L, nc        * sizeof *ptr);

        hybsys_free(sys);
        deallocate_aux_arrays(pconn, src, gflux);
    }
}
