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

    ok =       (nlhs == 2) && (nrhs == 6);
    ok = ok && mxIsDouble(prhs[0]);
    ok = ok && mxIsDouble(prhs[1]);
    ok = ok && (mxIsDouble(prhs[2]) || mxIsInt32(prhs[2]));
    ok = ok && (mxIsDouble(prhs[3]) || mxIsInt32(prhs[3]));
    ok = ok && mxIsDouble(prhs[4]);
    ok = ok && mxIsDouble(prhs[5]);

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
deallocate_aux_arrays(int *nconn, int *conn,
                      double *src, double *gflux, double *work)
/* ---------------------------------------------------------------------- */
{
    /* Apparently mxFree() makes no guarantee regarding NULL arguments.
     * Institute a belt-and-suspenders approach to releasing resources.
     */
    if (work  != NULL) { mxFree(work);  }
    if (gflux != NULL) { mxFree(gflux); }
    if (src   != NULL) { mxFree(src);   }
    if (conn  != NULL) { mxFree(conn);  }
    if (nconn != NULL) { mxFree(nconn); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int max_ncf, int nc, int ncf_tot,
                    int **nconn, int **conn,
                    double **src, double **gflux,
                    double **work)
/* ---------------------------------------------------------------------- */
{
    int    ret, *n, *c;
    double *s, *g, *w;

    n = mxMalloc(nc      * sizeof *n);
    c = mxMalloc(ncf_tot * sizeof *c);
    s = mxMalloc(nc      * sizeof *s);
    g = mxMalloc(ncf_tot * sizeof *g);
    w = mxMalloc(max_ncf * sizeof *w);

    if ((n == NULL) || (c == NULL) ||
        (s == NULL) || (g == NULL) || (w == NULL)) {
        deallocate_aux_arrays(n, c, s, g, w);

        *nconn = NULL;
        *conn  = NULL;
        *src   = NULL;
        *gflux = NULL;
        *work  = NULL;

        ret = 0;
    } else {
        *nconn = n;
        *conn  = c;
        *src   = s;
        *gflux = g;
        *work  = w;

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


/*
 * [v, p] = mex_compute_press_flux(BI, pi, nconn, conn, F, L)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, i, nc, max_ncf, ncf_tot, *nconn, *conn;
    double *src, *gflux, *work, *Binv, *pi, *ptr;
    struct hybsys *sys;

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        count_cf(prhs[2], &nc, &max_ncf, &ncf_tot);

        allocate_aux_arrays(max_ncf, nc, ncf_tot,
                            &nconn, &conn, &src, &gflux, &work);

        plhs[0] = mxCreateDoubleMatrix(ncf_tot, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(nc,      1, mxREAL);

        sys = hybsys_allocate(max_ncf, nc, ncf_tot);
        hybsys_init(max_ncf, ncf_tot, sys);

        ptr = mxGetPr(prhs[4]);
        memcpy(sys->F, ptr, ncf_tot * sizeof *sys->F);

        ptr = mxGetPr(prhs[5]);
        memcpy(sys->L, ptr, nc      * sizeof *sys->L);

        get_nconn(prhs[2], nconn);
        get_conn (prhs[3], conn);

        Binv = mxGetPr(prhs[0]);
        pi   = mxGetPr(prhs[1]);

        for (i = 0; i < nc; i++)      { src[i]   = 0.0; } /* No sources */
        for (i = 0; i < ncf_tot; i++) { gflux[i] = 0.0; } /* No gravity */

        hybsys_compute_press_flux(nc, nconn, conn, gflux, src, Binv, sys,
                                  pi, mxGetPr(plhs[1]), mxGetPr(plhs[0]),
                                  work, max_ncf);

        hybsys_free(sys);
        deallocate_aux_arrays(nconn, conn, src, gflux, work);
    }
}
