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
count_conns(const mxArray *M_pconn, int *nc, int *max_nconn, int *nconn_tot)
/* ---------------------------------------------------------------------- */
{
    int    c, nconn, *pi;
    double           *pd;

    *nc = mxGetNumberOfElements(M_pconn) - 1;

    *max_nconn = *nconn_tot = 0;

    if (mxIsDouble(M_pconn)) {
        pd = mxGetPr(M_pconn);

        for (c = 0; c < *nc; c++) {
            nconn = pd[c + 1] - pd[c];
            *max_nconn  = MAX(*max_nconn, nconn);
            *nconn_tot += nconn;
        }
    } else {
        pi = mxGetData(M_pconn);

        for (c = 0; c < *nc; c++) {
            nconn = pi[c + 1] - pi[c];
            *max_nconn  = MAX(*max_nconn, nconn);
            *nconn_tot += nconn;
        }
    }
}


/* ---------------------------------------------------------------------- */
static void
deallocate_aux_arrays(int *pconn, int *conn,
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
    if (pconn != NULL) { mxFree(pconn); }
}


/* ---------------------------------------------------------------------- */
static int
allocate_aux_arrays(int max_nconn, int nc, int nconn_tot,
                    int **pconn, int **conn,
                    double **src, double **gflux,
                    double **work)
/* ---------------------------------------------------------------------- */
{
    int    ret, *p, *c;
    double *s, *g, *w;

    p = mxMalloc((nc + 1)  * sizeof *p);
    c = mxMalloc(nconn_tot * sizeof *c);
    s = mxMalloc(nc        * sizeof *s);
    g = mxMalloc(nconn_tot * sizeof *g);
    w = mxMalloc(max_nconn * sizeof *w);

    if ((p == NULL) || (c == NULL) ||
        (s == NULL) || (g == NULL) || (w == NULL)) {
        deallocate_aux_arrays(p, c, s, g, w);

        *pconn = NULL;
        *conn  = NULL;
        *src   = NULL;
        *gflux = NULL;
        *work  = NULL;

        ret = 0;
    } else {
        *pconn = p;
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
copy_M_int_vector(const mxArray *M_a, int *a)
/* ---------------------------------------------------------------------- */
{
    size_t nel, i;

    int    *pi;
    double *pd;

    nel = mxGetNumberOfElements(M_a);

    if (mxIsDouble(M_a)) {
        pd = mxGetPr(M_a);

        for (i = 0; i < nel; i++) { a[i] = pd[i] - 1; }
    } else {
        pi = mxGetData(M_a);

        for (i = 0; i < nel; i++) { a[i] = pi[i] - 1; }
    }
}


/* ---------------------------------------------------------------------- */
static void
get_pconn(const mxArray *M_pconn, int *pconn)
/* ---------------------------------------------------------------------- */
{
    copy_M_int_vector(M_pconn, pconn);
}


/* ---------------------------------------------------------------------- */
static void
get_conn(const mxArray *M_conn, int *conn)
/* ---------------------------------------------------------------------- */
{
    copy_M_int_vector(M_conn, conn);
}


/*
 * [v, p] = mex_compute_press_flux(BI, pi, connPos, conns, F, L)
 */

/* ---------------------------------------------------------------------- */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ---------------------------------------------------------------------- */
{
    int ok, i, nc, max_nconn, nconn_tot, *pconn, *conn;
    double *src, *gflux, *work, *Binv, *pi, *ptr;
    struct hybsys *sys;

    ok = verify_args(nlhs, nrhs, prhs);

    if (ok) {
        count_conns(prhs[2], &nc, &max_nconn, &nconn_tot);

        allocate_aux_arrays(max_nconn, nc, nconn_tot,
                            &pconn, &conn, &src, &gflux, &work);

        plhs[0] = mxCreateDoubleMatrix(nconn_tot, 1, mxREAL);
        plhs[1] = mxCreateDoubleMatrix(nc,        1, mxREAL);

        sys = hybsys_allocate(max_nconn, nc, nconn_tot);
        hybsys_init(max_nconn, nconn_tot, sys);

        ptr = mxGetPr(prhs[4]);
        memcpy(sys->F, ptr, nconn_tot * sizeof *sys->F);

        ptr = mxGetPr(prhs[5]);
        memcpy(sys->L, ptr, nc        * sizeof *sys->L);

        get_pconn(prhs[2], pconn);
        get_conn (prhs[3], conn);

        Binv = mxGetPr(prhs[0]);
        pi   = mxGetPr(prhs[1]);

        for (i = 0; i < nc; i++)        { src[i]   = 0.0; } /* No sources */
        for (i = 0; i < nconn_tot; i++) { gflux[i] = 0.0; } /* No gravity */

        hybsys_compute_press_flux(nc, pconn, conn, gflux, src, Binv, sys,
                                  pi, mxGetPr(plhs[1]), mxGetPr(plhs[0]),
                                  work, max_nconn);

        hybsys_free(sys);
        deallocate_aux_arrays(pconn, conn, src, gflux, work);
    }
}
