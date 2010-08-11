#include <assert.h>
#include <string.h>
#include <mex.h>

#include "mrst_api.h"
#include "mimetic.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))


/* ------------------------------------------------------------------ */
static int
verify_faces_structure(mxArray *faces)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(faces, "neighbors") >= 0);
    ok = ok && (mxGetFieldNumber(faces, "areas"    ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "normals"  ) >= 0);
    ok = ok && (mxGetFieldNumber(faces, "centroids") >= 0);

    return ok;
}


/* ------------------------------------------------------------------ */
static int
verify_cells_structure(mxArray *cells)
/* ------------------------------------------------------------------ */
{
    /* Shallow structural inspection only.  Assume valid fields... */
    int ok;

    ok =       (mxGetFieldNumber(cells, "facePos"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "faces"    ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "volumes"  ) >= 0);
    ok = ok && (mxGetFieldNumber(cells, "centroids") >= 0);

    return ok;
}


/* ------------------------------------------------------------------ */
static int
verify_grid_structure(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    int nodes_ok = 0, faces_ok = 0, cells_ok = 0, field_no;

    mxArray *pm;

    if (mxIsStruct(G)) {
        nodes_ok = mxGetFieldNumber(G, "nodes") >= 0;

        field_no = mxGetFieldNumber(G, "faces");
        faces_ok = field_no >= 0;
        if (faces_ok) {
            pm = mxGetFieldByNumber(G, 0, field_no);
            faces_ok = verify_faces_structure(pm);
        }

        field_no = mxGetFieldNumber(G, "cells");
        cells_ok = field_no >= 0;
        if (cells_ok) {
            pm = mxGetFieldByNumber(G, 0, field_no);
            cells_ok = verify_cells_structure(pm);
        }
    }

    return nodes_ok && faces_ok && cells_ok;
}


/* ------------------------------------------------------------------ */
static int
verify_well_structure(const mxArray *W)
/* ------------------------------------------------------------------ */
{
    int well_ok = mxIsStruct(W);

    if (well_ok && (mxGetNumberOfElements(W) > 0)) {
        /* All STRUCT array elements have the same fields... */
        well_ok =            (mxGetFieldNumber(W, "cells") >= 0);
        well_ok = well_ok && (mxGetFieldNumber(W, "WI"   ) >= 0);
    }

    return well_ok;
}


/* ------------------------------------------------------------------ */
static int
args_ok(int nlhs, int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int grid_ok = 0, rock_ok = 0, well_ok;

    if (nrhs > 0) {
        grid_ok = verify_grid_structure(prhs[0]);
    }

    if (nrhs > 1) {
        /* rock must contain a field 'perm' */
        rock_ok = mxIsStruct(prhs[1]) &&
                  (mxGetFieldNumber(prhs[1], "perm") >= 0);
    }

    /* well structure */
    if (nrhs == 2) {
        well_ok = 1;
    } else {
        well_ok = verify_well_structure(prhs[2]);
    }

    return (nlhs == 3) && grid_ok && rock_ok && well_ok;
}


/* ------------------------------------------------------------------ */
static int
extract_face_data(const mxArray *G, int **fneighbour,
                  double **farea, double **fnormal, double **fcentroid)
/* ------------------------------------------------------------------ */
{
   *fneighbour = getFaceCellNeighbors(G);
   *farea      = getFaceAreas(G);
   *fnormal    = getFaceNormals(G);
   *fcentroid  = getFaceCentroids(G);

   return getNumberOfFaces(G);
}


/* ------------------------------------------------------------------ */
static void
deallocate_face_data(int *fneighbour, double *fnormal, double *fcentroid)
/* ------------------------------------------------------------------ */
{
    mxFree(fcentroid);
    mxFree(fnormal);
    mxFree(fneighbour);
}


/* ------------------------------------------------------------------ */
static void
copy_int_vector(const mxArray *M_v, int *v)
/* ------------------------------------------------------------------ */
{
    size_t n, i;

    int    *pi;
    double *pd;

    n = mxGetNumberOfElements(M_v);

    if (mxIsDouble(M_v)) {
        pd = mxGetPr(M_v);

        for (i = 0; i < n; i++) {
            mxAssert ((INT_MIN <= pd[i]) && (pd[i] <= INT_MAX),
                      "Data element outside range for INT32");

            v[i] = pd[i];
        }
    } else {
        pi = mxGetData(M_v);
        memcpy(v, pi, n * sizeof *v);
    }
}


/* ------------------------------------------------------------------ */
static void
extract_connection_data(const mxArray *M_nconn, const mxArray *M_conn,
                        const int nc, int **nconn, int **conn,
                        int *sum_nconn, int *sum_nconn2)
/* ------------------------------------------------------------------ */
{
    size_t i, n;

    n = mxGetNumberOfElements(M_conn);

    *nconn = mxMalloc(nc * sizeof *nconn);
    *conn  = mxMalloc(n  * sizeof *conn );

    copy_int_vector(M_nconn, *nconn);
    copy_int_vector(M_conn,  *conn );

    /* Adjust connections for 1-based indexing. */
    for (i = 0; i < n; i++) {
        (*conn)[i] -= 1;
    }

    /* Count connections */
    n = nc;
    *sum_nconn = *sum_nconn2 = 0;
    for (i = 0; i < n; i++) {
        *sum_nconn  += (*nconn)[i];
        *sum_nconn2 += (*nconn)[i] * (*nconn)[i];
    }
}


/* ------------------------------------------------------------------ */
static void
extract_grid_connections(const mxArray *G, int *nc, int *nf,
                         int **ncf, int **pconn, int **conn)
/* ------------------------------------------------------------------ */
{
    int i;

    *nc = getNumberOfCells(G);
    *nf = getNumberOfFaces(G);

    /* Yes, two calls to getCellFacePos() */
    *ncf   = getCellFacePos(G);
    *pconn = getCellFacePos(G);

    *conn  = getCellFaces(G);

    for (i = 0; i < *nc; i++) {
        (*ncf)[i] = (*ncf)[i + 1] - (*ncf)[i];
    }
}


/* ------------------------------------------------------------------ */
static void
add_well_connections(const mxArray *W, int nc, int nf,
                     int *nw, int *pconn, int **conn)
/* ------------------------------------------------------------------ */
{
    int    c, i, np, fld, neconn, w, p, n, dst, src;
    int    *pi, *cwork, *c2w, *w2c, *tmp;
    double *pd;

    mxArray *M_perfs;

    if ((*nw = mxGetNumberOfElements(W)) > 0) {
        fld = mxGetFieldNumber(W, "cells");
    }

    neconn = 0;                 /* Number of additional/extra conns */
    for (w = 0; w < *nw; w++) {
        neconn += mxGetNumberOfElements(mxGetFieldByNumber(W, w, fld));
    }

    if (neconn > 0) {
        cwork = mxCalloc(nc    ,  sizeof *cwork); /* all(cwork == 0) */
        c2w   = mxMalloc(neconn * sizeof *c2w  );
        w2c   = mxMalloc(neconn * sizeof *w2c  );

        i = 0;
        for (w = 0; w < *nw; w++) {
            M_perfs = mxGetFieldByNumber(W, w, fld);
            np = mxGetNumberOfElements(M_perfs);

            if (mxIsInt32(M_perfs)) {
                pi = mxGetData(M_perfs);

                for (p = 0; p < np; p++, i++) {
                    c = pi[p] - 1;

                    w2c[i] = c;
                    c2w[i] = w;

                    cwork[c]++;
                }
            } else {
                pd = mxGetPr(M_perfs);

                for (p = 0; p < np; p++, i++) {
                    c = pd[p] - 1;

                    w2c[i] = c;
                    c2w[i] = w;

                    cwork[c]++;
                }
            }
        }

        tmp = mxRealloc(*conn, (pconn[nc] + neconn) * sizeof *tmp);
        if (tmp != NULL) {
            n          = pconn[nc] - pconn[nc - 1];  /* # exisiting conns */
            pconn[nc] += neconn;                     /* ubnd, *LAST* bin */

            for (c = nc - 1; c > 0; c--) {
                dst = pconn[c + 1] - (n + cwork[c]);
                src = pconn[c + 0];

                /* Move existing connections to start of new bin.
                 * Note: Bins may overlap so use memmove(). */
                memmove(tmp + dst, tmp + src, n * sizeof *tmp);

                /* Set cell's pointer for new connections. */
                cwork[c] = dst + n;

                n        = pconn[c] - pconn[c - 1];
                pconn[c] = dst;
            }


            /* Assign well connections.
             * Note: Generally, new DOF is nf + # ext. wells + c2w[i].
             *       In this case, though, # ext. wells == 0 ... */
            for (i = 0; i < neconn; i++) {
                assert (cwork[w2c[i]] < pconn[w2c[i] + 1]);
                tmp[cwork[w2c[i]]++] = nf + c2w[i];
            }

            /* Finally, assign updated connection structure. */
            *conn = tmp;
        }
        
        mxFree(w2c);  mxFree(c2w);  mxFree(cwork);
    }
}


/* ------------------------------------------------------------------ */
static void
count_connections(int nc, const int *pconn,
                  int *max_nconn, int *sum_nconn2)
/* ------------------------------------------------------------------ */
{
    int c, nconn;

    *max_nconn = *sum_nconn2 = 0;

    for (c = 0; c < nc; c++) {
        nconn = pconn[c + 1] - pconn[c];

        *max_nconn   = MAX(*max_nconn, nconn);
        *sum_nconn2 += nconn * nconn;
    }
}


/* ------------------------------------------------------------------ */
static int
extract_cell_data(const mxArray *G,
                  double **ccentroids, double **cvolumes)
/* ------------------------------------------------------------------ */
{
    *ccentroids = getCellCentroids(G);
    *cvolumes   = getCellVolumes(G);

    return getNumberOfCells(G);
}


/* ------------------------------------------------------------------ */
static void
deallocate_cell_data(int *ncf, int *pconn, int *conn, double *ccentroids)
/* ------------------------------------------------------------------ */
{
    mxFree(ccentroids);
    mxFree(conn);
    mxFree(pconn);
    mxFree(ncf);
}


/* ------------------------------------------------------------------ */
static void
deallocate_perm_data(double *K)
/* ------------------------------------------------------------------ */
{
    mxFree(K);
}


/* ------------------------------------------------------------------ */
static void
set_connections_M(int nc, const int *pconn, const int *conn,
                  int *pconn_M, int *conn_M)
/* ------------------------------------------------------------------ */
{
    int i;

    for (i = 0; i <= nc       ; i++) { pconn_M[i] = pconn[i] + 1; }
    for (i = 0; i <  pconn[nc]; i++) { conn_M [i] = conn [i] + 1; }
}


/*
 * [BI, connPos, conns] = mex_ip_simple(G, rock)
 * [BI, connPos, conns] = mex_ip_simple(G, rock, W)
 */

/* ------------------------------------------------------------------ */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int d, nc, nf, nw = 0, max_nconn, sum_nconn2;
    int *pconn, *conn, *ncf, *fneighbour;
    double *fcentroid, *fnormal, *farea, *ccentroids, *cvolumes, *perm;
    char errmsg[1023 + 1];

    double *Binv;

    if (args_ok(nlhs, nrhs, prhs)) {
        extract_grid_connections(prhs[0], &nc, &nf, &ncf, &pconn, &conn);

        if (nrhs == 3) {
            /* [...] = mex_ip_simple(G, rock, W) */
            add_well_connections(prhs[2], nc, nf, &nw, pconn, &conn);
        } else {
            nw = 0;
        }

        count_connections(nc, pconn, &max_nconn, &sum_nconn2);

        extract_face_data(prhs[0], &fneighbour, &farea,
                          &fnormal, &fcentroid);

        extract_cell_data(prhs[0], &ccentroids, &cvolumes);

        d    = getNumberOfDimensions(prhs[0]);
        perm = getPermeability(mxGetField(prhs[1], 0, "perm"), d);

        /* Create return values */
        plhs[0] = mxCreateDoubleMatrix (sum_nconn2, 1,                mxREAL);
        plhs[1] = mxCreateNumericMatrix(nc + 1    , 1, mxINT32_CLASS, mxREAL);
        plhs[2] = mxCreateNumericMatrix(pconn[nc] , 1, mxINT32_CLASS, mxREAL);

        /* Compute IP for all cells (reservoir) */
        Binv    = mxGetPr(plhs[0]); /* mxCreateDoubleMatrix == ZEROS */

        /* mim_ip_simple_all() requires zeroed target array... */
        mim_ip_simple_all(nc, d, max_nconn, ncf, pconn, conn,
                          fneighbour, fcentroid, fnormal, farea,
                          ccentroids, cvolumes, perm, Binv);

        set_connections_M(nc, pconn, conn,
                          mxGetData(plhs[1]),
                          mxGetData(plhs[2]));

        deallocate_perm_data(perm);
        deallocate_cell_data(ncf, pconn, conn, ccentroids);
        deallocate_face_data(fneighbour, fnormal, fcentroid);
    } else {
        sprintf(errmsg,
                "Calling sequence is\n"
                "\t[BI, connPos, conn] = %s(G, rock)\nor\n"
                "\t[BI, connPos, conn] = %s(G, rock, W)\n",
                mexFunctionName(), mexFunctionName());

        mexErrMsgTxt(errmsg);
    }
}
#if 0
    int ncells, nfaces, max_ncf, sum_ncf, sum_ncf2, d, structure_ok;

    int *fneighbour;
    double *farea, *fnormal, *fcentroid;

    int *ncfaces, *nconn, *conn;
    double *ccentroids, *cvolumes;

    double *perm;
    double *Binv;

    structure_ok = verify_structural_consistency(nlhs, nrhs, prhs);

    if (structure_ok) {
        extract_grid_connections(

        d = getNumberOfDimensions(prhs[0]);

        nfaces = extract_face_data(prhs[0], &fneighbour, &farea,
                                   &fnormal, &fcentroid);

        ncells = extract_cell_data(prhs[0], prhs[2], prhs[3],
                                   &max_ncf, &sum_ncf, &sum_ncf2, &ncfaces,
                                   &nconn, &conn, &ccentroids, &cvolumes);

        perm = getPermeability(mxGetField(prhs[1], 0, "perm"), d);

        if ((ncells > 0) && (nfaces > 0)) {
            plhs[0] = mxCreateDoubleMatrix(sum_ncf2, 1, mxREAL);
            Binv    = mxGetPr(plhs[0]); /* mxCreateDoubleMatrix == ZEROS */

            /* mim_ip_simple_all() requires zeroed target array... */
            mim_ip_simple_all(ncells, d, max_ncf, ncfaces, nconn, conn,
                              fneighbour, fcentroid, fnormal, farea,
                              ccentroids, cvolumes, perm, Binv);
        } else {
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
        }

        deallocate_perm_data(perm);
        deallocate_cell_data(ncfaces, nconn, conn, ccentroids);
        deallocate_face_data(fneighbour, fnormal, fcentroid);
    } else {
        plhs[0] = mxCreateDoubleScalar(mxGetNaN());
    }
#endif
