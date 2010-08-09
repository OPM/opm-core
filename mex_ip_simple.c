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
    int nodes_ok, faces_ok, cells_ok;

    size_t field_no;
    mxArray *pm;

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

    return nodes_ok && faces_ok && cells_ok;
}


/* ------------------------------------------------------------------ */
static int
verify_structural_consistency(int nlhs, int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int rock_ok;
    int grid_ok;
    int conn_ok;

    mxAssert((nlhs == 1) && (nrhs == 4),
             "Must be called with exactly four inputs and one output.");

    mxAssert(mxIsStruct(prhs[0]) && mxIsStruct(prhs[1]),
             "First two inputs must be STRUCTs.");

    mxAssert((mxGetNumberOfElements(prhs[0]) == 1) &&
             (mxGetNumberOfElements(prhs[1]) == 1),
             "First two inputs must be 1-by-1 STRUCT arrays.");

    grid_ok = verify_grid_structure(prhs[0]);

    /* rock must contain a field 'perm' */
    rock_ok = mxGetFieldNumber(prhs[1], "perm") >= 0;

    /* connection structure */
    conn_ok = (mxIsDouble(prhs[2]) || mxIsInt32(prhs[2])) &&
              (mxIsDouble(prhs[3]) || mxIsInt32(prhs[3])) &&
              (mxGetNumberOfElements(prhs[2]) <
               mxGetNumberOfElements(prhs[3]));
    
    return grid_ok && rock_ok && conn_ok;
}


/* ------------------------------------------------------------------ */
static int
extract_face_data(const mxArray *G, int d, int **fneighbour,
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
static int
extract_cell_data(const mxArray *G,
                  const mxArray *M_nconn, const mxArray *M_conn,
                  int d, int *max_ncf, int *sum_nconn, int *sum_nconn2,
                  int **ncfaces, int **nconn, int **conn,
                  double **ccentroids, double **cvolumes)
/* ------------------------------------------------------------------ */
{
   *ncfaces = getCellFacePos(G);
   int ncells = getNumberOfCells(G);
   int i;

   extract_connection_data(M_nconn, M_conn, ncells,
                           nconn, conn, sum_nconn, sum_nconn2);
   
   *max_ncf = 0;

   for(i=0; i<ncells; ++i)
   {
      int n = (*ncfaces)[i+1] - (*ncfaces)[i];
      (*ncfaces)[i] = n;
      
      *max_ncf   = MAX(*max_ncf, n);
   }

   *ccentroids = getCellCentroids(G);
   *cvolumes   = getCellVolumes(G);
   
   return ncells;
}


/* ------------------------------------------------------------------ */
static void
deallocate_cell_data(int *ncfaces, int *nconn, int *conn, double *ccentroids)
/* ------------------------------------------------------------------ */
{
    mxFree(ccentroids);
    mxFree(conn);
    mxFree(nconn);
    mxFree(ncfaces);
}


/* ------------------------------------------------------------------ */
static void
deallocate_perm_data(double *K)
/* ------------------------------------------------------------------ */
{
    mxFree(K);
}


/*
 * BI = mex_ip_simple(G, rock, nconn, conn)
 */

/* ------------------------------------------------------------------ */
void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
/* ------------------------------------------------------------------ */
{
    int ncells, nfaces, max_ncf, sum_ncf, sum_ncf2, d, structure_ok;

    int *fneighbour;
    double *farea, *fnormal, *fcentroid;

    int *ncfaces, *nconn, *conn;
    double *ccentroids, *cvolumes;

    double *perm;
    double *Binv;

    structure_ok = verify_structural_consistency(nlhs, nrhs, prhs);

    if (structure_ok) {
        d = getNumberOfDimensions(prhs[0]);

        nfaces = extract_face_data(prhs[0], d,
                                   &fneighbour, &farea, &fnormal, &fcentroid);

        ncells = extract_cell_data(prhs[0], prhs[2], prhs[3], d,
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
}
