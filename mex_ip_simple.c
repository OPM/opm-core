#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "mex.h"
#include "matrix.h"


#define MAT_SIZE_T mwSignedIndex


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

    mxAssert((nlhs == 1) && (nrhs == 2),
             "Must be called with exactly two inputs and one output.");

    mxAssert(mxIsStruct(prhs[0]) && mxIsStruct(prhs[1]),
             "Inputs must be STRUCTs.");

    mxAssert((mxGetNumberOfElements(prhs[0]) == 1) &&
             (mxGetNumberOfElements(prhs[1]) == 1),
             "Inputs must be 1-by-1 STRUCT arrays.");

    grid_ok = verify_grid_structure(prhs[0]);

    /* rock must contain a field 'perm' */
    rock_ok = mxGetFieldNumber(prhs[1], "perm") >= 0;

    return grid_ok && rock_ok;
}


/* ------------------------------------------------------------------ */
static void
deallocate_face_data(int *fneighbour, double *fnormal, double *fcentroid)
/* ------------------------------------------------------------------ */
{
    free(fcentroid);
    free(fnormal);
    free(fneighbour);
}


/* ------------------------------------------------------------------ */
static int
allocate_face_data(int nfaces, int d, int **fneighbour,
                   double **fnormal, double **fcentroid)
/* ------------------------------------------------------------------ */
{
    int *neigh, ret;
    double *n, *c;

    neigh = malloc(2 * nfaces * sizeof *neigh);
    n     = malloc(d * nfaces * sizeof *n);
    c     = malloc(d * nfaces * sizeof *c);

    if ((neigh != NULL) && (n != NULL) && (c != NULL)) {
        *fneighbour = neigh;
        *fnormal    = n;
        *fcentroid  = c;

        ret = 1;
    } else {
        deallocate_face_data(neigh, n, c);

        ret = 0;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
static void
extract_double_neighbours(const mxArray *pn, int *neighbour)
/* ------------------------------------------------------------------ */
{
    int f, nf;
    double *n;

    nf = mxGetM (pn);
    n  = mxGetPr(pn);

    for (f = 0; f < nf; f++) {
        mxAssert ((n[f + 0*nf] <= INT_MAX) &&
                  (n[f + 1*nf] <= INT_MAX),
                  "Neighbour cell exceeds INT_MAX");

        neighbour[2*f + 0] = n[f + 0*nf] - 1;
        neighbour[2*f + 1] = n[f + 1*nf] - 1;
    }
}


/* ------------------------------------------------------------------ */
static void
extract_int_neighbours(const mxArray *pn, int *neighbour)
/* ------------------------------------------------------------------ */
{
    int f, nf;
    int *n;

    nf = mxGetM   (pn);
    n  = mxGetData(pn);

    for (f = 0; f < nf; f++) {
        mxAssert ((n[f + 0*nf] <= INT_MAX) &&
                  (n[f + 1*nf] <= INT_MAX),
                  "Neighbour cell exceeds INT_MAX");

        neighbour[2*f + 0] = n[f + 0*nf] - 1;
        neighbour[2*f + 1] = n[f + 1*nf] - 1;
    }
}


/* ------------------------------------------------------------------ */
static int
extract_face_data(const mxArray *faces, int d, int **fneighbour,
                  double **farea, double **fnormal, double **fcentroid)
/* ------------------------------------------------------------------ */
{
    int nfaces, alloc_ok, ret;

    mxArray *pm;

    int f, j;

    /* int    *neigh; */
    /* double *fneigh; */

    double *fn, *n;
    double *fc, *c;

    pm     = mxGetField(faces, 0, "neighbors");
    nfaces = mxGetM(pm);

    alloc_ok = allocate_face_data(nfaces, d, fneighbour,
                                  fnormal, fcentroid);
    if (alloc_ok) {
        if (mxIsDouble(pm)) {
            extract_double_neighbours(pm, *fneighbour);
        } else {
            mxAssert (mxIsInt32(pm),
                      "Neighbours is neither DOUBLE nor INT32");
            extract_int_neighbours(pm, *fneighbour);
        }

        pm = mxGetField(faces, 0, "areas");
        *farea = mxGetPr(pm);

        pm = mxGetField(faces, 0, "normals");
        fn = mxGetPr(pm);

        pm = mxGetField(faces, 0, "centroids");
        fc = mxGetPr(pm);

        /* neigh = *fneighbour; */
        n     = *fnormal;
        c     = *fcentroid;

        for (f = 0; f < nfaces; f++) {
            /* neigh[2*f + 0] = fneigh[f + 0*nfaces] - 1; */
            /* neigh[2*f + 1] = fneigh[f + 1*nfaces] - 1; */

            for (j = 0; j < d; j++) {
                n[j + f*d] = fn[f + j*nfaces];
                c[j + f*d] = fc[f + j*nfaces];
            }
        }

        ret = nfaces;
    } else {
        ret = -1;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
static void
deallocate_cell_data(int *ncfaces, int *cfaces, double *ccentroids)
/* ------------------------------------------------------------------ */
{
    free(ccentroids);
    free(cfaces);
    free(ncfaces);
}


/* ------------------------------------------------------------------ */
static int
allocate_cell_data(int ncells, int d, int ncellfaces,
                   int **ncfaces, int **cfaces, double **ccentroids)
/* ------------------------------------------------------------------ */
{
    int ret, *ncf, *cf;
    double *cc;

    ncf = malloc(ncells       * sizeof *ncf);
    cf  = malloc(ncellfaces   * sizeof *cf);
    cc  = malloc((ncells * d) * sizeof *cc);

    if ((ncf != NULL) && (cf != NULL) && (cc != NULL)) {
        *ncfaces    = ncf;
        *cfaces     = cf;
        *ccentroids = cc;

        ret = 1;
    } else {
        deallocate_cell_data(ncf, cf, cc);

        ret = 0;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
static void
extract_double_ncf(const mxArray *pfp, int *ncf,
                   int *max_ncf, int *sum_ncf, int *sum_ncf2)
/* ------------------------------------------------------------------ */
{
    int c, nc;
    double *fp;

    nc = mxGetNumberOfElements(pfp) - 1;
    fp = mxGetPr(pfp);

    *max_ncf = -1;  *sum_ncf = *sum_ncf2 = 0;
    for (c = 0; c < nc; c++) {
        ncf[c] = fp[c + 1] - fp[c];

        *max_ncf = MAX(ncf[c], *max_ncf);

        *sum_ncf  += ncf[c];
        *sum_ncf2 += ncf[c] * ncf[c];
    }
}



/* ------------------------------------------------------------------ */
static void
extract_int_ncf(const mxArray *pfp, int *ncf,
                int *max_ncf, int *sum_ncf, int *sum_ncf2)
/* ------------------------------------------------------------------ */
{
    int c, nc;
    int *fp;

    nc = mxGetNumberOfElements(pfp) - 1;
    fp = mxGetData(pfp);

    *max_ncf = -1;  *sum_ncf = *sum_ncf2 = 0;
    for (c = 0; c < nc; c++) {
        ncf[c] = fp[c + 1] - fp[c];

        *max_ncf = MAX(ncf[c], *max_ncf);

        *sum_ncf  += ncf[c];
        *sum_ncf2 += ncf[c] * ncf[c];
    }
}


/* ------------------------------------------------------------------ */
static void
extract_double_cf(const mxArray *pcf, int *cf)
/* ------------------------------------------------------------------ */
{
    int p, m;
    double *f;

    m = mxGetNumberOfElements(pcf);
    f = mxGetPr(pcf);

    for (p = 0; p < m; p++) {
        mxAssert (f[p] <= INT_MAX,
                  "cells.faces element exceeds INT_MAX");
        cf[p] = f[p] - 1;
    }
}


/* ------------------------------------------------------------------ */
static void
extract_int_cf(const mxArray *pcf, int *cf)
/* ------------------------------------------------------------------ */
{
    int p, m;
    int *f;

    m = mxGetNumberOfElements(pcf);
    f = mxGetData(pcf);

    for (p = 0; p < m; p++) {
        mxAssert (f[p] <= INT_MAX,
                  "cells.faces element exceeds INT_MAX");
        cf[p] = f[p] - 1;
    }
}


/* ------------------------------------------------------------------ */
static int
extract_cell_data(const mxArray *cells, int d,
                  int *max_ncf, int *sum_ncf, int *sum_ncf2,
                  int **ncfaces, int **cfaces,
                  double **ccentroids, double **cvolumes)
/* ------------------------------------------------------------------ */
{
    int ncells, ncellfaces, alloc_ok, ret;

    int c, i, j, pos;

    mxArray *pfp, *pcf, *pm;

    double *centroids;

    /* int *ncf, *cf; */
    double *cc;

    pfp = mxGetField(cells, 0, "facePos");
    pcf = mxGetField(cells, 0, "faces"  );

    ncells     = mxGetNumberOfElements(pfp) - 1;
    ncellfaces = mxGetNumberOfElements(pcf);

    alloc_ok = allocate_cell_data(ncells, d, ncellfaces,
                                  ncfaces, cfaces, ccentroids);
    if (alloc_ok) {
        if (mxIsDouble(pfp)) {
            extract_double_ncf(pfp, *ncfaces, max_ncf, sum_ncf, sum_ncf2);
        } else {
            mxAssert (mxIsInt32(pfp),
                      "facePos is neither DOUBLE nor INT32");
            extract_int_ncf(pfp, *ncfaces, max_ncf, sum_ncf, sum_ncf2);
        }

        if (mxIsDouble(pcf)) {
            extract_double_cf(pcf, *cfaces);
        } else {
            mxAssert (mxIsInt32(pcf),
                      "cells.faces is neither DOUBLE nor INT32");
            extract_int_cf(pcf, *cfaces);
        }
        /* facePos = mxGetPr(pfp); */
        /* faces   = mxGetPr(pcf); */

        pm = mxGetField(cells, 0, "volumes");
        *cvolumes = mxGetPr(pm);

        pm = mxGetField(cells, 0, "centroids");
        centroids = mxGetPr(pm);

        /* ncf = *ncfaces; */
        /* cf  = *cfaces; */
        cc  = *ccentroids;

        /* *max_ncf = -1; *sum_ncf = *sum_ncf2 = 0; pos = 0; */
        for (c = 0; c < ncells; c++) {
            /* ncf[c] = facePos[c + 1] - facePos[c]; */

            /* *max_ncf   = MAX(ncf[c], *max_ncf); */
            /* *sum_ncf  += ncf[c]; */
            /* *sum_ncf2 += ncf[c] * ncf[c]; */

            for (j = 0; j < d; j++) {
                cc[j + c*d] = centroids[c + j*ncells];
            }

            /* for (i = 0; i < ncf[c]; i++, pos++) { */
            /*     cf[pos] = faces[pos] - 1; */
            /* } */
        }

        ret = ncells;
    } else {
        ret = -1;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
static void
deallocate_perm_data(double *K)
/* ------------------------------------------------------------------ */
{
    free(K);
}


/* ------------------------------------------------------------------ */
static int
allocate_perm_data(int ncells, int d, double **K)
/* ------------------------------------------------------------------ */
{
    int ret;
    double *perm;

    perm = malloc(ncells * d * d * sizeof *perm);

    if (perm != NULL) {
        *K = perm;
        ret = 1;
    } else {
        ret = 0;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
static void
extract_perm_data(const mxArray *perm, int d, double **K)
/* ------------------------------------------------------------------ */
{
    int ncells, ncomp, alloc_ok;

    int c, i, off;

    double *k, *tensor;

    ncells = mxGetM(perm);
    ncomp  = mxGetN(perm);

    alloc_ok = allocate_perm_data(ncells, d, K);

    if (alloc_ok) {
        k = *K;
        for (i = 0; i < ncells * d * d; i++) {
            k[i] = 0.0;
        }

        tensor = mxGetPr(perm);

        if (ncomp == 1) {
            /* Isotropic (scalar) tensor */
            for (c = 0; c < ncells; c++) {
                off = c * d * d;
                for (i = 0; i < d; i++) {
                    k[i*(d + 1) + off] = tensor[c];
                }
            }
        } else if (ncomp == d) {
            /* Diagonal tensor */
            for (c = 0; c < ncells; c++) {
                off = c * d * d;
                for (i = 0; i < d; i++) {
                    k[i*(d + 1) + off] = tensor[c + i*ncells];
                }
            }
        } else if (d == 2) {
            /* Full 2D tensor */
            assert (ncomp == 3);

            for (c = 0; c < ncells; c++) {
                off = c * d * d;
                k[0 + off] = tensor[c + 0*ncells];
                k[1 + off] = tensor[c + 1*ncells];
                k[2 + off] = tensor[c + 1*ncells];
                k[3 + off] = tensor[c + 2*ncells];
            }
        } else {
            /* Full 3D tensor */
            assert ((d == 3) && (ncomp == 6));

            for (c = 0; c < ncells; c++) {
                off = c * d * d;

                k[0 + off] = tensor[c + 0*ncells];
                k[1 + off] = tensor[c + 1*ncells];
                k[2 + off] = tensor[c + 2*ncells];

                k[3 + off] = tensor[c + 1*ncells];
                k[4 + off] = tensor[c + 3*ncells];
                k[5 + off] = tensor[c + 4*ncells];

                k[6 + off] = tensor[c + 2*ncells];
                k[7 + off] = tensor[c + 4*ncells];
                k[8 + off] = tensor[c + 5*ncells];
            }
        }
    }
}


/* ------------------------------------------------------------------ */
static int
get_dimension(const mxArray *G)
/* ------------------------------------------------------------------ */
{
    mxArray *p1, *p2;

    p1 = mxGetField(G , 0, "nodes" );
    p2 = mxGetField(p1, 0, "coords");

    return mxGetN(p2);
}


/*
 * BI = mex_ip_simple(G, rock)
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

    int *ncfaces, *cfaces;
    double *ccentroids, *cvolumes;

    double *perm;
    double *Binv;

    structure_ok = verify_structural_consistency(nlhs, nrhs, prhs);

    if (structure_ok) {
        d = get_dimension(prhs[0]);

        nfaces = extract_face_data(mxGetField(prhs[0], 0, "faces"), d,
                                   &fneighbour, &farea, &fnormal, &fcentroid);

        ncells = extract_cell_data(mxGetField(prhs[0], 0, "cells"), d,
                                   &max_ncf, &sum_ncf, &sum_ncf2, &ncfaces,
                                   &cfaces, &ccentroids, &cvolumes);

        extract_perm_data(mxGetField(prhs[1], 0, "perm"), d, &perm);

        if ((ncells > 0) && (nfaces > 0)) {
            plhs[0] = mxCreateDoubleMatrix(sum_ncf2, 1, mxREAL);
            Binv    = mxGetPr(plhs[0]);

            mim_ip_simple_all(ncells, d, max_ncf, ncfaces, cfaces, fneighbour,
                              fcentroid, fnormal, farea, ccentroids, cvolumes,
                              perm, Binv);
        } else {
            plhs[0] = mxCreateDoubleScalar(mxGetNaN());
        }

        deallocate_perm_data(perm);
        deallocate_cell_data(ncfaces, cfaces, ccentroids);
        deallocate_face_data(fneighbour, fnormal, fcentroid);
    } else {
        plhs[0] = mxCreateDoubleScalar(mxGetNaN());
    }
}
