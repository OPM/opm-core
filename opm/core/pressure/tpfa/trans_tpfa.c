#include "config.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>


#ifdef __cplusplus
#include "TransTpfa.hpp"
#endif

/* ---------------------------------------------------------------------- */
/* htrans <- sum(C(:,i) .* K(cellNo,:) .* N(:,j), 2) ./ sum(C.*C, 2) */
/* ---------------------------------------------------------------------- */
void
tpfa_htrans_compute(struct UnstructuredGrid *G, const double *perm, double *htrans)
/* ---------------------------------------------------------------------- */
{
    #ifdef __cplusplus
    return tpfa_htrans_compute<UnstructuredGrid>(G, totmob, htrans, trans);
    #endif
    
    int    c, d, f, i, j;
    double s, dist, denom;

    double Kn[3];
    double *cc, *fc, *n;
    const double *K;

    MAT_SIZE_T nrows, ncols, ldA, incx, incy;
    double a1, a2;

    d = G->dimensions;

    nrows = ncols    = ldA = d;
    incx  = incy     = 1      ;
    a1    = 1.0;  a2 = 0.0    ;

    for (c = i = 0; c < G->number_of_cells; c++) {
        K  = perm + (c * d * d);
        cc = G->cell_centroids + (c * d);

        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];
            s = 2.0*(G->face_cells[2*f + 0] == c) - 1.0;

            n  = G->face_normals   + (f * d);
            fc = G->face_centroids + (f * d);

            dgemv_("No Transpose", &nrows, &ncols,
                   &a1, K, &ldA, n, &incx, &a2, &Kn[0], &incy);

            htrans[i] = denom = 0.0;
            for (j = 0; j < d; j++) {
                dist = fc[j] - cc[j];

                htrans[i] += s * dist * Kn[j];
                denom     +=     dist * dist;
            }

            assert (denom > 0);
            htrans[i] /= denom;
            htrans[i]  = fabs(htrans[i]);
        }
    }
}


/* ---------------------------------------------------------------------- */
void
tpfa_trans_compute(struct UnstructuredGrid *G, const double *htrans, double *trans)
/* ---------------------------------------------------------------------- */
{
    #ifdef __cplusplus
    return tpfa_trans_compute<UnstructuredGrid>(G, totmob, htrans, trans);
    #endif
    
    int c, i, f;

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 0.0;
    }

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            trans[f] += 1.0 / htrans[i];
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 1.0 / trans[f];
    }
}


/* ---------------------------------------------------------------------- */
void
tpfa_eff_trans_compute(struct UnstructuredGrid       *G,
                       const double *totmob,
                       const double *htrans,
                       double       *trans)
/* ---------------------------------------------------------------------- */
{
    #ifdef __cplusplus
    return tpfa_eff_trans_compute<UnstructuredGrid>(G, totmob, htrans, trans);
    #endif
    
    int c, i, f;

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 0.0;
    }

    for (c = i = 0; c < G->number_of_cells; c++) {
        for (; i < G->cell_facepos[c + 1]; i++) {
            f = G->cell_faces[i];

            trans[f] += 1.0 / (totmob[c] * htrans[i]);
        }
    }

    for (f = 0; f < G->number_of_faces; f++) {
        trans[f] = 1.0 / trans[f];
    }
}
