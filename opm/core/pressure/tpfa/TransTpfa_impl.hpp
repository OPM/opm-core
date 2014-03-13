#include "config.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/pressure/tpfa/trans_tpfa.h>
#include <opm/core/grid/GridHelpers.hpp>

namespace Dune
{
class CpGrid;
}

namespace Opm
{
namespace UgGridHelpers
{
int dimensions(const Dune::CpGrid&);

double faceArea(const Dune::CpGrid&, int);
}
}

namespace
{
const double* multiplyFaceNormalWithArea(const Dune::CpGrid& grid, int face_index, const double* in)
{
    int d=Opm::UgGridHelpers::dimensions(grid);
    double* out=new double[d];
    double area=Opm::UgGridHelpers::faceArea(grid, face_index);
    
    for(int i=0;i<d;++i)
        out[i]=in[i]*area;
    return out;
}

inline const double* multiplyFaceNormalWithArea(const UnstructuredGrid& grid, int face_index, const double* in)
{
    return in;
}

inline void maybeFreeFaceNormal(const Dune::CpGrid&, const double* array)
{
    delete[] array;
}

inline void maybeFreeFaceNormal(const UnstructuredGrid&, const double*)
{}
}

/* ---------------------------------------------------------------------- */
/* htrans <- sum(C(:,i) .* K(cellNo,:) .* N(:,j), 2) ./ sum(C.*C, 2) */
/* ---------------------------------------------------------------------- */
template<class Grid>
void
tpfa_htrans_compute(const Grid* G, const double *perm, double *htrans)
/* ---------------------------------------------------------------------- */
{
    using namespace Opm::UgGridHelpers;
    int    d, j;
    double s, dist, denom;

    double Kn[3];
    typename CellCentroidTraits<Grid>::IteratorType cc = beginCellCentroids(*G);
    typename Cell2FacesTraits<Grid>::Type c2f = cell2Faces(*G);
    typename FaceCellTraits<Grid>::Type face_cells = faceCells(*G);
    
    const double *n;
    const double *K;

    MAT_SIZE_T nrows, ncols, ldA, incx, incy;
    double a1, a2;

    d = dimensions(*G);

    nrows = ncols    = ldA = d;
    incx  = incy     = 1      ;
    a1    = 1.0;  a2 = 0.0    ;

    for (int c =0, i = 0; c < numCells(*G); c++) {
        K  = perm + (c * d * d);
        
        typedef typename Cell2FacesTraits<Grid>::Type::row_type FaceRow;
        FaceRow faces = c2f[c];
        
        for(typename FaceRow::const_iterator f=faces.begin(), end=faces.end();
            f!=end; ++f, ++i)
        {
            s = 2.0*(face_cells(*f, 0) == c) - 1.0;
            n = faceNormal(*G, *f);
            const double* nn=multiplyFaceNormalWithArea(*G, *f, n);
            const double* fc = &(faceCentroid(*G, *f)[0]);
            dgemv_("No Transpose", &nrows, &ncols,
                   &a1, K, &ldA, nn, &incx, &a2, &Kn[0], &incy);
            maybeFreeFaceNormal(*G, nn);
            
            htrans[i] = denom = 0.0;
            for (j = 0; j < d; j++) {
                dist = fc[j] - getCoordinate(cc, j);

                htrans[i] += s * dist * Kn[j];
                denom     +=     dist * dist;
            }

            assert (denom > 0);
            htrans[i] /= denom;
            htrans[i]  = fabs(htrans[i]);
        }
        // Move to next cell centroid.
        cc = increment(cc, 1, d);
    }
}


/* ---------------------------------------------------------------------- */
template<class Grid>
void
tpfa_trans_compute(const Grid* G, const double *htrans, double *trans)
/* ---------------------------------------------------------------------- */
{
    using namespace Opm::UgGridHelpers;

    for (int f = 0; f < numFaces(*G); f++) {
        trans[f] = 0.0;
    }

    typename Cell2FacesTraits<Grid>::Type c2f = cell2Faces(*G);

    for (int c = 0, i = 0; c < numCells(*G); c++) {
        typedef typename Cell2FacesTraits<Grid>::Type::row_type FaceRow;
        FaceRow faces = c2f[c];
        
        for(typename FaceRow::const_iterator f=faces.begin(), end=faces.end();
            f!=end; ++f, ++i)
        {
            trans[*f] += 1.0 / htrans[i];
        }
    }

    for (int f = 0; f < numFaces(*G); f++) {
        trans[f] = 1.0 / trans[f];
    }
}


/* ---------------------------------------------------------------------- */
 template<class Grid>
void
tpfa_eff_trans_compute(const Grid*        G,
                       const double *totmob,
                       const double *htrans,
                       double       *trans)
/* ---------------------------------------------------------------------- */
{
    using namespace Opm::UgGridHelpers;

    for (int f = 0; f < numFaces(*G); f++) {
        trans[f] = 0.0;
    }

    typename Cell2FacesTraits<Grid>::Type c2f = cell2Faces(*G);

    for (int c = 0, i = 0; c < numCells(*G); c++) {
        typedef typename Cell2FacesTraits<Grid>::Type::row_type FaceRow;
        FaceRow faces = c2f[c];
        
        for(typename FaceRow::const_iterator f=faces.begin(), end=faces.end();
            f!=end; ++f, ++i)
        {
            trans[*f] += 1.0 / (totmob[c] * htrans[i]);
        }
    }

             for (int f = 0; f < numFaces(*G); f++) {
        trans[f] = 1.0 / trans[f];
    }
}
