/* "API" */
#include <mex.h>
#include "mrst_api.h"

/* ------------------------------------------------------------------ */
static mxArray*
getField(const mxArray *a, char *field, char *subfield)
/* ------------------------------------------------------------------ */
{
   if (subfield)
   {
      return mxGetField(mxGetField(a , 0, field), 0, subfield);
   }
   else
   {
      return mxGetField(a , 0, field);
   }
}

/* ------------------------------------------------------------------ */
static int *
extractIntMatrix(const mxArray *a)
/* ------------------------------------------------------------------ */
{
   int n  = mxGetNumberOfElements(a);
   int *q = mxMalloc(n * sizeof *q);
   if (q != NULL)
   {
      if (mxIsInt32(a))
      {
         int *p = mxGetData(a);
         int i,j;
         for (i=0; i<n; ++i)
         {
            q[i] = p[i]-1;
         }      
      }
      else if(mxIsDouble(a))
      {
         double *p = mxGetPr(a);
         int i,j;
         for (i=0; i<n; ++i)
         {
            mxAssert ((1 <= p[i]) && (p[i] <= INT_MAX),
                      "Matrix entry exceeds INT_MAX");
            q[i] = p[i]-1;
         }      
      }
   }
   
   return q;
}

/* ------------------------------------------------------------------ */
static int *
extractIntMatrixTranspose(const mxArray *a)
/* ------------------------------------------------------------------ */
{
   int *p = extractIntMatrix(a);
   int  M = mxGetM(a);
   int  N = mxGetN(a);
   
   int *q = mxMalloc(M * N * sizeof *q);
   if (q != NULL)
   {
      int i,j;
      for(i=0; i<M; ++i)
      {
         for(j=0; j<N; ++j)
         {
            q[i*N+j] = p[i+M*j];
         }
      }
   }
   mxFree(p);
   return q;
}

/* ------------------------------------------------------------------ */
static double *
extractDoubleMatrixTranspose(const mxArray *a)
/* ------------------------------------------------------------------ */
{
   int  M    = mxGetM(a);
   int  N    = mxGetN(a);
   double *q = mxMalloc(M * N * sizeof *q);
   if (q != NULL)
   {
      double *p = mxGetPr(a);
      int i,j;
      for(i=0; i<M; ++i)
      {
         for(j=0; j<N; ++j)
         {
            q[i*N+j] = p[i+M*j];
         }
      }
   }
   return q;
}

/* ------------------------------------------------------------------ */
int  
getNumberOfDimensions(const mxArray *G)
/* ------------------------------------------------------------------ */
{   
   return mxGetN(getField(G, "nodes", "coords"));
}

/* ------------------------------------------------------------------ */
void 
getLocal2GlobalCellMap(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   mxAssert(0, "Not implemented!");
}

/* ------------------------------------------------------------------ */
int 
getNumberOfNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetM(getField(G, "nodes", "coords"));
}

/* ------------------------------------------------------------------ */
double *
getNodeCoordinates(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractDoubleMatrixTranspose(getField(G, "nodes", "coords"));
#if 0
   mxArray *p1, *p2;

   p1 = mxGetField(G , 0, "nodes" );
   p2 = mxGetField(p1, 0, "coords");
   
   const int n = getNumberOfNodes(G);
   const int d = getNumberOfDimensions(G);

   double *v = mxMalloc(n * d * sizeof *v);
   if (v != NULL)
   {

      double *tmp = mxGetPr(p2);
      int i,j;
      for (i=0; i<n; ++i)
      {
         for(j=0; j<d; ++j)
         {
            v[d*i+j] = tmp[i + n*j];
         }
      }
   }
   return v;
#endif
}

/* ------------------------------------------------------------------ */
int  
getNumberOfFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetNumberOfElements(getField(G, "faces", "nodePos"))-1;
}

/* ------------------------------------------------------------------ */
int *
getFaceNodePos(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractIntMatrix(getField(G, "faces", "nodePos"));
}

/* ------------------------------------------------------------------ */
int 
getNumberOfFaceNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
     return mxGetNumberOfElements(getField(G, "faces", "nodes"));
}

/* ------------------------------------------------------------------ */
int *
getFaceNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
     return extractIntMatrix(getField(G, "faces", "nodes"));
}


/* ------------------------------------------------------------------ */
int *
getFaceCellNeighbors(const mxArray *G)
/* ------------------------------------------------------------------ */
{
     return extractIntMatrixTranspose(getField(G, "faces", "neighbors"));
}

/* ------------------------------------------------------------------ */
double * 
getFaceAreas(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetPr(getField(G, "faces", "areas"));
}

/* ------------------------------------------------------------------ */
double *
getFaceNormals(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractDoubleMatrixTranspose(getField(G, "faces", "normals"));
}

/* ------------------------------------------------------------------ */
double *
getFaceCentroids(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractDoubleMatrixTranspose(getField(G, "faces", "centroids"));
}

/* ------------------------------------------------------------------ */
int  
getNumberOfCells(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetNumberOfElements(getField(G, "cells", "facePos"))-1;   
}

/* ------------------------------------------------------------------ */
int *getCellFacePos(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractIntMatrix(getField(G, "cells", "facePos"));
}

/* ------------------------------------------------------------------ */
int getNumberOfCellFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetNumberOfElements(getField(G, "cells", "faces"))-1;   
}

/* ------------------------------------------------------------------ */
int *getCellFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractIntMatrix(getField(G, "cells", "faces"));
}

/* ------------------------------------------------------------------ */
double *
getCellVolumes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetPr(getField(G, "cells", "volumes"));
}



/* ------------------------------------------------------------------ */
double *
getCellCentroids(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractDoubleMatrixTranspose(getField(G, "cells", "centroids"));
}



/*
 *
 *
 *
 */


/* ------------------------------------------------------------------ */
static int
allocate_perm_data(int ncells, int d, double **K)
/* ------------------------------------------------------------------ */
{
    int ret;
    double *perm;

    perm = mxMalloc(ncells * d * d * sizeof *perm);

    if (perm != NULL) {
        *K = perm;
        ret = 1;
    } else {
        ret = 0;
    }

    return ret;
}


/* ------------------------------------------------------------------ */
double *
getPermeability(const mxArray *perm, int d)
/* ------------------------------------------------------------------ */
{
    int ncells, ncomp, alloc_ok;

    int c, i, off;

    double *k, *tensor;

    ncells = mxGetM(perm);
    ncomp  = mxGetN(perm);

    alloc_ok = allocate_perm_data(ncells, d, &k);

    if (alloc_ok) {
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
           mxAssert (ncomp == 3, "");

            for (c = 0; c < ncells; c++) {
                off = c * d * d;
                k[0 + off] = tensor[c + 0*ncells];
                k[1 + off] = tensor[c + 1*ncells];
                k[2 + off] = tensor[c + 1*ncells];
                k[3 + off] = tensor[c + 2*ncells];
            }
        } else {
            /* Full 3D tensor */
            mxAssert ((d == 3) && (ncomp == 6), "");

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
    return k;
}
