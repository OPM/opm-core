/* "API" */
#include <stdlib.h>
#include <mex.h>
#include "mrst_api.h"

#define Assert mxAssert

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
   int n = mxGetNumberOfElements(a);
   int *q         = malloc(n * sizeof *q);
   if (q != NULL)
   {
      if (mxIsInt32(a))
      {
         int *p = mxGetData(a);
         int i,j;
         for (i=0; i<n; ++i)
         {
            mxAssert (p[i] <= INT_MAX,
                      "Matrix entry exceeds INT_MAX");
            q[i] = p[i]-1;
         }      
      }
      else if(mxIsDouble(a))
      {
         double *p = mxGetPr(a);
         int i,j;
         for (i=0; i<n; ++i)
         {
            mxAssert (p[i] <= INT_MAX,
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
   int *q = extractIntMatrix(a);
   int  M = mxGetM(a);
   int  N = mxGetN(a);

   int i,j;
   for(i=0; i<M; ++i)
   {
      for(j=0; j<i; ++j)
      {
         int tmp  = q[i+M*j];
         q[i+M*j] = q[j+N*i];
         q[j+N*i] = tmp;
      }
   }
   return q;
}

/* ------------------------------------------------------------------ */
int  getNumberOfDimensions(const mxArray *G)
/* ------------------------------------------------------------------ */
{   
   return mxGetN(getField(G, "nodes", "coords"));
}

/* ------------------------------------------------------------------ */
void getLocal2GlobalCellMap(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   Assert(0, "Not implemented!");
}

/* ------------------------------------------------------------------ */
int getNumberOfNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetM(getField(G, "nodes", "coords"));
}

/* ------------------------------------------------------------------ */
double *getNodeCoordinates(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   mxArray *p1, *p2;

   p1 = mxGetField(G , 0, "nodes" );
   p2 = mxGetField(p1, 0, "coords");
   
   const int n = getNumberOfNodes(G);
   const int d = getNumberOfDimensions(G);

   double *v = malloc(n * d * sizeof *v);
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
}

/* ------------------------------------------------------------------ */
int  getNumberOfFaces(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return mxGetNumberOfElements(getField(G, "faces", "nodePos"))-1;
}

/* ------------------------------------------------------------------ */
int *getFaceNodePos(const mxArray *G)
/* ------------------------------------------------------------------ */
{
   return extractIntMatrix(getField(G, "faces", "nodePos"));
}

/* ------------------------------------------------------------------ */
int getNumberOfFaceNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
     return mxGetNumberOfElements(getField(G, "faces", "nodes"));
}

/* ------------------------------------------------------------------ */
int *getFaceNodes(const mxArray *G)
/* ------------------------------------------------------------------ */
{
     return extractIntMatrix(getField(G, "faces", "nodes"));
}


/* ------------------------------------------------------------------ */
int *getFaceCellNeighbors(const mxArray *G)
/* ------------------------------------------------------------------ */
{
     return extractIntMatrixTranspose(getField(G, "faces", "neighbors"));
}

/* ------------------------------------------------------------------ */
void getFaceAreas(const mxArray *G, double **v)
/* ------------------------------------------------------------------ */
{
     *v = mxGetPr(getField(G, "faces", "areas"));
}

/* ------------------------------------------------------------------ */
void getFaceNormals(const mxArray *G, double **v)
/* ------------------------------------------------------------------ */
{
     *v = mxGetPr(getField(G, "faces", "normals"));
}

/* ------------------------------------------------------------------ */
void getFaceCentroids(const mxArray *G, double **v)
/* ------------------------------------------------------------------ */
{
     *v = mxGetPr(getField(G, "faces", "centroids"));
}

/* ------------------------------------------------------------------ */
int  getNumberOfCells(const mxArray *G)
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
void getCellVolumes(const mxArray *G, double **v)
/* ------------------------------------------------------------------ */
{
     *v = mxGetPr(getField(G, "cells", "volumes"));
}



/* ------------------------------------------------------------------ */
void getCellCentroids(const mxArray *G, double **v)
/* ------------------------------------------------------------------ */
{
     *v = mxGetPr(getField(G, "cells", "centroids"));
}

