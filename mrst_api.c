/* "API" */
#include <stdlib.h>
#include <mex.h>
#include "mrst_api.h"


int  getNumberOfDimensions(const mxArray *G)
{   
   mxArray *p1, *p2;

   p1 = mxGetField(G , 0, "nodes" );
   p2 = mxGetField(p1, 0, "coords");
   
   return mxGetN(p2);
}

void getLocal2GlobalCellMap(const mxArray *G)
{
}

int getNumberOfNodes(const mxArray *G)
{
   mxArray *p1, *p2;

   p1 = mxGetField(G , 0, "nodes" );
   p2 = mxGetField(p1, 0, "coords");
   
   return mxGetM(p2);
}

double *getNodeCoordinates(const mxArray *G)
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

int  getNumberOfFaces(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "nodePos");
   
     return mxGetNumberOfElements(p2)-1;
}

int *getFaceNodePos(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "nodePos");
   
     mxAssert (mxIsInt32(p2),
               "faces.nodePos must be INT32");
     int *nodepos = (int*) mxGetData(p2);

     int n  = getNumberOfFaces(G)+1;
     int *v = malloc(n * sizeof *v);
     if (v != NULL)
     {
        int i;
        for (i=0; i<n; ++i)
        {
           v[i] = nodepos[i]-1;
        }
     }
     return v;
}

int getNumberOfFaceNodes(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "nodes");
     return mxGetNumberOfElements(p2);
}

int *getFaceNodes(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "nodes");

     mxAssert (mxIsInt32(p2),
               "faces.nodes must be INT32");

     int *facenodes = (int*) mxGetData(p2);
     int n          = getNumberOfFaceNodes(G);
     int *v         = malloc(n * sizeof *v);
     if (v != NULL)
     {
        int i;
        for (i=0; i<n; ++i)
        {
           v[i] = facenodes[i]-1;
        }
     }
     return v;
}

int *getFaceCellNeighbors(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "neighbors");

     mxAssert (mxIsInt32(p2),
               "faces.neighbors must be INT32");
     int *neighbors = (int*) mxGetData(p2);
     int n  = 2*getNumberOfFaces(G);
     int *v = malloc(n * sizeof *v);
     if (v != NULL)
     {
        int i;
        for (i=0; i<n; ++i)
        {
           v[i] = neighbors[i]-1;
        }
     }
     return v;
}

void getFaceAreas(const mxArray *G, double **v)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "areas");
   
     *v = mxGetPr(p2);
}

void getFaceNormals(const mxArray *G, double **v)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "normals");
   
     *v = mxGetPr(p2);
}

void getFaceCentroids(const mxArray *G, double **v)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "faces" );
     p2 = mxGetField(p1, 0, "centroids");
   
     *v = mxGetPr(p2);
}

int  getNumberOfCells(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "cells" );
     p2 = mxGetField(p1, 0, "facePos");
   
     return mxGetNumberOfElements(p2)-1;   
}

int *getCellFacePos(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "cells" );
     p2 = mxGetField(p1, 0, "facePos");
   
     mxAssert (mxIsInt32(p2),
               "cells.facePos must be INT32");
     int *facepos = (int*) mxGetData(p2);
     int n = getNumberOfCells(G)+1;
     int *v = malloc(n * sizeof *v);
     if (v != NULL)
     {
        int i;
        for (i=0; i<n; ++i)
        {
           v[i] = facepos[i]-1;
        }
     }
     return v;
}

int getNumberOfCellFaces(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "cells" );
     p2 = mxGetField(p1, 0, "faces");
   
     return mxGetNumberOfElements(p2)-1;   
}

int *getCellFaces(const mxArray *G)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "cells" );
     p2 = mxGetField(p1, 0, "faces");
   
     mxAssert (mxIsInt32(p2),
               "cells.faces must be INT32");
     int *cellfaces = (int*) mxGetData(p2);
     int n = getNumberOfCellFaces(G);
     int *v = malloc(n * sizeof *v);
     if (v != NULL)
     {
        int i;
        for (i=0; i<n; ++i)
        {
           v[i] = cellfaces[i]-1;
        }
     }
     return v;
}

void getCellVolumes(const mxArray *G, double **v)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "cells" );
     p2 = mxGetField(p1, 0, "volumes");
   
     *v = mxGetPr(p2);
}



void getCellCentroids(const mxArray *G, double **v)
{
     mxArray *p1, *p2;

     p1 = mxGetField(G , 0, "cells" );
     p2 = mxGetField(p1, 0, "centroids");
   
     *v = mxGetPr(p2);
}

