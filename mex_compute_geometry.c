#include <stdlib.h>

#include <mex.h>
#include "mrst_api.h"
#include "geometry.h"


static void compute_geometry(const mxArray *G, double *a, double *fc, double *fn, double *cc, double *cv)
{

   /* Grid topology: */
   const int d          = getNumberOfDimensions(G);  mxAssert(d==3, "Sorry, only support for 3D grids.");
   int      *nodepos    = getFaceNodePos(G);
   int      *facenodes  = getFaceNodes(G);
   int      *cellfaces  = getCellFaces(G);
   int      *facepos    = getCellFacePos(G);
   double   *coords     = getNodeCoordinates(G);

   /* Compute face geometry */
   const int nf = getNumberOfFaces(G);
   compute_face_geometry(d, coords, nf, nodepos, facenodes, fn, fc, a);


   /* Compute cell geometry */
   const int nc = getNumberOfCells(G);
   int i;

   compute_cell_geometry(d, coords, nf, nodepos, facenodes, fc, nc, facepos, cellfaces, cc, cv);


   /* Clean up */
   mxFree(coords); mxFree(facenodes); mxFree(nodepos); mxFree(cellfaces); mxFree(facepos);
}

void
mexFunction(int nlhs,       mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
   const mxArray *G = prhs[0];

   mxArray *fa = mxCreateDoubleMatrix(getNumberOfFaces(G),   1,                      mxREAL);
   mxArray *fc = mxCreateDoubleMatrix(getNumberOfDimensions(G), getNumberOfFaces(G), mxREAL);
   mxArray *fn = mxCreateDoubleMatrix(getNumberOfDimensions(G), getNumberOfFaces(G), mxREAL);
   mxArray *cc = mxCreateDoubleMatrix(getNumberOfDimensions(G), getNumberOfCells(G), mxREAL);
   mxArray *cv = mxCreateDoubleMatrix(getNumberOfCells(G),   1,                      mxREAL);
   
   compute_geometry(G, mxGetPr(fa), mxGetPr(fc), mxGetPr(fn), mxGetPr(cc), mxGetPr(cv));

   mxAssert(nlhs==5, "Requires 5 return arguments");
   plhs[0] = fa;
   plhs[1] = fc;
   plhs[2] = fn;
   plhs[3] = cc;
   plhs[4] = cv;
}
