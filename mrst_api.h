#ifndef MRST_API_H_INCLUDED
#define MRST_API_H_INCLUDED

/* 
 *  "API" to MRST grid : implements access to raw C vectors.
 */

int  getNumberOfDimensions (const mxArray *G);
void getLocal2GlobalCellMap(const mxArray *G);

/* Node coordinates */
int     getNumberOfNodes      (const mxArray *G);
double *getNodeCoordinates(const mxArray *G); /* copy */

/* Face topology */
int     getNumberOfFaces      (const mxArray *G);
int     getNumberOfFaceNodes  (const mxArray *G);
int    *getFaceNodePos        (const mxArray *G);    /* copy */
int    *getFaceNodes          (const mxArray *G);    /* copy */
int    *getFaceCellNeighbors  (const mxArray *G);    /* copy */

/* Face geometry */
double *getFaceAreas          (const mxArray *G);
double *getFaceNormals        (const mxArray *G);    /* copy */
double *getFaceCentroids      (const mxArray *G);    /* copy */

/* Cell topology */
int     getNumberOfCells      (const mxArray *G);
int     getNumberOfCellFaces  (const mxArray *G);
int    *getCellFacePos        (const mxArray *G);    /* copy */
int    *getCellFaces          (const mxArray *G);    /* copy */

/* Cell geometry */
double *getCellVolumes        (const mxArray *G);
double *getCellCentroids      (const mxArray *G);    /* copy */


/* Rock properties */
double *
getPermeability(const mxArray *perm, int d);    /* copy */

#endif /* MRST_API_H_INCLUDED */
