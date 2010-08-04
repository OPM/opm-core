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
void    getFaceAreas          (const mxArray *G, double **v);
void    getFaceNormals        (const mxArray *G, double **v);
void    getFaceCentroids      (const mxArray *G, double **v);

/* Cell topology */
int     getNumberOfCells      (const mxArray *G);
int     getNumberOfCellFaces  (const mxArray *G);
int    *getCellFacePos        (const mxArray *G);    /* copy */
int    *getCellFaces          (const mxArray *G);    /* copy */

/* Cell geometry */
void    getCellVolumes        (const mxArray *G, double **v);
void    getCellCentroids      (const mxArray *G, double **v);

#endif /* MRST_API_H_INCLUDED */
