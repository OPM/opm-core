#ifndef MIMETIC_GEOMETRY_H_INCLUDED
#define MIMETIC_GEOMETRY_H_INCLUDED

void compute_face_geometry(int ndims, double *coords, int nfaces, 
                           int *nodepos, int *facenodes, 
                           double *fnormals, double *fcentroids, 
                           double *fareas);
void compute_cell_geometry(int ndims, double *coords, int nfaces,
                           int *nodepos, int *facenodes, 
                           double *fcentroids, int ncells, 
                           int *facepos, int *cellfaces, 
                           double *ccentroids, double *cvolumes);

#endif /* MIMETIC_GEOMETRY_H_INCLUDED */
