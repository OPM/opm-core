#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <mex.h>

#include "preprocess.h"
#include "mxgrdecl.h"


void fill_grid(mxArray **out, struct processed_grid *grid)
{
  const char *names[] = {"nodes", "faces", "cells", "cellFaces", "faceNodes"};
  mxArray *G = mxCreateStructMatrix(1,1,5,names);

  int i,j;
  double *ptr;


  /* nodes */
  const char *n2[] = {"num", "coords"};
  mxArray *nodes  = mxCreateStructMatrix(1,1,2,n2);
  mxSetField(nodes, 0, "num", mxCreateDoubleScalar(grid->number_of_nodes));

  mxArray *nodecoords = mxCreateDoubleMatrix(grid->number_of_nodes, 3, mxREAL);
  ptr = mxGetPr(nodecoords);
  
  for (j=0;j<3;++j){
    for(i=0; i<grid->number_of_nodes; ++i){
      ptr[i+grid->number_of_nodes*j] = grid->node_coordinates[3*i+j];
    }
  }
  mxSetField(nodes, 0, "coords", nodecoords);
  mxSetField(G, 0, "nodes", nodes);


  /* faces */
  const char *n3[] = {"num", "neighbors", "numNodes", "nodePos", "tag"};
  mxArray *faces = mxCreateStructMatrix(1,1,5,n3);


  mxSetField(faces, 0, "num", mxCreateDoubleScalar(grid->number_of_faces));

  mxArray *faceneighbors = mxCreateDoubleMatrix(grid->number_of_faces, 2, mxREAL);
  /* int i,j; */
  ptr = mxGetPr(faceneighbors);
  for(j=0; j<2; ++j){
    for (i=0; i<grid->number_of_faces; ++i){
      int ix = grid->face_neighbors[2*i+j];
      if (ix == -1){
	ptr[i+grid->number_of_faces*j] = 0;
      }else{
	ptr[i+grid->number_of_faces*j] = ix+1;/* grid->local_cell_index[ix]+1; */
      }
    }
  }
  mxSetField(faces, 0, "neighbors", faceneighbors);

  mxArray *numnodes = mxCreateDoubleMatrix(grid->number_of_faces,   1, mxREAL);
  mxArray *nodepos  = mxCreateDoubleMatrix(grid->number_of_faces+1, 1, mxREAL);
  double *ptr2 = mxGetPr(nodepos);
  ptr2[0] = 1;
  ptr = mxGetPr(numnodes);
  for (i=0; i<grid->number_of_faces; ++i){
    ptr[i] = grid->face_ptr[i+1]-grid->face_ptr[i];
    ptr2[i+1] = ptr2[i] + ptr[i]; 
  }
  mxSetField(faces, 0, "numNodes", numnodes);
  mxSetField(faces, 0, "nodePos",  nodepos);

  mxArray *tags = mxCreateDoubleMatrix(grid->number_of_faces, 1, mxREAL);
  ptr = mxGetPr(tags);
  for (i = 0; i < grid->number_of_faces; ++i) {
      ptr[i] = grid->face_tag[i] + 1;
  }
  mxSetField(faces, 0, "tag", tags);

  mxSetField(G, 0, "faces", faces);
  

  const char *n4[] = {"num", "numFaces", "facePos", "indexMap"};
  mxArray *cells = mxCreateStructMatrix(1,1,3,n4);

  mxSetField(cells, 0, "num", mxCreateDoubleScalar(grid->number_of_cells));
  mxArray *map = mxCreateDoubleMatrix(grid->number_of_cells, 1, mxREAL);
  ptr = mxGetPr(map);
  for(i=0; i<grid->number_of_cells; ++i){
    ptr[i] = grid->local_cell_index[i]+1;
  }
  mxSetField(cells, 0, "indexMap", map);
    

  mxArray *numfaces = mxCreateDoubleMatrix(grid->number_of_cells, 1, mxREAL);
  mxArray *facepos = mxCreateDoubleMatrix(grid->number_of_cells+1, 1, mxREAL);
  ptr = mxGetPr(numfaces);
  ptr2 = mxGetPr(facepos);
  ptr2[0] = 1;
  for(i=0; i<grid->number_of_cells; ++i){ 
    ptr[i] = 0.0;
  }
  for (i=0; i<2*grid->number_of_faces; ++i){
    int c=grid->face_neighbors[i];
    if(c != -1) {
      ptr[c]++;
    }
  }
  for(i=0; i<grid->number_of_cells; ++i){ 
    ptr2[i+1] = ptr2[i] + ptr[i];
  }
  
  mxSetField(cells, 0, "numFaces", numfaces);
  mxSetField(cells, 0, "facePos", facepos);

  mxSetField(G, 0, "cells", cells);



  int *counter = calloc(grid->number_of_cells, sizeof(*counter));
  int num_half_faces = 0;
  for(i=0; i<grid->number_of_cells; ++i){
    counter[i] = num_half_faces;
    num_half_faces += ptr[i];
  }
  
  mxArray *cellfaces = mxCreateDoubleMatrix(num_half_faces, 1, mxREAL);
  ptr = mxGetPr(cellfaces);


  for (i=0; i<grid->number_of_faces; ++i){
    int c1 = grid->face_neighbors[2*i];
    int c2 = grid->face_neighbors[2*i+1];
    if(c1 != -1) ptr[counter[c1]++] = i+1;
    if(c2 != -1) ptr[counter[c2]++] = i+1;
  }

  mxSetField(G, 0, "cellFaces", cellfaces);


  int n = grid->face_ptr[grid->number_of_faces];

  mxArray *facenodes = mxCreateDoubleMatrix(n, 1, mxREAL);
  ptr = mxGetPr(facenodes);
  for (i=0; i<n; ++i) {
    ptr[i] = grid->face_nodes[i]+1;
  }

  mxSetField(G, 0, "faceNodes", facenodes);
  free(counter);


  out[0] = G;

}


/* Gateway routine for Matlab mex function.              */
/*-------------------------------------------------------*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Set up data passed from Matlab */
  struct grdecl g;
  struct processed_grid out;
  double tolerance = 0.0;

  mx_init_grdecl(&g, prhs[0]);
  if (nrhs == 2){
    
    tolerance = mxGetScalar (prhs[1]);
  }
  
  process_grdecl(&g, tolerance, &out);


  if (plhs >0){
    /* write to matlab struct */
    fill_grid(plhs, &out);
  }
  

  /* Free whatever was allocated in initGrdecl. */
  free_processed_grid(&out);
}
