#include <stdlib.h>
#include "newinterface.h"

   

static int *compute_cell_facepos(grid_t *g)
{
   int i,j,k;
   int *facepos = malloc((g->number_of_cells + 1) * sizeof *facepos);
   int *fcells  = g->face_cells;
   
   for (i=0; i<g->number_of_cells; ++i) {
      facepos [i] = 0;
   }
   
   for (i=0; i<2*g->number_of_faces; ++i) {
         if (*fcells != -1) {
            (facepos[*fcells])++;
         }
         fcells++;
   }
   
   /* cumsum */
   j=0;
   for (i=0; i<g->number_of_cells; ++i) {
      k = j + facepos[i];
      facepos[i] = j;
      j = k;
   }
   facepos[i] = j;

   return facepos;
}


static int *compute_cell_faces(grid_t *g)
{
   int *cfaces = malloc(g->cell_facepos[g->number_of_cells] * sizeof *cfaces);
   int *work   = malloc(g->number_of_cells * sizeof *work);
   int *fcells = g->face_cells;
   int i,k,cell;
   for(i=0; i<g->number_of_cells; ++i) {
      work[i] = 0;
   }
  
   for (i=0; i<g->number_of_faces; ++i) {
      for (k=0;k<2; ++k) {
         
         if (*fcells != -1) {
            cell = *fcells;
            cfaces[g->cell_facepos[cell] + work[cell]] = i;
            work[cell]++;
         }
         fcells++;
      }
   }
   free(work);   

   return cfaces;
}

void preprocess         (const struct grdecl   *in, 
                         double                tol, 
                         cornerpoint_grid_t    *out)
{
   struct processed_grid pg;
   process_grdecl(in, tol, &pg);

   /* 
    *  General grid interface 
    */
   out->dimensions = 3;

   out->number_of_nodes  = pg.number_of_nodes;
   out->number_of_faces  = pg.number_of_faces;
   out->number_of_cells  = pg.number_of_cells;

   out->node_coordinates = pg.node_coordinates;

   out->face_nodes       = pg.face_nodes;
   out->face_nodepos     = pg.face_ptr;
   out->face_cells       = pg.face_neighbors;
   
   out->face_centroids   = NULL;
   out->face_normals     = NULL;
   out->face_areas       = NULL;
   
   /* NB: compute_cell_facepos must be called before compute_cell_faces */
   out->cell_facepos     = compute_cell_facepos((grid_t*) out);
   out->cell_faces       = compute_cell_faces  ((grid_t*) out);
   out->cell_centroids   = NULL;
   out->cell_volumes     = NULL;


   /* 
    *  Cornerpoint grid interface 
    */
   out->cartdims[0]      = pg.dimensions[0];
   out->cartdims[1]      = pg.dimensions[1];
   out->cartdims[2]      = pg.dimensions[2];

   out->face_tag         = pg.face_tag;
   out->number_of_nodes_on_pillars = pg.number_of_nodes_on_pillars;
   out->cartesian_cell_index = pg.local_cell_index;
}
   
void free_cornerpoint_grid(cornerpoint_grid_t    *g)
{
   free_grid((grid_t*) g);

   free(g->face_tag);
   free(g->cartesian_cell_index);
}

