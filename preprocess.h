#ifndef PREPROCESS_H
#define PREPROCESS_H


/* Input structure holding raw cornerpoint spec. */
struct grdecl{
  int           dims[3];
  const double *coord;
  const double *zcorn;
  const int    *actnum;
};


/* Output structure holding grid topology */
struct processed_grid{
  int    number_of_faces;
  int    *face_nodes;       /* Nodes numbers of each face sequentially.           */
  int    *face_ptr;         /* Start position for each face in face_nodes.        */
  int    *face_neighbors;   /* Global cell numbers.  2 ints per face sequentially */
  
  int    number_of_nodes;      
  double *node_coordinates; /* 3 doubles per node, sequentially                   */
  
  int    number_of_cells;   /* number of active cells                             */
  int    *local_cell_index; /* Global to local map                                */
};


void process_grdecl     (const struct grdecl   *g, 
			 double                tol, 
			 struct processed_grid *out);
void free_processed_grid(struct processed_grid *g);

#endif
