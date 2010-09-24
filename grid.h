#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED


#ifdef __cplusplus
extern "C" {
#endif


/*   GRID_TOPOLOGY and GRID_GEOMETRY must be at the beginning of every
 *   grid type.  
 *
 *
 */
#define GRID_TOPOLOGY                           \
   int    dimensions;                           \
                                                \
   int    number_of_cells;                      \
   int    number_of_faces;                      \
   int    number_of_nodes;                      \
                                                \
   int    *face_nodes;                          \
   int    *face_nodepos;                        \
   int    *face_cells;                          \
                                                \
   int    *cell_faces;                          \
   int    *cell_facepos;                        \
   
#define GRID_GEOMETRY                           \
   double *node_coordinates;                    \
                                                \
   double *face_centroids;                      \
   double *face_areas;                          \
   double *face_normals;                        \
                                                \
   double *cell_centroids;                      \
   double *cell_volumes;                        \
   
   
typedef struct {
   GRID_TOPOLOGY
   GRID_GEOMETRY
} grid_t;


void free_grid           (grid_t *g);
void alloc_grid_geometry (grid_t *g);
void print_grid_summary  (grid_t *g);

#ifdef __cplusplus
}
#endif

#endif /* GRID_H_INCLUDED */
