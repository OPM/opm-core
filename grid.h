/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

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
