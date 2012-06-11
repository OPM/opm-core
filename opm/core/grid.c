/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/grid.h>
#include <stdlib.h>


void
destroy_grid(struct UnstructuredGrid *g)
{
    if (g!=NULL)
    {
	free(g->face_nodes);
	free(g->face_nodepos);
	free(g->face_cells);
	free(g->cell_facepos);
	free(g->cell_faces);

	free(g->node_coordinates);
	free(g->face_centroids);
	free(g->face_areas);
	free(g->face_normals);
	free(g->cell_centroids);
	free(g->cell_volumes);

	free(g->global_cell);
	free(g->cell_facetag);
    }

    free(g);
}


struct UnstructuredGrid *
create_grid_empty(void)
{
    struct UnstructuredGrid *G, g = { 0 };

    G = malloc(1 * sizeof *G);

    if (G != NULL) {
	*G = g;
    }

    return G;
}


struct UnstructuredGrid *
allocate_grid(size_t ndims     ,
	      size_t ncells    ,
	      size_t nfaces    ,
	      size_t nfacenodes,
	      size_t ncellfaces,
	      size_t nnodes    )
{
    size_t nel;
    struct UnstructuredGrid *G;

    G = create_grid_empty();

    if (G != NULL) {
	/* Node fields ---------------------------------------- */
	nel                 = nnodes * ndims;
	G->node_coordinates = malloc(nel * sizeof *G->node_coordinates);

	/* Face fields ---------------------------------------- */
	nel               = nfacenodes;
	G->face_nodes     = malloc(nel * sizeof *G->face_nodes);

	nel               = nfaces + 1;
	G->face_nodepos   = malloc(nel * sizeof *G->face_nodepos);

	nel               = 2 * nfaces;
	G->face_cells     = malloc(nel * sizeof *G->face_cells);

	nel               = nfaces * ndims;
	G->face_centroids = malloc(nel * sizeof *G->face_centroids);

	nel               = nfaces * ndims;
	G->face_normals   = malloc(nel * sizeof *G->face_normals);

	nel               = nfaces * 1;
	G->face_areas     = malloc(nel * sizeof *G->face_areas);


	/* Cell fields ---------------------------------------- */
	nel               = ncellfaces;
	G->cell_faces     = malloc(nel * sizeof *G->cell_faces);

	nel               = ncells + 1;
	G->cell_facepos   = malloc(nel * sizeof *G->cell_facepos);

	nel               = ncells * ndims;
	G->cell_centroids = malloc(nel * sizeof *G->cell_centroids);

	nel               = ncells * 1;
	G->cell_volumes   = malloc(nel * sizeof *G->cell_volumes);

	if ((G->node_coordinates == NULL) ||
	    (G->face_nodes       == NULL) ||
	    (G->face_nodepos     == NULL) ||
	    (G->face_cells       == NULL) ||
	    (G->face_centroids   == NULL) ||
	    (G->face_normals     == NULL) ||
	    (G->face_areas       == NULL) ||
	    (G->cell_faces       == NULL) ||
	    (G->cell_facepos     == NULL) ||
	    (G->cell_centroids   == NULL) ||
	    (G->cell_volumes     == NULL)  )
	    {
		destroy_grid(G);
		G = NULL;
	    }
    }

    return G;
}
