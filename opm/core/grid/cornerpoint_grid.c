/*===========================================================================
//
// File: cgridinterface.c
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
*/

#include <assert.h>
#include <stdlib.h>

#include <opm/core/grid/cornerpoint_grid.h>
#include <opm/core/grid/cpgpreprocess/geometry.h>
#include <opm/core/grid/cpgpreprocess/preprocess.h>
#include <opm/core/grid.h>


static int
fill_cell_topology(struct processed_grid  *pg,
                   struct UnstructuredGrid *g )
{
    int    f, c1, c2, tag;
    size_t c, nc, nhf;

    nc = g->number_of_cells;

    g->cell_facepos = malloc((nc + 1) * sizeof *g->cell_facepos);

    if (g->cell_facepos != NULL) {
        /* Allocate and initialise compressed cell-to-face topology. */

        for (c = 0; c < nc + 1; c++) { g->cell_facepos[c] = 0; }

        for (f = 0; f < g->number_of_faces; f++) {
            c1 = g->face_cells[2*f + 0];
            c2 = g->face_cells[2*f + 1];

            if (c1 >= 0) { g->cell_facepos[c1 + 1] += 1; }
            if (c2 >= 0) { g->cell_facepos[c2 + 1] += 1; }
        }

        for (c = 1; c <= nc; c++) {
            g->cell_facepos[0] += g->cell_facepos[c];
            g->cell_facepos[c]  = g->cell_facepos[0] - g->cell_facepos[c];
        }

        nhf = g->cell_facepos[0];
        g->cell_facepos[0] = 0;

        g->cell_faces   = malloc(nhf * sizeof *g->cell_faces  );
        g->cell_facetag = malloc(nhf * sizeof *g->cell_facetag);

        if ((g->cell_faces == NULL) || (g->cell_facetag == NULL)) {
            free(g->cell_facetag);  g->cell_facetag = NULL;
            free(g->cell_faces);    g->cell_faces   = NULL;
            free(g->cell_facepos);  g->cell_facepos = NULL;
        }
    }

    if (g->cell_facepos != NULL) {
        /* Compute final cell-to-face mapping and half-face tags.
         *
         * Process relies on preprocess() producing neighbourship
         * definitions for which the normals point (logically) in the
         * positive I,J,K directions *and* from ->face_cells[2*f+0] to
         * ->face_cells[2*f+1] (i.e., from first to second cell of
         * interface 'f'--be it internal or outer).
         *
         * For instance, a "LEFT" connection (pg->face_tag==LEFT==0)
         * for which the normal points into the cell (i.e., when
         * ->face_cells[2*f+1] >= 0), is a half-face of type 0.
         *
         * Simlarly, a "TOP" connection (pg->face_tag==TOP==2) for
         * which the normal points out of the cell (i.e., when
         * ->face_cells[2*f+0] >= 0), is a half-face of type 5. */

        for (f = 0; f < g->number_of_faces; f++) {
            tag = 2 * pg->face_tag[f];    /* [0, 2, 4] */
            c1  = g->face_cells[2*f + 0];
            c2  = g->face_cells[2*f + 1];

            if (c1 >= 0) {
                g->cell_faces   [ g->cell_facepos[c1 + 1] ] = f;
                g->cell_facetag [ g->cell_facepos[c1 + 1] ] = tag + 1;

                g->cell_facepos[c1 + 1] += 1;
            }

            if (c2 >= 0) {
                g->cell_faces   [ g->cell_facepos[c2 + 1] ] = f;
                g->cell_facetag [ g->cell_facepos[c2 + 1] ] = tag + 0;

                g->cell_facepos[c2 + 1] += 1;
            }
        }
    }

    return g->cell_facepos != NULL;
}

static int
allocate_geometry(struct UnstructuredGrid *g)
{
    int ok;
    size_t nc, nf, nd;

    assert (g->dimensions == 3);

    nc = g->number_of_cells;
    nf = g->number_of_faces;
    nd = 3;

    g->face_areas     = malloc(nf * 1  * sizeof *g->face_areas);
    g->face_centroids = malloc(nf * nd * sizeof *g->face_centroids);
    g->face_normals   = malloc(nf * nd * sizeof *g->face_normals);

    g->cell_volumes   = malloc(nc * 1  * sizeof *g->cell_volumes);
    g->cell_centroids = malloc(nc * nd * sizeof *g->cell_centroids);

    ok  = g->face_areas     != NULL;
    ok += g->face_centroids != NULL;
    ok += g->face_normals   != NULL;

    ok += g->cell_volumes   != NULL;
    ok += g->cell_centroids != NULL;

    return ok == 5;
}


void compute_geometry(struct UnstructuredGrid *g)
{
    assert (g != NULL);
    if (g!=NULL)
    {
        assert (g->dimensions == 3);

        assert (g->face_centroids != NULL);
        assert (g->face_normals   != NULL);
        assert (g->face_areas     != NULL);
        assert (g->cell_centroids != NULL);
        assert (g->cell_volumes   != NULL);

        compute_face_geometry(g->dimensions  , g->node_coordinates,
                              g->number_of_faces, g->face_nodepos,
                              g->face_nodes, g->face_normals,
                              g->face_centroids, g->face_areas);

        compute_cell_geometry(g->dimensions, g->node_coordinates,
                              g->face_nodepos, g->face_nodes,
                              g->face_cells, g->face_normals,
                              g->face_centroids, g->number_of_cells,
                              g->cell_facepos, g->cell_faces,
                              g->cell_centroids, g->cell_volumes);
    }
}


struct UnstructuredGrid *
preprocess (const struct grdecl *in, double tol)
{
    struct UnstructuredGrid *g;
   int                      ok;
   struct processed_grid    pg;

   g = malloc(1 * sizeof *g);
   if (g == NULL)
   {
       return NULL;
   }

   process_grdecl(in, tol, &pg);

   /*
    *  Convert "struct processed_grid" to "struct UnstructuredGrid".
    *
    *  In particular, convey resource ownership from 'pg' to 'g'.
    *  Consequently, memory resources obtained in process_grdecl()
    *  will be released in free_grid().
    */
   g->dimensions = 3;

   g->number_of_nodes  = pg.number_of_nodes;
   g->number_of_faces  = pg.number_of_faces;
   g->number_of_cells  = pg.number_of_cells;

   g->node_coordinates = pg.node_coordinates;

   g->face_nodes       = pg.face_nodes;
   g->face_nodepos     = pg.face_ptr;
   g->face_cells       = pg.face_neighbors;

   /* Explicitly relinquish resource references conveyed to 'g'.  This
    * is needed to avoid creating dangling references in the
    * free_processed_grid() call. */
   pg.node_coordinates = NULL;
   pg.face_nodes       = NULL;
   pg.face_ptr         = NULL;
   pg.face_neighbors   = NULL;

   /* Initialise subsequently allocated fields to a defined state lest
    * we free() random pointers in free_grid() if either of the
    * fill_cell_topology() or allocate_geometry() functions fail. */
   g->face_centroids   = NULL;
   g->face_normals     = NULL;
   g->face_areas       = NULL;

   g->cell_centroids   = NULL;
   g->cell_volumes     = NULL;

   g->global_cell      = NULL;

   /* allocate and fill g->cell_faces/g->cell_facepos and
    * g->cell_facetag as well as the geometry-related fields. */
   ok =       fill_cell_topology(&pg, g);
   ok = ok && allocate_geometry(g);

   if (!ok)
   {
       free_grid(g);
       g = NULL;
   }
   else
   {

       compute_geometry(g);

       g->cartdims[0]      = pg.dimensions[0];
       g->cartdims[1]      = pg.dimensions[1];
       g->cartdims[2]      = pg.dimensions[2];

       g->global_cell      = pg.local_cell_index;

       /* Explicitly relinquish resource references conveyed to 'g'.
        * This is needed to avoid creating dangling references in the
        * free_processed_grid() call. */
       pg.local_cell_index = NULL;
   }

   free_processed_grid(&pg);

   return g;
}
