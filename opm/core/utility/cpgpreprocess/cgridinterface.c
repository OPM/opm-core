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

#include <opm/core/utility/cpgpreprocess/cgridinterface.h>
#include <opm/core/utility/cpgpreprocess/geometry.h>
#include <opm/core/utility/cpgpreprocess/preprocess.h>
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

        g->cell_faces = malloc(nhf * sizeof *g->cell_faces);

        /* struct UnstructuredGrid member */
        g->cell_facetag  = malloc(nhf * sizeof *g->cell_facetag );
        
        if ((g->cell_faces == NULL) || (g->cell_facetag == NULL)) {
            free(g->cell_facetag);  g->cell_facetag = NULL;
            free(g->cell_faces);    g->cell_faces   = NULL;
            free(g->cell_facepos);  g->cell_facepos = NULL;
        }
    }

    if (g->cell_facepos != NULL) {
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

void free_grid(struct UnstructuredGrid *g)
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
    */
   g->dimensions = 3;

   g->number_of_nodes  = pg.number_of_nodes;
   g->number_of_faces  = pg.number_of_faces;
   g->number_of_cells  = pg.number_of_cells;

   g->node_coordinates = pg.node_coordinates;

   g->face_nodes       = pg.face_nodes;
   g->face_nodepos     = pg.face_ptr;
   g->face_cells       = pg.face_neighbors;

   g->face_centroids   = NULL;
   g->face_normals     = NULL;
   g->face_areas       = NULL;

   g->cell_centroids   = NULL;
   g->cell_volumes     = NULL;

   /* Put ->global_cell in defined and harmless state to prevent
    * freeing a random pointer in case of failing to allocate geometry
    * resources. */
   g->global_cell = NULL;

   /* allocate and fill g->cell_faces/g->cell_facepos and
    * g->cell_facetag */
   fill_cell_topology(&pg, g);

   ok = allocate_geometry(g);
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
   }

   return g;
}
