#include <assert.h>
#include <stdlib.h>

#include "geometry.h"
#include "cgridinterface.h"


#if 0
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
#endif

static int
fill_cell_topology(struct processed_grid  *pg,
                   struct CornerpointGrid *G )
{
    int    f, c1, c2, tag;
    size_t c, nc, nhf;

    struct UnstructuredGrid *g;

    g = (struct UnstructuredGrid *) G;

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

        /* struct CornerpointGrid member */
        G->cface_tag  = malloc(nhf * sizeof *G->cface_tag );

        if ((g->cell_faces == NULL) || (G->cface_tag == NULL)) {
            free(G->cface_tag);     G->cface_tag    = NULL;
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
                g->cell_faces[ g->cell_facepos[c1 + 1] ] = f;

                /* struct CornerpointGrid member (!) */
                G->cface_tag [ g->cell_facepos[c1 + 1] ] = tag + 1;

                g->cell_facepos[c1 + 1] += 1;
            }

            if (c2 >= 0) {
                g->cell_faces[ g->cell_facepos[c2 + 1] ] = f;

                /* struct CornerpointGrid member (!) */
                G->cface_tag [ g->cell_facepos[c2 + 1] ] = tag + 0;

                g->cell_facepos[c2 + 1] += 1;
            }
        }
    }

    return g->cell_facepos != NULL;
}

void preprocess         (const struct grdecl   *in,
                         double                tol,
                         struct CornerpointGrid *G)
{
   struct processed_grid    pg;
   struct UnstructuredGrid *base;

   base = (struct UnstructuredGrid *) G;

   process_grdecl(in, tol, &pg);

   /*
    *  General grid interface
    */
   base->dimensions = 3;

   base->number_of_nodes  = pg.number_of_nodes;
   base->number_of_faces  = pg.number_of_faces;
   base->number_of_cells  = pg.number_of_cells;

   base->node_coordinates = pg.node_coordinates;

   base->face_nodes       = pg.face_nodes;
   base->face_nodepos     = pg.face_ptr;
   base->face_cells       = pg.face_neighbors;

   base->face_centroids   = NULL;
   base->face_normals     = NULL;
   base->face_areas       = NULL;

   fill_cell_topology(&pg, G);

   base->cell_centroids   = NULL;
   base->cell_volumes     = NULL;


   /*
    *  Cornerpoint grid interface
    */
   G->cartdims[0]      = pg.dimensions[0];
   G->cartdims[1]      = pg.dimensions[1];
   G->cartdims[2]      = pg.dimensions[2];

   free(pg.face_tag);

   G->index_map = pg.local_cell_index;
}

void free_cornerpoint_grid(struct CornerpointGrid *G)
{
    free(G->grid.face_nodes);
    free(G->grid.face_nodepos);
    free(G->grid.face_cells);
    free(G->grid.cell_facepos);
    free(G->grid.cell_faces);

    free(G->grid.node_coordinates);
    free(G->grid.face_centroids);
    free(G->grid.face_areas);
    free(G->grid.face_normals);
    free(G->grid.cell_centroids);
    free(G->grid.cell_volumes);

    free(G->index_map);
    free(G->cface_tag);
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

void compute_geometry(struct CornerpointGrid *G)
{
    int ok;

    struct UnstructuredGrid *g;

    g = (struct UnstructuredGrid *) G;

    ok = allocate_geometry(g);

    if (ok) {
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
