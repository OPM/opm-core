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

#include "config.h"
#include <opm/core/grid.h>

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


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
        /* Grid fields ---------------------------------------- */
        G->dimensions       = ndims;
        G->number_of_cells  = ncells;
        G->number_of_faces  = nfaces;
        G->number_of_nodes  = nnodes;

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

        G->cell_facetag   = malloc(nel * sizeof *G->cell_facetag);

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
            (G->cell_facetag     == NULL) ||
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


#define GRID_NMETA      6
#define GRID_NDIMS      0
#define GRID_NCELLS     1
#define GRID_NFACES     2
#define GRID_NNODES     3
#define GRID_NFACENODES 4
#define GRID_NCELLFACES 5


static void
input_error(FILE *fp, const char * const err)
{
    int save_errno = errno;

    if (ferror(fp)) {
        fprintf(stderr, "%s: %s\n", err, strerror(save_errno));
        clearerr(fp);
    }
    else if (feof(fp)) {
        fprintf(stderr, "%s: End-of-file\n", err);
    }

    errno = save_errno;
}


static struct UnstructuredGrid *
allocate_grid_from_file(FILE *fp, int *has_tag, int *has_indexmap)
{
    struct UnstructuredGrid *G;

    int           save_errno;
    unsigned long tmp;
    size_t        dimens[GRID_NMETA], i;

    save_errno = errno;

    i = 0;
    while ((i < GRID_NMETA) && (fscanf(fp, " %lu", &tmp) == 1)) {
        dimens[i] = tmp;

        i += 1;
    }

    if (i == GRID_NMETA) {
        if (fscanf(fp, "%d %d", has_tag, has_indexmap) == 2) {
            G = allocate_grid(dimens[GRID_NDIMS]     ,
                              dimens[GRID_NCELLS]    ,
                              dimens[GRID_NFACES]    ,
                              dimens[GRID_NFACENODES],
                              dimens[GRID_NCELLFACES],
                              dimens[GRID_NNODES]    );

            if (G != NULL) {
                if (! *has_tag) {
                    free(G->cell_facetag);
                    G->cell_facetag = NULL;
                }

                if (*has_indexmap) {
                    G->global_cell =
                        malloc(dimens[GRID_NCELLS] * sizeof *G->global_cell);

                    /* Allocation failure checked elsewhere. */
                }

                G->number_of_cells = (int) dimens[GRID_NCELLS];
                G->number_of_faces = (int) dimens[GRID_NFACES];
                G->number_of_nodes = (int) dimens[GRID_NNODES];
                G->dimensions      = (int) dimens[GRID_NDIMS];

                i = 0;
                while ((i < dimens[GRID_NDIMS]) &&
                       (fscanf(fp, "%d", & G->cartdims[ i ]) == 1)) {
                    i += 1;
                }

                if (i < dimens[GRID_NDIMS]) {
                    input_error(fp, "Unable to read Cartesian dimensions");

                    destroy_grid(G);
                    G = NULL;
                }
                else {
                    /* Account for dimens[GRID_DIMS] < 3 */
                    size_t n = (sizeof G->cartdims) / (sizeof G->cartdims[0]);
                    for (; i < n; i++) { G->cartdims[ i ] = 1; }
                }
            }
        }
        else {
            input_error(fp, "Unable to read grid predicates");

            G = NULL;
        }
    }
    else {
        input_error(fp, "Unable to read grid dimensions");

        G = NULL;
    }

    errno = save_errno;

    return G;
}


static int
read_grid_nodes(FILE *fp, struct UnstructuredGrid *G)
{
    int    save_errno;
    size_t i, n;

    save_errno = errno;

    n  = G->dimensions;
    n *= G->number_of_nodes;

    i = 0;
    while ((i < n) &&
           (fscanf(fp, " %lf", & G->node_coordinates[ i ]) == 1)) {
        i += 1;
    }

    if (i < n) {
        input_error(fp, "Unable to read node coordinates");
    }

    errno = save_errno;

    return i == n;
}


static int
read_grid_faces(FILE *fp, struct UnstructuredGrid *G)
{
    int    save_errno, ok;
    size_t nf, nfn, i;

    save_errno = errno;

    nf = G->number_of_faces;

    /* G->face_nodepos */
    i = 0;
    while ((i < nf + 1) &&
           (fscanf(fp, " %d", & G->face_nodepos[ i ]) == 1)) {
        i += 1;
    }
    ok = i == nf + 1;

    if (! ok) {
        input_error(fp, "Unable to read node indirection array");
    }
    else {
        /* G->face_nodes */
        nfn = G->face_nodepos[ nf ];

        i = 0;
        while ((i < nfn) && (fscanf(fp, " %d", & G->face_nodes[ i ]) == 1)) {
            i += 1;
        }

        ok = i == nfn;
        if (! ok) {
            input_error(fp, "Unable to read face-nodes");
        }
    }

    if (ok) {
        /* G->face_cells */
        i = 0;
        while ((i < 2 * nf) && (fscanf(fp, " %d", & G->face_cells[ i ]) == 1)) {
            i += 1;
        }

        ok = i == 2 * nf;
        if (! ok) {
            input_error(fp, "Unable to read neighbourship");
        }
    }

    if (ok) {
        /* G->face_areas */
        i = 0;
        while ((i < nf) && (fscanf(fp, " %lf", & G->face_areas[ i ]) == 1)) {
            i += 1;
        }

        ok = i == nf;
        if (! ok) {
            input_error(fp, "Unable to read face areas");
        }
    }

    if (ok) {
        /* G->face_centroids */
        size_t n;

        n  = G->dimensions;
        n *= nf;

        i = 0;
        while ((i < n) && (fscanf(fp, " %lf", & G->face_centroids[ i ]) == 1)) {
            i += 1;
        }

        ok = i == n;
        if (! ok) {
            input_error(fp, "Unable to read face centroids");
        }
    }

    if (ok) {
        /* G->face_normals */
        size_t n;

        n  = G->dimensions;
        n *= nf;

        i = 0;
        while ((i < n) && (fscanf(fp, " %lf", & G->face_normals[ i ]) == 1)) {
            i += 1;
        }

        ok = i == n;
        if (! ok) {
            input_error(fp, "Unable to read face normals");
        }
    }

    errno = save_errno;

    return ok;
}


static int
read_grid_cells(FILE *fp, int has_tag, int has_indexmap,
                struct UnstructuredGrid *G)
{
    int    save_errno, ok;
    size_t nc, ncf, i;

    save_errno = errno;

    nc = G->number_of_cells;

    /* G->cell_facepos */
    i = 0;
    while ((i < nc + 1) && (fscanf(fp, " %d", & G->cell_facepos[ i ]) == 1)) {
        i += 1;
    }
    ok = i == nc + 1;

    if (! ok) {
        input_error(fp, "Unable to read face indirection array");
    }
    else {
        /* G->cell_faces (and G->cell_facetag if applicable) */
        ncf = G->cell_facepos[ nc ];
        i   = 0;

        if (has_tag) {
            assert (G->cell_facetag != NULL);

            while ((i < ncf) &&
                   (fscanf(fp, " %d %d",
                           & G->cell_faces  [ i ],
                           & G->cell_facetag[ i ]) == 2)) {
                i += 1;
            }
        }
        else {
            while ((i < ncf) &&
                   (fscanf(fp, " %d", & G->cell_faces[ i ]) == 1)) {
                i += 1;
            }
        }

        ok = i == ncf;
        if (! ok) {
            input_error(fp, "Unable to read cell-faces");
        }
    }

    if (ok) {
        /* G->global_cell if applicable */
        if (has_indexmap) {
            i = 0;

            if (G->global_cell != NULL) {
                while ((i < nc) &&
                       (fscanf(fp, " %d", & G->global_cell[ i ]) == 1)) {
                    i += 1;
                }
            }
            else {
                int discard;

                while ((i < nc) && (fscanf(fp, " %d", & discard) == 1)) {
                    i += 1;
                }
            }
        }
        else {
            assert (G->global_cell == NULL);
            i = nc;
        }

        ok = i == nc;
        if (! ok) {
            input_error(fp, "Unable to read global cellmap");
        }
    }

    if (ok) {
        /* G->cell_volumes */
        i = 0;
        while ((i < nc) && (fscanf(fp, " %lf", & G->cell_volumes[ i ]) == 1)) {
            i += 1;
        }

        ok = i == nc;
        if (! ok) {
            input_error(fp, "Unable to read cell volumes");
        }
    }

    if (ok) {
        /* G->cell_centroids */
        size_t n;

        n  = G->dimensions;
        n *= nc;

        i = 0;
        while ((i < n) && (fscanf(fp, " %lf", & G->cell_centroids[ i ]) == 1)) {
            i += 1;
        }

        ok = i == n;
        if (! ok) {
            input_error(fp, "Unable to read cell centroids");
        }
    }

    errno = save_errno;

    return ok;
}


struct UnstructuredGrid *
read_grid(const char *fname)
{
    struct UnstructuredGrid *G;
    FILE                    *fp;

    int save_errno;
    int has_tag, has_indexmap, ok;

    save_errno = errno;

    fp = fopen(fname, "rt");
    if (fp != NULL) {
        G = allocate_grid_from_file(fp, & has_tag, & has_indexmap);

        ok = G != NULL;

        if (ok) { ok = read_grid_nodes(fp, G); }
        if (ok) { ok = read_grid_faces(fp, G); }
        if (ok) { ok = read_grid_cells(fp, has_tag, has_indexmap, G); }

        if (! ok) {
            destroy_grid(G);
            G = NULL;
        }

        fclose(fp);
    }
    else {
        G = NULL;
    }

    errno = save_errno;

    return G;
}
