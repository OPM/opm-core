/*===========================================================================
//
// File: cart_grid.c
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
//==========================================================================*/


/*
  Copyright 2011 SINTEF ICT, Applied Mathematics.
  Copyright 2011 Statoil ASA.

  This file is part of the Open Porous Media Project (OPM).

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

#include <stdlib.h>
#include <stdio.h>
#include "cart_grid.h"

static struct UnstructuredGrid *allocate_cart_grid_3d(int nx, int ny, int nz);
static void fill_cart_topology_3d(struct UnstructuredGrid *G, int nx, int ny, int nz);
static void fill_cart_geometry_3d(struct UnstructuredGrid *G, int nx, int ny, int nz);

struct UnstructuredGrid*
create_cart_grid_3d(int nx, int ny, int nz)
{
    struct UnstructuredGrid *G;

    G = allocate_cart_grid_3d(nx, ny, nz);
    if (G == NULL)
    {
        fprintf(stderr, "Cannot allocate space for grid.\n");
        exit(EXIT_FAILURE);
    }

    fill_cart_topology_3d(G, nx, ny, nz);
    fill_cart_geometry_3d(G, nx, ny, nz);

    return G;
}

/* --------------------------------------------------------------------- */

static struct UnstructuredGrid *allocate_cart_grid_2d(int nx, int ny);
static void fill_cart_topology_2d(struct UnstructuredGrid *G, int nx, int nz);
static void fill_cart_geometry_2d(struct UnstructuredGrid *G, int nx, int nz);

struct UnstructuredGrid*
create_cart_grid_2d(int nx, int ny)
{
    struct UnstructuredGrid *G;

    G = allocate_cart_grid_2d(nx, ny);
    if (G == NULL)
    {
        fprintf(stderr, "Cannot allocate space for grid.\n");
        exit(EXIT_FAILURE);
    }

    fill_cart_topology_2d(G, nx, ny);
    fill_cart_geometry_2d(G, nx, ny);

    return G;
}

/* --------------------------------------------------------------------- */

struct UnstructuredGrid *
create_tensor_grid_2d(int nx, int ny, double x[], double y[])
{
    int    i,j;

    double *coord;
    struct UnstructuredGrid *G;

    G = allocate_cart_grid_2d(nx, ny);

    if (G == NULL)
    {
        fprintf(stderr, "Cannot allocate space for grid.\n");
        exit(EXIT_FAILURE);
    }
    
    fill_cart_topology_2d(G, nx, ny);

    coord = G->node_coordinates;
    for (j=0; j<ny+1; ++j) {
        for (i=0; i<nx+1; ++i) {
            *coord++ = x[i];
            *coord++ = y[j];
        }
    }
    
    return G;
}

/* --------------------------------------------------------------------- */

struct UnstructuredGrid *
create_tensor_grid_3d(int nx, int ny, int nz, double x[], double y[], double z[])
{
    int    i,j,k;

    double *coord;
    struct UnstructuredGrid *G;

    G = allocate_cart_grid_3d(nx, ny, nz);

    if (G == NULL)
    {
        fprintf(stderr, "Cannot allocate space for grid.\n");
        exit(EXIT_FAILURE);
    }
    
    fill_cart_topology_3d(G, nx, ny, nz);

    coord = G->node_coordinates;
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx+1; ++i) {
                *coord++ = x[i];
                *coord++ = y[j];
                *coord++ = z[k];
            }
        }
    }

    return G;
}

/* --------------------------------------------------------------------- */

void
destroy_cart_grid(struct UnstructuredGrid *G)
{
    if (G != NULL) {
        free(G->node_coordinates);

        free(G->face_nodes);
        free(G->face_nodepos);
        free(G->face_cells);
        free(G->face_centroids);
        free(G->face_normals);
        free(G->face_areas);

        free(G->cell_faces);
        free(G->cell_facepos);
        free(G->cell_centroids);
        free(G->cell_volumes);
    }

    free(G);
}

/* --------------------------------------------------------------------- */
/* Static functions follow:                                              */
/* --------------------------------------------------------------------- */

static struct UnstructuredGrid*
allocate_cart_grid_3d(int nx, int ny, int nz)
{
    struct UnstructuredGrid *G;
    int Nx, Ny, Nz;
    int nxf, nyf, nzf;
    G = malloc(1 * sizeof *G);
    if (G != NULL)
    {

        G->dimensions = 3;

        Nx  = nx+1;
        Ny  = ny+1;
        Nz  = nz+1;

        nxf = Nx*ny*nz;
        nyf = nx*Ny*nz;
        nzf = nx*ny*Nz;

        G->number_of_cells  = nx*ny*nz;
        G->number_of_faces  = nxf+nyf+nzf;
        G->number_of_nodes  = Nx*Ny*Nz;

        G->node_coordinates = malloc(G->number_of_nodes * 3 * sizeof *(G->node_coordinates));

        G->face_nodes       = malloc(G->number_of_faces * 4 * sizeof *(G->face_nodes));
        G->face_nodepos     = malloc((G->number_of_faces+1) * sizeof *(G->face_nodepos));
        G->face_cells       = malloc(G->number_of_faces * 2 * sizeof *(G->face_cells));
        G->face_centroids   = malloc(G->number_of_faces * 3 * sizeof *(G->face_centroids));
        G->face_normals     = malloc(G->number_of_faces * 3 * sizeof *(G->face_normals));
        G->face_areas       = malloc(G->number_of_faces * 1 * sizeof *(G->face_areas));

        G->cell_faces       = malloc(G->number_of_cells * 6 * sizeof *(G->cell_faces));
        G->cell_facepos     = malloc((G->number_of_cells+1) * sizeof *(G->cell_facepos));
        G->cell_centroids   = malloc(G->number_of_cells * 3 * sizeof *(G->cell_centroids));
        G->cell_volumes     = malloc(G->number_of_cells * 1 * sizeof *(G->cell_volumes));

        if ((G->face_nodes       == NULL ) ||
            (G->face_nodepos     == NULL ) || 
            (G->face_cells       == NULL ) || 
            (G->face_centroids   == NULL ) || 
            (G->face_normals     == NULL ) || 
            (G->face_areas       == NULL ) ||                       
            (G->cell_faces       == NULL ) || 
            (G->cell_facepos     == NULL ) || 
            (G->cell_centroids   == NULL ) || 
            (G->cell_volumes     == NULL )  )
        {
            destroy_cart_grid(G);
            G = NULL;
        }
    }

    return G;
}

/* --------------------------------------------------------------------- */

static void
fill_cart_topology_3d(struct UnstructuredGrid *G, int nx, int ny, int nz)
{
    int Nx, Ny, Nz;
    int nxf, nyf, nzf;
    int i,j,k;

    int *cfaces, *cfacepos, *fnodes, *fnodepos, *fcells;

    Nx  = nx+1;
    Ny  = ny+1;
    Nz  = nz+1;
    
    nxf = Nx*ny*nz;
    nyf = nx*Ny*nz;
    nzf = nx*ny*Nz;

    cfaces     = G->cell_faces;
    cfacepos   = G->cell_facepos;
    cfacepos[0] = 0;
    for (k=0; k<nz; ++k)  {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *cfaces++ = i+  Nx*(j+  ny* k   );
                *cfaces++ = i+1+Nx*(j+  ny* k   );
                *cfaces++ = i+  nx*(j+  Ny* k   )  +nxf;
                *cfaces++ = i+  nx*(j+1+Ny* k   )  +nxf;
                *cfaces++ = i+  nx*(j+  ny* k   )  +nxf+nyf;
                *cfaces++ = i+  nx*(j+  ny*(k+1))  +nxf+nyf;

                cfacepos[1] = cfacepos[0]+6;
                ++cfacepos;
            }
        }
    }


    fnodes     = G->face_nodes;
    fnodepos   = G->face_nodepos;
    fcells     = G->face_cells;

    /* Faces with x-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx+1; ++i) {
                *fnodes++ = i+Nx*(j   + Ny * k   );
                *fnodes++ = i+Nx*(j+1 + Ny * k   );
                *fnodes++ = i+Nx*(j+1 + Ny *(k+1));
                *fnodes++ = i+Nx*(j   + Ny *(k+1));
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (i==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (i == nx) {
                    *fcells++ =  i-1+nx*(j+ny*k);
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i-1 + nx*(j+ny*k);
                    *fcells++ =  i   + nx*(j+ny*k);
                }
            }
        }
    }
    /* Faces with y-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx; ++i) {
                *fnodes++ = i+    Nx*(j + Ny * k   );
                *fnodes++ = i   + Nx*(j + Ny *(k+1));
                *fnodes++ = i+1 + Nx*(j + Ny *(k+1));
                *fnodes++ = i+1 + Nx*(j + Ny * k   );
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (j==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (j == ny) {
                    *fcells++ =  i+nx*(j-1+ny*k);
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i+nx*(j-1+ny*k);
                    *fcells++ =  i+nx*(j+ny*k);
                }
            }
        }
    }
    /* Faces with z-normal */
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *fnodes++ = i+    Nx*(j   + Ny * k);
                *fnodes++ = i+1 + Nx*(j   + Ny * k);
                *fnodes++ = i+1 + Nx*(j+1 + Ny * k);
                *fnodes++ = i+    Nx*(j+1 + Ny * k);
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (k==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (k == nz) {
                    *fcells++ =  i+nx*(j+ny*(k-1));
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i+nx*(j+ny*(k-1));
                    *fcells++ =  i+nx*(j+ny*k);
                }
            }
        }
    }    
}

/* --------------------------------------------------------------------- */

static void
fill_cart_geometry_3d(struct UnstructuredGrid *G, int nx, int ny, int nz)
{
    int Nx, Ny, Nz;
    int nxf, nyf, nzf;
    int i,j,k;

    double *coord, *ccentroids, *cvolumes;
    double *fnormals, *fcentroids, *fareas;

    Nx  = nx+1;
    Ny  = ny+1;
    Nz  = nz+1;

    nxf = Nx*ny*nz;
    nyf = nx*Ny*nz;
    nzf = nx*ny*Nz;


    ccentroids = G->cell_centroids;
    cvolumes   = G->cell_volumes;
    for (k=0; k<nz; ++k)  {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *ccentroids++ = i+0.5;
                *ccentroids++ = j+0.5;
                *ccentroids++ = k+0.5;

                *cvolumes++ = 1;
            }
        }
    }


    fnormals   = G->face_normals;
    fcentroids = G->face_centroids;
    fareas     = G->face_areas;

    /* Faces with x-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx+1; ++i) {
                *fnormals++ = 1;
                *fnormals++ = 0;
                *fnormals++ = 0;

                *fcentroids++ = i;
                *fcentroids++ = j+0.5;
                *fcentroids++ = k+0.5;

                *fareas++ = 1;
            }
        }
    }
    /* Faces with y-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx; ++i) {
                *fnormals++ = 0;
                *fnormals++ = 1;
                *fnormals++ = 0;

                *fcentroids++ = i+0.5;
                *fcentroids++ = j;
                *fcentroids++ = k+0.5;

                *fareas++ = 1;
            }
        }
    }
    /* Faces with z-normal */
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *fnormals++ = 0;
                *fnormals++ = 0;
                *fnormals++ = 1;

                *fcentroids++ = i+0.5;
                *fcentroids++ = j+0.5;
                *fcentroids++ = k;

                *fareas++ = 1;
            }
        }
    }

    coord = G->node_coordinates;
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx+1; ++i) {
                *coord++ = i;
                *coord++ = j;
                *coord++ = k;
            }
        }
    }
}

/* --------------------------------------------------------------------- */

struct UnstructuredGrid*
create_cart_grid(int nx, int ny, int nz)
{
    int    i,j,k;
    int    nxf, nyf, nzf;
    int    Nx, Ny, Nz;

    struct UnstructuredGrid *G;
    double *coord, *ccentroids, *cvolumes;
    double *fnormals, *fcentroids, *fareas;

    int    *fnodes, *fnodepos, *fcells, *cfaces, *cfacepos;

    G = malloc(1 * sizeof *G);
    if (G == NULL)
    {
        fprintf(stderr, "Cannot allocate space for grid.\n");
        exit(EXIT_FAILURE);
    }

    G->dimensions = 3;

    Nx  = nx+1;
    Ny  = ny+1;
    Nz  = nz+1;

    nxf = Nx*ny*nz;
    nyf = nx*Ny*nz;
    nzf = nx*ny*Nz;

    G->number_of_cells  = nx*ny*nz;
    G->number_of_faces  = nxf+nyf+nzf;
    G->number_of_nodes  = Nx*Ny*Nz;

    G->node_coordinates = malloc(G->number_of_nodes * 3 * sizeof *(G->node_coordinates));

    G->face_nodes       = malloc(G->number_of_faces * 4 * sizeof *(G->face_nodes));
    G->face_nodepos     = malloc((G->number_of_faces+1) * sizeof *(G->face_nodepos));
    G->face_cells       = malloc(G->number_of_faces * 2 * sizeof *(G->face_cells));
    G->face_centroids   = malloc(G->number_of_faces * 3 * sizeof *(G->face_centroids));
    G->face_normals     = malloc(G->number_of_faces * 3 * sizeof *(G->face_normals));
    G->face_areas       = malloc(G->number_of_faces * 1 * sizeof *(G->face_areas));

    G->cell_faces       = malloc(G->number_of_cells * 6 * sizeof *(G->cell_faces));
    G->cell_facepos     = malloc((G->number_of_cells+1) * sizeof *(G->cell_facepos));
    G->cell_centroids   = malloc(G->number_of_cells * 3 * sizeof *(G->cell_centroids));
    G->cell_volumes     = malloc(G->number_of_cells * 1 * sizeof *(G->cell_volumes));

    if ((G->face_nodes       == NULL ) ||
        (G->face_nodepos     == NULL ) || 
        (G->face_cells       == NULL ) || 
        (G->face_centroids   == NULL ) || 
        (G->face_normals     == NULL ) || 
        (G->face_areas       == NULL ) ||                       
        (G->cell_faces       == NULL ) || 
        (G->cell_facepos     == NULL ) || 
        (G->cell_centroids   == NULL ) || 
        (G->cell_volumes     == NULL )  )
    {
        fprintf(stderr, "Cannot allocate space for grid.\n");
        exit(EXIT_FAILURE);
    }

    cfaces     = G->cell_faces;
    cfacepos   = G->cell_facepos;
    ccentroids = G->cell_centroids;
    cvolumes   = G->cell_volumes;
    cfacepos[0] = 0;
    for (k=0; k<nz; ++k)  {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *cfaces++ = i+  Nx*(j+  ny* k   );
                *cfaces++ = i+1+Nx*(j+  ny* k   );
                *cfaces++ = i+  nx*(j+  Ny* k   )  +nxf;
                *cfaces++ = i+  nx*(j+1+Ny* k   )  +nxf;
                *cfaces++ = i+  nx*(j+  ny* k   )  +nxf+nyf;
                *cfaces++ = i+  nx*(j+  ny*(k+1))  +nxf+nyf;

                cfacepos[1] = cfacepos[0]+6;
                ++cfacepos;

                *ccentroids++ = i+0.5;
                *ccentroids++ = j+0.5;
                *ccentroids++ = k+0.5;

                *cvolumes++ = 1;
            }
        }
    }


    fnodes     = G->face_nodes;
    fnodepos   = G->face_nodepos;
    fcells     = G->face_cells;
    fnormals   = G->face_normals;
    fcentroids = G->face_centroids;
    fareas     = G->face_areas;

    /* Faces with x-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx+1; ++i) {
                *fnodes++ = i+Nx*(j   + Ny * k   );
                *fnodes++ = i+Nx*(j+1 + Ny * k   );
                *fnodes++ = i+Nx*(j+1 + Ny *(k+1));
                *fnodes++ = i+Nx*(j   + Ny *(k+1));
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (i==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (i == nx) {
                    *fcells++ =  i-1+nx*(j+ny*k);
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i-1 + nx*(j+ny*k);
                    *fcells++ =  i   + nx*(j+ny*k);
                }

                *fnormals++ = 1;
                *fnormals++ = 0;
                *fnormals++ = 0;

                *fcentroids++ = i;
                *fcentroids++ = j+0.5;
                *fcentroids++ = k+0.5;

                *fareas++ = 1;
            }
        }
    }
    /* Faces with y-normal */
    for (k=0; k<nz; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx; ++i) {
                *fnodes++ = i+    Nx*(j + Ny * k   );
                *fnodes++ = i   + Nx*(j + Ny *(k+1));
                *fnodes++ = i+1 + Nx*(j + Ny *(k+1));
                *fnodes++ = i+1 + Nx*(j + Ny * k   );
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (j==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (j == ny) {
                    *fcells++ =  i+nx*(j-1+ny*k);
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i+nx*(j-1+ny*k);
                    *fcells++ =  i+nx*(j+ny*k);
                }

                *fnormals++ = 0;
                *fnormals++ = 1;
                *fnormals++ = 0;

                *fcentroids++ = i+0.5;
                *fcentroids++ = j;
                *fcentroids++ = k+0.5;

                *fareas++ = 1;
            }
        }
    }
    /* Faces with z-normal */
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny; ++j) {
            for (i=0; i<nx; ++i) {
                *fnodes++ = i+    Nx*(j   + Ny * k);
                *fnodes++ = i+1 + Nx*(j   + Ny * k);
                *fnodes++ = i+1 + Nx*(j+1 + Ny * k);
                *fnodes++ = i+    Nx*(j+1 + Ny * k);
                fnodepos[1] = fnodepos[0] + 4;
                ++fnodepos;
                if (k==0) {
                    *fcells++ = -1;
                    *fcells++ =  i+nx*(j+ny*k);
                }
                else if (k == nz) {
                    *fcells++ =  i+nx*(j+ny*(k-1));
                    *fcells++ = -1;
                }
                else {
                    *fcells++ =  i+nx*(j+ny*(k-1));
                    *fcells++ =  i+nx*(j+ny*k);
                }

                *fnormals++ = 0;
                *fnormals++ = 0;
                *fnormals++ = 1;

                *fcentroids++ = i+0.5;
                *fcentroids++ = j+0.5;
                *fcentroids++ = k;

                *fareas++ = 1;
            }
        }
    }

    coord = G->node_coordinates;
    for (k=0; k<nz+1; ++k) {
        for (j=0; j<ny+1; ++j) {
            for (i=0; i<nx+1; ++i) {
                *coord++ = i;
                *coord++ = j;
                *coord++ = k;
            }
        }
    }

    return G;
}

/* --------------------------------------------------------------------- */

static struct UnstructuredGrid*
allocate_cart_grid_2d(int nx, int ny)
{
    int    nxf, nyf;
    int    Nx, Ny;

    struct UnstructuredGrid *G;

    G = malloc(1 * sizeof *G);
    if (G != NULL)
    {
        G->dimensions = 2;

        Nx  = nx+1;
        Ny  = ny+1;

        nxf = Nx*ny;
        nyf = nx*Ny;

        G->number_of_cells  = nx*ny;
        G->number_of_faces  = nxf+nyf;
        G->number_of_nodes  = Nx*Ny;

        G->node_coordinates = malloc(G->number_of_nodes * 2 * sizeof *(G->node_coordinates));

        G->face_nodes       = malloc(G->number_of_faces * 2 * sizeof *(G->face_nodes));
        G->face_nodepos     = malloc((G->number_of_faces+1) * sizeof *(G->face_nodepos));
        G->face_cells       = malloc(G->number_of_faces * 2 * sizeof *(G->face_cells));
        G->face_centroids   = malloc(G->number_of_faces * 2 * sizeof *(G->face_centroids));
        G->face_normals     = malloc(G->number_of_faces * 2 * sizeof *(G->face_normals));
        G->face_areas       = malloc(G->number_of_faces * 1 * sizeof *(G->face_areas));

        G->cell_faces       = malloc(G->number_of_cells * 4 * sizeof *(G->cell_faces));
        G->cell_facepos     = malloc((G->number_of_cells+1) * sizeof *(G->cell_facepos));
        G->cell_centroids   = malloc(G->number_of_cells * 2 * sizeof *(G->cell_centroids));
        G->cell_volumes     = malloc(G->number_of_cells * 1 * sizeof *(G->cell_volumes));

        if ((G->face_nodes       == NULL ) ||
            (G->face_nodepos     == NULL ) || 
            (G->face_cells       == NULL ) || 
            (G->face_centroids   == NULL ) || 
            (G->face_normals     == NULL ) || 
            (G->face_areas       == NULL ) ||                       
            (G->cell_faces       == NULL ) || 
            (G->cell_facepos     == NULL ) || 
            (G->cell_centroids   == NULL ) || 
            (G->cell_volumes     == NULL )  )
        {
            destroy_cart_grid(G);
            G = NULL;
        }
    }

    return G;
}

/* --------------------------------------------------------------------- */

static void
fill_cart_topology_2d(struct UnstructuredGrid *G, int nx, int ny)
{
    int    i,j;
    int    nxf;
    int    Nx;

    int    *fnodes, *fnodepos, *fcells, *cfaces, *cfacepos;

    cfaces     = G->cell_faces;
    cfacepos   = G->cell_facepos;

    Nx  = nx + 1;
    nxf = Nx * ny;

    cfacepos[0] = 0;
    for (j=0; j<ny; ++j) {
        for (i=0; i<nx; ++i) {
            *cfaces++ = i+  Nx*j;
            *cfaces++ = i+  nx*j  +nxf;
            *cfaces++ = i+1+Nx*j;
            *cfaces++ = i+  nx*(j+1)+nxf;
            
            cfacepos[1] = cfacepos[0]+4;
            ++cfacepos;
        }
    }



    fnodes     = G->face_nodes;
    fnodepos   = G->face_nodepos;
    fcells     = G->face_cells;
    
    /* Faces with x-normal */
    for (j=0; j<ny; ++j) {
        for (i=0; i<nx+1; ++i) {
            *fnodes++ = i+Nx*j;
            *fnodes++ = i+Nx*j+1;
            fnodepos[1] = fnodepos[0] + 2;
            ++fnodepos;
            if (i==0) {
                *fcells++ = -1;
                *fcells++ =  i+nx*j;
            }
            else if (i == nx) {
                *fcells++ =  i-1+nx*j;
                *fcells++ = -1;
            }
            else {
                *fcells++ =  i-1 + nx*j;
                *fcells++ =  i   + nx*j;
            }
        }
    }

    /* Faces with y-normal */
    for (j=0; j<ny+1; ++j) {
        for (i=0; i<nx; ++i) {
            *fnodes++ = i+    Nx*j;
            *fnodes++ = i+1 + Nx*j;
            fnodepos[1] = fnodepos[0] + 2;
            ++fnodepos;
            if (j==0) {
                *fcells++ = -1;
                *fcells++ =  i+nx*j;
            }
            else if (j == ny) {
                *fcells++ =  i+nx*j;
                *fcells++ = -1;
            }
            else {
                *fcells++ =  i+nx*(j-1);
                *fcells++ =  i+nx*j;
            }
        }
    }
}

/* --------------------------------------------------------------------- */

static void
fill_cart_geometry_2d(struct UnstructuredGrid *G, int nx, int ny)
{
    int    i,j;
    int    nxf, nyf;
    int    Nx, Ny;

    double *coord, *ccentroids, *cvolumes;
    double *fnormals, *fcentroids, *fareas;

    Nx  = nx+1;
    Ny  = ny+1;

    nxf = Nx*ny;
    nyf = nx*Ny;


    ccentroids = G->cell_centroids;
    cvolumes   = G->cell_volumes;

    for (j=0; j<ny; ++j) {
        for (i=0; i<nx; ++i) {
            *ccentroids++ = i+0.5;
            *ccentroids++ = j+0.5;
            
            *cvolumes++ = 1;
        }
    }



    fnormals   = G->face_normals;
    fcentroids = G->face_centroids;
    fareas     = G->face_areas;
    
    /* Faces with x-normal */
    for (j=0; j<ny; ++j) {
        for (i=0; i<nx+1; ++i) {
            *fnormals++ = 1;
            *fnormals++ = 0;
            
            *fcentroids++ = i;
            *fcentroids++ = j+0.5;
            
            *fareas++ = 1;
        }
    }

    /* Faces with y-normal */
    for (j=0; j<ny+1; ++j) {
        for (i=0; i<nx; ++i) {
            *fnormals++ = 0;
            *fnormals++ = 1;
            
            *fcentroids++ = i+0.5;
            *fcentroids++ = j;
            
            *fareas++ = 1;
        }
    }
    
    coord = G->node_coordinates;
    for (j=0; j<ny+1; ++j) {
        for (i=0; i<nx+1; ++i) {
            *coord++ = i;
            *coord++ = j;
        }
    }
}





#if UNIT_TEST
int main()
{
    
    grid_t *g = create_cart_grid_2d(2,2);
    int i;
    int k;
    for(i=0; i<g->number_of_cells; ++i)
    {
        fprintf(stderr, "%d: ", i);
        for (k=g->cell_facepos[i]; k<g->cell_facepos[i+1]; ++k)
        {
            fprintf(stderr, "%d ", g->cell_faces[k]); 
        }
        fprintf(stderr, "\n");
    }
    destroy_cart_grid(g);
    return 0;
}
#endif
