/*
  Copyright 2010, 2011, 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_GRID_HEADER_INCLUDED
#define OPM_GRID_HEADER_INCLUDED

#include <stddef.h>

/**
 * \file
 *
 * Main OPM-Core grid data structure along with destructor and default
 * constructor.
 */

#ifdef __cplusplus
extern "C" {
#endif

/*
---- synopsis of grid.h ----

struct UnstructuredGrid
{
    int    dimensions;
    int    number_of_cells;
    int    number_of_faces;
    int    number_of_nodes;
    int    *face_nodes;
    int    *face_nodepos;
    int    *face_cells;
    int    *cell_faces;
    int    *cell_facepos;
    double *node_coordinates;
    double *face_centroids;
    double *face_areas;
    double *face_normals;
    double *cell_centroids;
    double *cell_volumes;
    int    *global_cell;
    int     cartdims[3];
    int    *cell_facetag;
};

void destroy_grid(struct UnstructuredGrid *g);

struct UnstructuredGrid *
allocate_grid(size_t ndims     ,
              size_t ncells    ,
              size_t nfaces    ,
              size_t nfacenodes,
              size_t ncellfaces,
              size_t nnodes    );

 ---- end of synopsis of grid.h ----
*/

/**
   Data structure for an unstructured grid, unstructured meaning that
   any cell may have an arbitrary number of adjacent cells. The struct
   contains both topological and geometrical data.

   The grid consists of a set of cells, which are assumed to partion
   the grid domain. A face is defined as the nonempty intersection of
   (the closure of) two grid cells (the grid is a cell complex).

   The topology information is limited to some adjacency relations
   between cells, faces and nodes only. The data structure does not
   contain any information pertaining to edges (except in 2d, where
   edges are the same as faces).

   The geometry information is limited to centroids, areas/volumes and
   normals.
 */
struct UnstructuredGrid
{
    /**
       The topological and geometrical dimensionality of the
       grid. Note that we do not support grids that are embedded in
       higher-dimensional spaces, such as 2d grids embedded in 3d.
       This number must be 2 or 3.
    */
    int    dimensions;

    /** The number of cells in the grid. */
    int    number_of_cells;
    /** The number of faces in the grid. */
    int    number_of_faces;
    /** The number of nodes in the grid. */
    int    number_of_nodes;

    /**
       Contains for each face, the indices of its adjacent nodes.
       The size of the array is equal to the sum over all faces of
       each face's number of adjacent nodes, which also is equal to
       face_nodepos[number_of_faces].
    */
    int    *face_nodes;
    /**
       For a face f, face_nodepos[f] contains the starting index
       for f's nodes in the face_nodes array.
       The size of the array is equal to (number_of_faces + 1).
    */
    int    *face_nodepos;
    /**
       For a face f, face_cells[2*f] and face_cells[2*f + 1] contain
       the cell indices of the cells adjacent to f. The number -1
       indicates the outer boundary.
       The order is significant, as it gives the orientation: if
       face_cells[2*f] == a and face_cells[2*f + 1] == b, f is
       oriented from a to b.  The inverse of this mapping is stored in
       cell_faces and cell_facepos.
       The size of the array is equal to (2*number_of_faces).
    */
    int    *face_cells;

    /**
       Contains for each cell, the indices of its adjacent faces.
       The size of the array is equal to the sum over all cells of
       each cell's number of adjacent faces, which also is equal to
       cell_facepos[number_of_cells].
    */
    int    *cell_faces;
    /**
       For a face c, cell_facepos[c] contains the starting index
       for c's faces in the cell_faces array.
       The size of the array is equal to (number_of_cells + 1).
    */
    int    *cell_facepos;

    /**
       Node coordinates, stored consecutively for each node. That is,
       for a node i, node_coordinates[dimensions*i + d] contains the
       d'th coordinate of node i.
       The size of the array is equal to (dimensions*number_of_nodes).
    */
    double *node_coordinates;

    /**
       Exact or approximate face centroids, stored consecutively for each face. That is,
       for a face f, face_centroids[dimensions*f + d] contains the
       d'th coordinate of f's centroid.
       The size of the array is equal to (dimensions*number_of_faces).
    */
    double *face_centroids;
    /**
       Exact or approximate face areas.
       The size of the array is equal to number_of_faces.
    */
    double *face_areas;
    /**
       Exact or approximate face normals, stored consecutively for
       each face. That is, for a face f, face_normals[dimensions*f + d]
       contains the d'th coordinate of f's normal.
       The size of the array is equal to (dimensions*number_of_faces).

       IMPORTANT: the normals are not normalized to have unit length!
       They are assumed to always have length equal to the
       corresponding face's area.
    */
    double *face_normals;

    /**
       Exact or approximate cell centroids, stored consecutively for each cell. That is,
       for a cell c, cell_centroids[dimensions*c + d] contains the
       d'th coordinate of c's centroid.
       The size of the array is equal to (dimensions*number_of_cells).
    */
    double *cell_centroids;
    /**
       Exact or approximate cell volumes.
       The size of the array is equal to number_of_cells.
    */
    double *cell_volumes;


    /**
       If non-null, this array contains the logical cartesian indices
       (in a lexicographic ordering) of each cell.
       In that case, the array size is equal to number_of_cells.

       This field is intended for grids that have a (possibly
       degenerate) logical cartesian structure, for example
       cornerpoint grids.
    */
    int    *global_cell;

    /**
       Contains the size of the logical cartesian structure (if any) of the grid.

       This field is intended for grids that have a (possibly
       degenerate) logical cartesian structure, for example
       cornerpoint grids.
    */
    int     cartdims[3];
    /**
       If non-null, this array contains a number for cell-face
       adjacency indicating the face's position with respect to the
       cell, in a logical cartesian sense. The tags are in [0, ..., 5]
       meaning [I-, I+, J-, J+, K-, K+], where I, J, K are the logical
       cartesian principal directions.
       The structure of this array is identical to cell_faces, and
       cell_facepos indices into cell_facetag as well.

       If non-null, the array size is equal to
       cell_facepos[number_of_cells].

       This field is intended for grids that have a (possibly
       degenerate) logical cartesian structure, for example
       cornerpoint grids.
    */
    int    *cell_facetag;
};

/**
   Destroy and deallocate an UnstructuredGrid and all its data.

   This function assumes that all arrays of the UnstructuredGrid (if
   non-null) have been individually allocated by malloc(). They will
   be deallocated with free().
 */
void destroy_grid(struct UnstructuredGrid *g);

/**
   Allocate and initialise an empty UnstructuredGrid.

   This is the moral equivalent of a C++ default constructor.

   \return Fully formed UnstructuredGrid with all fields zero or
   <code>NULL</code> as appropriate.  <code>NULL</code> in case of
   allocation failure.
 */
struct UnstructuredGrid *
create_grid_empty(void);

/**
   Allocate and initialise an UnstructuredGrid where pointers are set
   to location with correct size.

   \param[in] ndims      Number of physical dimensions.
   \param[in] ncells     Number of cells.
   \param[in] nfaces     Number of faces.
   \param[in] nfacenodes Size of mapping from faces to nodes.
   \param[in] ncellfaces Size of mapping from cells to faces.
                         (i.e., the number of `half-faces')
   \param[in] nnodes     Number of Nodes.

   \return Fully formed UnstructuredGrid with all fields except
   <code>global_cell</code> allocated, but not filled with meaningful
   values.  <code>NULL</code> in case of allocation failure.
 */
struct UnstructuredGrid *
allocate_grid(size_t ndims     ,
              size_t ncells    ,
              size_t nfaces    ,
              size_t nfacenodes,
              size_t ncellfaces,
              size_t nnodes    );


/**
 * Import a grid from a character representation stored in file.
 *
 * @param[in] fname File name.
 * @return Fully formed UnstructuredGrid with all fields allocated and filled.
 * Returns @c NULL in case of allocation failure.
 */
struct UnstructuredGrid *
read_grid(const char *fname);

#ifdef __cplusplus
}
#endif

#endif /* OPM_GRID_HEADER_INCLUDED */
