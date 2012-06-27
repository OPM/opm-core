/*===========================================================================
//
// File: preprocess.h
//
// Created: Fri Jun 19 08:43:04 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//==========================================================================*/

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_PREPROCESS_HEADER
#define OPENRS_PREPROCESS_HEADER


#ifdef __cplusplus
extern "C" {
#endif

    /* Input structure holding raw cornerpoint spec. */
    struct grdecl{
        int           dims[3];
        const double *coord;
        const double *zcorn;
        const int    *actnum;
        const double *mapaxes;  /* 6 Element rotation vector - can be NULL. */
    };

    /* Constant:     I     J     K    */
    enum face_tag { LEFT, BACK, TOP };

    /* Output structure holding grid topology */
    struct processed_grid{
        int m;                    /** Upper bound on "number_of_faces" */
        int n;                    /** Upper bound on "number_of_nodes" */
        int    dimensions[3];     /* Cartesian dimension                                */
        int    number_of_faces;
        int    *face_nodes;       /* Nodes numbers of each face sequentially.           */
        int    *face_ptr;         /* Start position for each face in face_nodes.        */
        int    *face_neighbors;   /* Global cell numbers.  2 ints per face sequentially */
        enum face_tag *face_tag;

        int    number_of_nodes;
        int    number_of_nodes_on_pillars; /** Total number of unique cell vertices that lie on pillars. */
        double *node_coordinates; /* 3 doubles per node, sequentially                   */

        int    number_of_cells;   /* number of active cells                             */
        int    *local_cell_index; /* Global to local map                                */
    };


    void process_grdecl     (const struct grdecl   *g,
                             double                tol,
                             struct processed_grid *out);
    void free_processed_grid(struct processed_grid *g);

#ifdef __cplusplus
}
#endif


#endif /* OPENRS_PREPROCESS_HEADER */


/* Local Variables:    */
/* c-basic-offset:4    */
/* End:                */
