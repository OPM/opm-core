/*===========================================================================
//
// File: preprocess.h
//
// Created: Fri Jun 19 08:43:04 2009
//
// Author: Jostein R. Natvig <Jostein.R.Natvig@sintef.no>
//
// $Date: 2010-08-27 19:12:16 +0200 (Fri, 27 Aug 2010) $
//
// $Revision: 930 $
//
//==========================================================================*/

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

#ifndef OPM_CORNERPOINT_GRID_HEADER_INCLUDED
#define OPM_CORNERPOINT_GRID_HEADER_INCLUDED

/**
 * \file
 * Routines to form a complete UnstructuredGrid from a corner-point
 * specification.
 */

#include <opm/core/grid.h>
#include <opm/core/grid/cpgpreprocess/preprocess.h>

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * Construct grid representation from corner-point specification of a
     * particular geological model.
     *
     * Pinched cells will be removed irrespective of any explicit "active" map
     * in the geological model input specification.  Geometric primitives such
     * as cell barycenters (i.e., centroids), volumes and interface areas are
     * computed internally using function compute_geometry().  The caller does
     * not need to compute this information separately.
     *
     * @param[in] in  Corner-point specification.  If "actnum" is NULL, then the
     *                specification is interpreted as if all cells are initially
     *                active.
     *
     * @param[in] tol Absolute tolerance of node-coincidence.
     *
     * @return Fully formed grid data structure that manages the grid defined by
     * the input corner-point specification. Must be destroyed using function
     * destroy_grid().
     */
    struct UnstructuredGrid *
    create_grid_cornerpoint(const struct grdecl *in, double tol);


    /**
     * Compute derived geometric primitives in a grid.
     *
     * This function computes values for each of the following quantities
     * - Quantities pertaining to interfaces (connections, faces)
     *   -# Barycenters (centroids), <CODE>g->dimensions</CODE> scalars per face
     *      stored sequentially in </CODE>g->face_centroids</CODE>.
     *   -# Areas, one scalar per face stored sequentially in
     *      <CODE>g->face_areas</CODE>.
     *   -# Normals, <CODE>g->dimensions</CODE> scalars per face stored
     *      sequentially in <CODE>g->face_normals</CODE>.  The Euclidian norm of
     *      each normal is equal to the corresponding face's area.
     *
     * - Quantities pertaining to cells (volumes)
     *   -# Barycenters (centroids), <CODE>g->dimensions</CODE> scalars per cell
     *      stored sequentially in <CODE>g->cell_centroids</CODE>.
     *   -# Volumes, one scalar per cell stored sequentially in
     *      <CODE>g->cell_volumes</CODE>.
     *
     * These fields must be allocated prior to calling compute_geometry().
     *
     * @param[in,out] g Grid structure.
     */
    void compute_geometry(struct UnstructuredGrid *g);
    
#ifdef __cplusplus
}
#endif

#endif /* OPM_CORNERPOINT_GRID_HEADER_INCLUDED */
