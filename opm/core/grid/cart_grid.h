/*===========================================================================
//
// File: cart_grid.h
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

#ifndef OPM_CART_GRID_H_HEADER
#define OPM_CART_GRID_H_HEADER

/**
 * \file
 * Routines to construct fully formed grid structures from a simple Cartesian
 * (i.e., tensor product) description.
 *
 * The cells are lexicographically ordered with the @c i index cycling the most
 * rapidly, followed by the @c j index and then, in three space dimensions, the
 * @c k (`layer') index as the least rapidly cycling index.
 */

#ifdef __cplusplus
extern "C" {
#endif



/**
 * Form geometrically Cartesian grid in two space dimensions with equally
 * sized cells.
 *
 * @param[in] nx Number of cells in @c x direction.
 * @param[in] ny Number of cells in @c y direction.
 * @param[in] dx Length, in meters, of each cell's @c x extent.
 * @param[in] dy Length, in meters, of each cell's @c y extent.
 *
 * @return Fully formed grid structure containing valid geometric primitives.
 * Must be destroyed using function destroy_grid().
 */
struct UnstructuredGrid *
create_grid_cart2d(int nx, int ny, double dx, double dy);


/**
 * Form geometrically Cartesian grid in three space dimensions with unit-sized
 * cells.
 *
 * @param[in] nx Number of cells in @c x direction.
 * @param[in] ny Number of cells in @c y direction.
 * @param[in] nz Number of cells in @c z direction.
 *
 * @return Fully formed grid structure containing valid geometric primitives.
 * Must be destroyed using function destroy_grid().
 */
struct UnstructuredGrid *
create_grid_cart3d(int nx, int ny, int nz);


/**
 * Form geometrically Cartesian grid in three space dimensions with equally
 * sized cells.
 *
 * Each cell has physical size (volume) \f$\mathit{dx}\times \mathit{dy}\times
 * \mathit{dz}\f$.
 *
 * @param[in] nx Number of cells in @c x direction.
 * @param[in] ny Number of cells in @c y direction.
 * @param[in] nz Number of cells in @c z direction.
 *
 * @param[in] dx Length, in meters, of each cell's @c x extent.
 * @param[in] dy Length, in meters, of each cell's @c y extent.
 * @param[in] dz Length, in meters, of each cell's @c z extent.
 *
 * @return Fully formed grid structure containing valid geometric primitives.
 * Must be destroyed using function destroy_grid().
 */
struct UnstructuredGrid *
create_grid_hexa3d(int    nx, int    ny, int    nz,
                   double dx, double dy, double dz);


/**
 * Form tensor product (Cartesian) grid in two space dimensions.
 *
 * The size (volume) of cell \f$(i,j)\f$ is
 * \f[
 * v_{ij} = (x_{i+1} - x_i)\cdot (y_{j+1} - y_j)
 * \f]
 * Similar relations hold for the cell and interface centroids as well as the
 * interface areas and normal vectors. In other words, cell \f$(i,j)\f$ is the
 * convex hull bounded by the tensor product of nodes \f$x_i\f$, \f$x_{i+1}\f$,
 * \f$y_j\f$, and \f$y_{j+1}\f$.
 *
 * @param[in] nx Number of cells in @c x direction.
 * @param[in] ny Number of cells in @c y direction.
 *
 * @param[in] x  Position along @c x axis of each grid line with constant @c x
 *               coordinate.  Array of size <CODE>nx + 1</CODE>.
 * @param[in] y  Position along @c y axis of each grid line with constant @c y
 *               coordinate.  Array of size <CODE>ny + 1</CODE>.
 *
 * @return Fully formed grid structure containing valid geometric primitives.
 * Must be destroyed using function destroy_grid().
 */
struct UnstructuredGrid *
create_grid_tensor2d(int           nx, int           ny,
                     const double *x , const double *y );


/**
 * Form tensor product (i.e., topologically Cartesian) grid in three space
 * dimensions--possibly with a variable top-layer topography.
 *
 * If @c depthz is @c NULL, then geometric information such as volumes and
 * centroids is calculated from analytic expressions.  Otherwise, these values
 * are computed using function compute_geometry().
 *
 * @param[in] nx Number of cells in @c x direction.
 * @param[in] ny Number of cells in @c y direction.
 * @param[in] nz Number of cells in @c z direction.
 *
 * @param[in] x  Position along @c x axis of each grid line with constant @c x
 *               coordinate.  Array of size <CODE>nx + 1</CODE>.
 * @param[in] y  Position along @c y axis of each grid line with constant @c y
 *               coordinate.  Array of size <CODE>ny + 1</CODE>.
 * @param[in] z  Distance (depth) from top-layer measured along the @c z axis of
 *               each grid line with constant @c z coordinate.  Array of size
 *               <CODE>nz + 1</CODE>.
 *
 * @param[in] depthz
 *               Top-layer topography specification.  If @c NULL, interpreted as
 *               horizontal top-layer at <CODE>z=0</CODE>.  Otherwise, must be
 *               an array of size <CODE>(nx + 1) * (ny + 1)</CODE>, ordered
 *               lexicographically, that defines the depth of each top-layer
 *               pillar vertex.
 *
 * @return Fully formed grid structure containing valid geometric primitives.
 * Must be destroyed using function destroy_grid().
 */
struct UnstructuredGrid *
create_grid_tensor3d(int           nx    ,
                     int           ny    ,
                     int           nz    ,
                     const double *x     ,
                     const double *y     ,
                     const double *z     ,
                     const double *depthz);
#ifdef __cplusplus
}
#endif
#endif  /* OPM_CART_GRID_H_HEADER */
