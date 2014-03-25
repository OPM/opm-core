/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2014 Dr. Blatt - HPC-Simulation-Software & Services

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

#ifndef OPM_TRANSTPFA_HEADER_INCLUDED
#define OPM_TRANSTPFA_HEADER_INCLUDED

/**
 * \file
 * Routines to assist in the calculation of two-point transmissibilities.
 */

/**
 * Calculate static, one-sided transmissibilities for use in the two-point flux
 * approximation method.
 *
 * The one-sided transmissibilities are defined by the formula
 * \f[
 * t_i = \frac{\vec{n}_f \mathsf{K}_c \vec{c}_c}{\lVert \vec{c}_c \rVert^2}
 * \f]
 * in which @c i is the half-face index corresponding to the cell-face index
 * pair <CODE>(c,f)</CODE> and \f$\vec{c}_{cf} = \Bar{x}_f - \Bar{x}_c\f$ is the
 * centroid difference vector.
 *
 * @param[in]  G       Grid.
 * @param[in]  perm    Permeability.  One symmetric, positive definite tensor
 *                     per grid cell.
 * @param[out] htrans  One-sided transmissibilities.  Array of size at least
 *                     <CODE>G->cell_facepos[ G->number_of_cells  ]</CODE>.
 */
template<class Grid>
void
tpfa_htrans_compute(const Grid   *G      ,
                    const double *perm  ,
                    double       *htrans);

/**
 * Compute two-point transmissibilities from one-sided transmissibilities.
 *
 * The two-point transmissibilities are given by the simple, near-harmonic
 * average (save a factor of two)
 * \f[
 * \mathsf{T}_f = \big(\frac{1}{t_1} + \frac{1}{t_2}\big)^{-1}
 *              = \frac{t_1t_2}{t_1 + t_2}
 * \f]
 * in which \f$t_1\f$ and \f$t_2\f$ are the one-sided transmissibilities that
 * connect the neighbouring cells of face @c f.
 *
 * @param[in]  G       Grid.
 * @param[in]  htrans  One-sided transmissibilities as defined by function
 *                     tpfa_htrans_compute().
 * @param[out] trans   Interface, two-point transmissibilities.  Array of size
 *                     at least <CODE>numFaces(G)</CODE>.
 */
template<class Grid>
void
tpfa_trans_compute(const Grid   *G      ,
                   const double *htrans,
                   double       *trans );

/**
 * Calculate effective two-point transmissibilities from one-sided, total
 * mobility weighted, transmissibilities.
 *
 * Specifically, compute the following product
 * \f[
 * \mathsf{T}_f = \big(\frac{1}{\lambda_1t_1} + \frac{1}{\lambda_2t_2}\big)^{-1}
 *              = \lambda_1\lambda_2 \frac{t_1t_2}{t_1 + t_2}
 * \f]
 * in which \f$t_1\f$ and \f$t_2\f$ are the one-sided, static transmissibility
 * values connecting the cells of face @c f and \f$\lambda_1\f$ and
 * \f$\lambda_2\f$ denote the total mobilities of the respective cells.
 *
 * @param[in]  G      Grid.
 * @param[in]  totmob Total mobilities. One positive scalar value for each cell.
 * @param[in]  htrans One-sided transmissibilities as defined by function
 *                    tpfa_htrans_compute().
 * @param[out] trans  Effective, two-point transmissibilities.  Array of size at
 *                    least <CODE>numFaces(G)</CODE>.
 */
template<class Grid>
void
tpfa_eff_trans_compute(const Grid   *G      ,
                       const double *totmob,
                       const double *htrans,
                       double       *trans );

#include "TransTpfa_impl.hpp"
#endif  /* OPM_TRANS_TPFA_HEADER_INCLUDED */
