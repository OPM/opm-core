/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_FSH_HEADER_INCLUDED
#define OPM_FSH_HEADER_INCLUDED

/**
 * \file
 * Routines and data structures to support the construction and
 * formation of hybridized pressure solvers based on Schur
 * complement reductions.
 *
 * Pressure solvers based on this strategy will be structured
 * according to the following scheme
 * -# Construct @c FSH data object suitable for the particular
 *    problem using either of the functions cfsh_construct() or
 *    ifsh_construct() for compressible or incompressible flow,
 *    respectively
 * -# Compute static discretisation quantities, for instance
 *    using functions mim_ip_simple_all() and mim_ip_compute_gpress()
 * -# For each time step or non-linear iteration:
 *    -# Compute dynamic discretisation quantities incorporating
 *       effects of multiple phases and non-linear feedback.
 *    -# Assemble system of simultaneous linear equations using
 *       functions cfsh_assemble() or ifsh_assemble()
 *    -# Solve the resulting system of linear equations, available
 *       in the @c A and @c b objects of the @c FSH data object,
 *       using some linear solver software.  The solution should
 *       be placed in the @c x object of the @c FSH object.
 *    -# Reconstruct derived quantities such as cell pressures and
 *       interface fluxes using function fsh_press_flux().
 *       Function fsh_press_flux() relies on the solution to the
 *       linear system being stored in the @c x object.
 * -# Release resources using function fsh_destroy() at end of
 *    simulation.
 */

#include <opm/core/grid.h>
#include <opm/core/pressure/legacy_well.h>
#include <opm/core/pressure/flow_bc.h>

#ifdef __cplusplus
extern "C" {
#endif


struct CSRMatrix;
struct fsh_impl;

/**
 * Main data structure of hybridized pressure solvers based on Schur
 * complement reductions.  Mainly intended to present a common view
 * of a Schur complement system of simultaneous linear equations
 * and to hold the solution of said system.
 */
struct fsh_data {
    /**
     * Maximum number of connections in any grid cell,
     * \f[
     * \mathit{max\_ngconn} = \max_c \{ n_c \}
     * \f]
     * in which \f$n_c\f$ denotes the number connections
     * (i.e., faces) of cell \f$c\f$.
     */
    int               max_ngconn;

    /**
     * Sum of squared number of connections in all grid cells,
     * \f[
     * \mathit{sum\_ngconn2} = \sum_c n_c^2.
     * \f]
     */
    size_t            sum_ngconn2;

    /* Linear system */
    struct CSRMatrix *A;        /**< Coefficient matrix */
    double           *b;        /**< System RHS */
    double           *x;        /**< Solution */

    /** Private implementational details. */
    struct fsh_impl  *pimpl;
};



/**
 * Dispose of all memory associated to <CODE>FSH</CODE> object.
 *
 * @param[in,out] h <CODE>FSH</CODE> object.  Following a call
 *                  to function fsh_destroy(), the pointer is
 *                  invalid.
 */
void
fsh_destroy(struct fsh_data *h);


/**
 * Construct compressible hybrid flow-solver data object for a
 * given grid and well pattern.
 *
 * @param[in] G Grid.
 * @param[in] W Well topology.  Ignored.
 * @return Fully formed data object suitable for use in a
 * compressible pressure solver.  @c NULL in case of construction
 * failure.
 */
struct fsh_data *
cfsh_construct(struct UnstructuredGrid *G, well_t *W);



/**
 * Form Schur-complement system of simultaneous linear equations
 * arising in compressible flow using a hybridized formulation.
 *
 * Upon returning from function cfsh_assemble(), the resulting
 * system of simultaneous linear equations is stored in
 * <CODE>h->A</CODE> and <CODE>h->b</CODE>.
 *
 * @param[in]     bc     Boundary conditions.
 * @param[in]     src    Explicit source terms.
 * @param[in]     Binv   Inverse of block-diagonal matrix \f$B\f$
 *                       Typically computed using functions
 *                       mim_ip_simple_all() and
 *                       mim_ip_mobility_update().
 * @param[in]     Biv    \f$B^{-1}v\f$.
 * @param[in]     P      Compressible accumulation term.
 * @param[in]     gpress Gravity pressure.
 * @param[in]     wctrl  Well controls.  Ignored.
 * @param[in]     WI     Well indices.  Ignored.
 * @param[in]     BivW   \f$B^{-1}v\f$ for wells.  Ignored.
 * @param[in]     wdp    Gravity pressure along well track.  Ignored.
 * @param[in,out] h      Data object.
 */
void
cfsh_assemble(struct FlowBoundaryConditions *bc,
              const double    *src,
              const double    *Binv,
              const double    *Biv,
              const double    *P,
              const double    *gpress,
              well_control_t  *wctrl,
              const double    *WI,
              const double    *BivW,
              const double    *wdp,
              struct fsh_data *h);


/**
 * Construct incompressible hybrid flow-solver data object for a
 * given grid and well pattern.
 *
 * @param G Grid.
 * @param W Well topology.
 * @return Fully formed data object suitable for use in an
 * incompressible pressure solver.  @c NULL in case of construction
 * failure.
 */
struct fsh_data *
ifsh_construct(struct UnstructuredGrid *G, well_t *W);


/**
 * Form Schur-complement system of simultaneous linear equations
 * arising in compressible flow using a hybridized formulation.
 *
 * Upon returning from function cfsh_assemble(), the resulting
 * system of simultaneous linear equations is stored in
 * <CODE>h->A</CODE> and <CODE>h->b</CODE>.
 *
 * @param[in]     bc     Boundary conditions.
 * @param[in]     src    Explicit source terms.
 * @param[in]     Binv   Inverse of block-diagonal matrix \f$B\f$
 *                       Typically computed using functions
 *                       mim_ip_simple_all() and
 *                       mim_ip_mobility_update().
 * @param[in]     gpress Gravity pressure.
 * @param[in]     wctrl  Well controls.
 * @param[in]     WI     Well indices.
 * @param[in]     wdp    Gravity pressure along well track.
 * @param[in,out] h      Data object.
 */
void
ifsh_assemble(struct FlowBoundaryConditions *bc,
              const double    *src,
              const double    *Binv,
              const double    *gpress,
              well_control_t  *wctrl,
              const double    *WI,
              const double    *wdp,
              struct fsh_data *h);




/**
 * Compute cell pressures (cpress) and interface fluxes (fflux) from
 * current solution of system of linear equations, <CODE>h->x</CODE>.
 * Back substitution process, projected half-contact fluxes.
 *
 * @param[in]  G      Grid.
 * @param[in]  Binv   Inverse of block-diagonal matrix \f$B\f$
 *                    Must coincide with the matrix used to
 *                    form the system of linear equations
 *                    currently stored in the data object.
 * @param[in]  gpress Gravity pressure.  Must coincide with
 *                    the array used to form the system of
 *                    linear equations.
 * @param[in]  h      Data object.
 * @param[out] cpress Cell pressure.
 * @param[out] fflux  Interface fluxes.
 * @param[out] wpress Well pressure.
 * @param[out] wflux  Well perforation fluxes.
 */
void
fsh_press_flux(struct UnstructuredGrid *G,
               const double *Binv, const double *gpress,
               struct fsh_data *h,
               double *cpress, double *fflux,
               double *wpress, double *wflux);

#ifdef __cplusplus
}
#endif


#endif  /* OPM_FSH_HEADER_INCLUDED */
