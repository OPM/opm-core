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
#define OPM_FHS_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/well.h>
#include <opm/core/pressure/flow_bc.h>

#ifdef __cplusplus
extern "C" {
#endif


/***************************************************************/
/* Data type common to compressible and incompressible solver. */
/***************************************************************/

struct CSRMatrix;
struct fsh_impl;

/** Contains the linear system for assembly, as well as internal data
 *  for the assembly routines.
 */
struct fsh_data {
    /* Let \f$n_i\f$ be the number of connections/faces of grid cell
     * number \f$i\f$. Then max_ngconn = \f$\max_i n_i\f$
     */
    int               max_ngconn;
    /* With n_i as above, sum_ngconn2 = \f$\sum_i n_i^2\f$ */
    size_t            sum_ngconn2;

    /* Linear system */
    struct CSRMatrix *A;        /* Coefficient matrix */
    double           *b;        /* System RHS */
    double           *x;        /* Solution */

    /* Private implementational details. */
    struct fsh_impl  *pimpl;
};



/** Destroys the fsh data object */
void
fsh_destroy(struct fsh_data *h);


/*********************************/
/* Compressible solver routines. */
/*********************************/


/** Constructs compressible hybrid flow-solver data object for a
 *  given grid and well pattern.
 */
struct fsh_data *
cfsh_construct(struct UnstructuredGrid *G, well_t *W);



/** Assembles the hybridized linear system for face pressures.
 */
void
cfsh_assemble(flowbc_t        *bc,
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


/***********************************/
/* Incompressible solver routines. */
/***********************************/
/** Constructs incompressible hybrid flow-solver data object for a
 *  given grid and well pattern.
 *
 *  @param G The grid
 *  @param W The wells
 */
struct fsh_data *
ifsh_construct(struct UnstructuredGrid *G, well_t *W);



/** Assembles the hybridized linear system for face pressures.
 *
 * This routine produces no output, other than changing the linear
 * system embedded in the ifsh_data object.
 * @param bc Boundary conditions.
 * @param src Per-cell source terms (volume per second). Positive
 *            values flow are sources, negative values are sinks.
 * @param Binv The cell-wise effective inner products to employ in
 *             assembly. This should be an array of length equal to
 *             sum_ngconn2 of the ifsh_data object. For each cell i,
 *             there are \f$n_i^2\f$ entries, giving the inner product for
 *             that cell. The inner products may for example be
 *             computed by the functions of mimetic.h.
 * @param gpress Effective gravity terms. This should be an array of length
 *               \f$\sum_i n_i\f$. For each cell, the \f$n_i\f$ elements
 *               corresponding to cell \f$i\f$ should be given by
 *               \f$\omega g \cdot (f_c - c_c)\f$ where the symbols
 *               represent the fractional-flow-weighted densities,
 *               the gravity vector, face centroid and cell centroid.
 * @param wctrl \TODO
 * @param WI \TODO
 * @param wdp \TODO
 * @param h The fsh_data object to use (and whose linear system will
 *          be modified). Must already be constructed.
 */
void
ifsh_assemble(flowbc_t        *bc,
              const double    *src,
              const double    *Binv,
              const double    *gpress,
              well_control_t  *wctrl,
              const double    *WI,
              const double    *wdp,
              struct fsh_data *h);




/**********************************/
/* Common postprocessing routine. */
/**********************************/

/** Computes cell pressures, face fluxes, well pressures and well
 * fluxes from face pressures.
 *
 * @param G The grid.
 * @param h The fsh_data object. You must have called [ic]fsh_assemble()
 *          prior to this, and solved the embedded linear system of
 *          this object before you call fsh_press_flux().
 * @param cpress[out] Cell pressures.
 * @param fflux[out] Oriented face fluxes.
 * @param wpress[out] \TODO
 * @param wflux[out] \TODO
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
