/*
 * Copyright (c) 2010 SINTEF ICT, Applied Mathematics
 */

#ifndef IFSH_H_INCLUDED
#define IFHS_H_INCLUDED

#include "grid.h"
#include "well.h"
#include "flow_bc.h"
#include "sparse_sys.h"


#ifdef __cplusplus
extern "C" {
#endif

/** @file Incompressible flow solver, using hybridization.
 *
 *  These functions implements assembly of a hybridized linear
 *  system for face-pressures in incompressible two-phase flow.
 *  A routine for back-substitution that computes cell pressures
 *  and face fluxes is also included.
 */



struct ifsh_impl;

/** Contains the linear system for assembly, as well as internal data
 *  for the assembly routines.
 */
struct ifsh_data {
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
    struct ifsh_impl *pimpl;
};



/** Constructs the ifsh_data object for a given grid and well
 *  pattern.
 *  @param G The grid
 *  @param W The wells
 */
struct ifsh_data *
ifsh_construct(grid_t *G, well_t *W);



/** Destroys the ifsh_data object */
void
ifsh_destroy(struct ifsh_data *h);



/** Assembles the hybridized linear system for face pressures.
 *
 * This routine produces no output, other than changing the linear
 * system embedded in the ifsh_data object.
 * @param bc Boundary conditions.
 * @param src Per-cell source terms (volume per second). Positive
 *            values flow are sources, negative values are sinks.
 * @param Binv The cell-wise inner products to employ in
 *             assembly. This should be an array of length equal to
 *             sum_ngconn2 of the ifsh_data object. For each cell i,
 *             there are \f$n_i^2\f$ entries, giving the inner product for
 *             that cell. The inner products may for example be
 *             computed by the functions of mimetic.h.
 * @param gpress Gravity terms. This should be an array of length
 *               \f$\sum_i n_i\f$. For each cell, the \f$n_i\f$ elements
 *               corresponding to cell \f$i\f$ should be given by
 *               \f$g \cdot (f_c - c_c)\f$ where the symbols represent
 *               the gravity vector, face centroid and cell centroid.
 * @param wctrl \TODO
 * @param WI \TODO
 * @param wdp \TODO
 * @param totmob Cell-wise total mobilities to use for this assembly.
 * @param omega Cell-wise phase densities weighted by fractional flow.
 * @param h The ifsh_data object to use (and whose linear system will
 *          be modified). Must already be constructed.
 */
void
ifsh_assemble(flowbc_t         *bc,
              double           *src,
              double           *Binv,
              double           *gpress,
              well_control_t   *wctrl,
              double           *WI,
              double           *wdp,
              double           *totmob, /* \sum_i \lambda_i */
              double           *omega,  /* \sum_i \rho_i f_i */
              struct ifsh_data *h);

/** Computes cell pressures, face fluxes, well pressures and well
 * fluxes from face pressures.
 *
 * @param G The grid.
 * @param h The ifsh_data object. You must have called ifsh_assemble()
 *          prior to this, and solved the embedded linear system of
 *          this object before you call ifsh_press_flux().
 * @param cpress[out] Cell pressures.
 * @param fflux[out] Oriented face fluxes.
 * @param wpress[out] \TODO
 * @param wflux[out] \TODO
 */
void
ifsh_press_flux(grid_t *G, struct ifsh_data *h,
                double *cpress, double *fflux,
                double *wpress, double *wflux);


#ifdef __cplusplus
}
#endif


#endif  /* IFSH_H_INCLUDED */
