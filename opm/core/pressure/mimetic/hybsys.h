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

#ifndef OPM_HYBSYS_HEADER_INCLUDED
#define OPM_HYBSYS_HEADER_INCLUDED

/**
 * \file
 * Routines and data structures to manage local contributions to a
 * global system of simultaneous linear equations arising from a
 * Schur complement reduction of an original block system.
 *
 * Specifically, these data structures and related routines compute
 * and store the elemental (cell-based) contributions of the Schur
 * complement reduction of the block system of simultaneous linear
 * equations
 * \f[
 * \begin{pmatrix}
 * B              & C_1 & D \\
 * C_2^\mathsf{T} & P   & 0 \\
 * D^\mathsf{T}   & 0   & 0
 * \end{pmatrix}
 * \begin{pmatrix}
 * v \\ -p \\ \pi
 * \end{pmatrix} = \begin{pmatrix}
 * G \\ g \\ h
 * \end{pmatrix}
 * \f]
 * in which \f$G\f$ accounts for effects of gravity.  The traditional
 * Schurcomplement reduction (block Gaussian elimination) then produces
 * the equivalent system of simultaneous linear equations
 * \f[
 * \begin{pmatrix}
 * B & C_1 &  D   \\
 * 0 & -L  & -F_2 \\
 * 0 &  0  &  A
 * \end{pmatrix}
 * \begin{pmatrix}
 * v \\ -p \\ \pi
 * \end{pmatrix} = \begin{pmatrix}
 * G \\ g - C_2^\mathsf{T}B^{-1}G \\ b
 * \end{pmatrix}.
 * \f]
 * Here, the matrix \f$A\f$ and the right hand side vector \f$b\f$ are given
 * by
 * \f[
 * \begin{aligned}
 * A &= D^\mathsf{T}B^{-1}D - F_1^\mathsf{T}L^{-1}F_2 \\
 * b &= D^\mathsf{T}B^{-1}G +
 *      F_1^\mathsf{T}L^{-1}(g - C_2^\mathsf{T}B^{-1}G) - h,
 * \end{aligned}
 * \f]
 * and the component matrices \f$F_1\f$, \f$F_2\f$, and \f$L\f$ are given
 * by
 * \f[
 *   F_1 = C_1^\mathsf{T}B^{-1}D, \quad
 *   F_2 = C_2^\mathsf{T}B^{-1}D, \quad
 *   L   = C_2^\mathsf{T}B^{-1}C_1 - P.
 * \f]
 * In the case of incompressible flow, the matrix \f$C_2\f$ is the same
 * as \f$C_1\f$ and \f$P=0\f$ whence the coefficient matrix \f$A\f$ of
 * the Schur complement system \f$A\pi=b\f$ is symmetric.
 *
 * A great deal of simplification arises from the simple characterisation
 * of the \f$C_1\f$ and \f$D\f$ matrices.  Specifically,
 * \f[
 * (C_1)_{ij} = \begin{cases}
 * 1, &\quad i\in\{\mathit{pconn}_j, \dots, \mathit{pconn}_{j+1}-1\}, \\
 * 0, &\quad \text{otherwise},
 * \end{cases}
 * \f]
 * and
 * \f[
 * (D)_{ij} = \begin{cases}
 * 1, &\quad \mathit{conn}_i = j, \\
 * 0, &\quad \text{otherwise}.
 * \end{cases}
 * \f]
 * When viewed in the context of a single cell, then the \f$D\f$ matrix
 * is, effectively, the identity with the \f$\mathit{conn}\f$ array
 * simply affecting a finite-element style redistribution (assembly)
 * of the local contributions.  This module leverages that property
 * extensively.
 */

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Elemental contributions (from cells) to block system of simultaneous
 * linear equations.  Mixes quantities of single cells (@c r,
 * @c S and @c one) with those that cater to all cells (@c L,
 * @c q, @c F1, and--possibly--@c F2).
 */
struct hybsys {
    double *L;    /**< \f$C_2^\mathsf{T}B^{-1}C - P\f$, all cells */
    double *q;    /**< \f$g - F_2 G\f$, all cells */
    double *F1;   /**< \f$C_1^\mathsf{T}B^{-1}\f$, all cells */
    double *F2;   /**< \f$C_2^\mathsf{T}B^{-1}\f$, all cells*/
    double *r;    /**< Data buffer for system right-hand side, single cell */
    double *S;    /**< Data buffer system matrix, single cell */
    double *one;  /**< \f$(1,1,\dots,1)^\mathsf{T}\f$, single cell */
};


/**
 * Elemental contributions (from wells) to block system of simultaneous
 * linear equations.  Mixes quantities of single cell connections (@c r,
 * @c w2r, @c r2w, and @c w2w) and those that pertain to all well
 * connections (perforations) in concert (@c F1 and @c F2).
 */
struct hybsys_well {
    double *F1;    /**< \f$C_1^\mathsf{T}B^{-1}\f$, all connections. */
    double *F2;    /**< \f$C_2^\mathsf{T}B^{-1}\f$, all connections. */
    double *r;     /**< Data buffer for system right-hand side, single cell. */

    double *w2r;   /**< Well-to-reservoir connection strength, single cell. */
    double *r2w;   /**< Reservoir-to-well connection strength, single cell. */
    double *w2w;   /**< Aggregate well-to-well connection strength. */

    double *data;  /**< Linear storage array.  Structure undisclosed. */
};


/**
 * Allocate a hybrid system management structure suitable for discretising
 * a symmetric (i.e., incompressible) flow problem on a grid model of
 * given size.
 *
 * @param[in] max_nconn Maximum number of single cell faces.
 * @param[in] nc        Total number of grid cells.
 * @param[in] nconn_tot Aggregate number of cell faces for all cells.
 * @return Fully formed hybrid system management structure if successful or
 * @c NULL in case of allocation failure.
 */
struct hybsys *
hybsys_allocate_symm(int max_nconn, int nc, int nconn_tot);


/**
 * Allocate a hybrid system management structure suitable for discretising
 * an unsymmetric (i.e., compressible) flow problem on a grid model of
 * given size.
 *
 * @param[in] max_nconn Maximum number of single cell faces.
 * @param[in] nc        Total number of grid cells.
 * @param[in] nconn_tot Aggregate number of cell faces for all cells.
 * @return Fully formed hybrid system management structure if successful or
 * @c NULL in case of allocation failure.
 */
struct hybsys *
hybsys_allocate_unsymm(int max_nconn, int nc, int nconn_tot);


/**
 * Allocate a hybrid system management structure suitable for discretising
 * an incompressible (i.e., symmetric) well flow problem on a grid model
 * of given size.
 *
 * @param[in] max_nconn Maximum number of single cell faces.
 * @param[in] nc        Total number of grid cells.
 * @param[in] cwpos     Indirection array that defines each cell's
 *                      connecting wells.  Values typically computed
 *                      using function derive_cell_wells().
 * @return Fully formed hybrid system management structure if successful or
 * @c NULL in case of allocation failure.
 */
struct hybsys_well *
hybsys_well_allocate_symm(int max_nconn, int nc, int *cwpos);


/**
 * Allocate a hybrid system management structure suitable for discretising
 * a compressible (i.e., unsymmetric) well flow problem on a grid model
 * of given size.
 *
 * @param[in] max_nconn Maximum number of single cell faces.
 * @param[in] nc        Total number of grid cells.
 * @param[in] cwpos     Indirection array that defines each cell's
 *                      connecting wells.  Values typically computed
 *                      using function derive_cell_wells().
 * @return Fully formed hybrid system management structure if successful
 * or @c NULL in case of allocation failure.
 */
struct hybsys_well *
hybsys_well_allocate_unsymm(int max_nconn, int nc, int *cwpos);


/**
 * Dispose of memory resources previously obtained through one of the
 * allocation functions, hybsys_allocate_symm() or
 * hybsys_allocate_unsymm().
 *
 * Following a call to hybsys_free(), the input pointer is no longer
 * valid.  <CODE>hybsys_free(NULL)</CODE> does nothing.
 *
 * @param[in,out] sys Previously allocated hybrid system management
 *                    structure (or @c NULL).
 */
void
hybsys_free(struct hybsys *sys);

/**
 * Dispose of memory resources previously obtained through one of the
 * allocation functions, hybsys_well_allocate_symm() or
 * hybsys_well_allocate_unsymm().
 *
 * Following a call to hybsys_well_free(), the input pointer is
 * no longer valid.  <CODE>hybsys_well_free(NULL)</CODE> does nothing.
 *
 * @param[in,out] wsys Previously allocated hybrid system management
 *                     structure (or @c NULL).
 */
void
hybsys_well_free(struct hybsys_well *wsys);


/**
 * Perform post-construction dynamic initialisation of system
 * structure obtained from function hybsys_allocate_symm() or
 * hybsys_allocate_unsymm().
 *
 * @param[in] max_nconn Maximum number of single cell faces.
 *                      Must coincide with the equally named
 *                      parameter of functions hybsys_allocate_symm()
 *                      or hybsys_allocate_unsymm().
 * @param[in,out] sys   Previously allocated hybrid system management
 *                      structure.
 */
void
hybsys_init(int max_nconn, struct hybsys *sys);

/**
 * Compute elemental (per-cell) contributions to symmetric Schur
 * system of simultaneous linear equations.
 *
 * This function assumes that the coefficient matrix of the hybrid
 * system of linear equations is that of the introduction with the
 * additional provision that \f$C_1=C_2=C\f$ and that \f$P=0\f$.
 * In other words, this function assumes that the coefficient matrix
 * is of the form
 * \f[
 * \begin{pmatrix}
 * B            & C & D \\
 * C^\mathsf{T} & 0 & 0 \\
 * D^\mathsf{T} & 0 & 0
 * \end{pmatrix}.
 * \f]
 * This function fills the @c F1 and @c L fields of the management
 * structure.
 *
 * @param[in]     nc    Total number of grid cells.
 * @param[in]     pconn Cell-to-face start pointers.
 * @param[in]     Binv  Inverse inner product results, usually
 *                      computed using mim_ip_simple_all() and
 *                      mim_ip_mobility_update().
 * @param[in,out] sys   Hybrid system management structure allocated
 *                      using hybsys_allocate_symm() and initialised
 *                      using hybsys_init().
 */
void
hybsys_schur_comp_symm(int nc, const int *pconn,
                       const double *Binv, struct hybsys *sys);


/**
 * Compute elemental (per-cell) contributions to unsymmetric Schur
 * system of simultaneous linear equations.
 *
 * This function assumes that the coefficient matrix of the hybrid
 * system of linear equations is that of the introduction with the
 * additional provision that \f$C_2=C_1-V\f$.  In other words, this
 * function assumes that the coefficient matrix is of the form
 * \f[
 * \begin{pmatrix}
 * B                & C & D \\
 * (C-V)^\mathsf{T} & P & 0 \\
 * D^\mathsf{T}     & 0 & 0
 * \end{pmatrix}.
 * \f]
 * This matrix arises in the ``\f$v^2\f$'' phase compressibility
 * formulation of the compressible black-oil model. This function
 * fills the @c F1, @c F2 and @c L fields of the management structure.
 *
 * @param[in]     nc    Total number of grid cells.
 * @param[in]     pconn Cell-to-face start pointers.
 * @param[in]     Binv  Inverse inner product results, usually
 *                      computed using mim_ip_simple_all() and
 *                      mim_ip_mobility_update().
 * @param[in]     BIV   \f$B^{-1}v\f$ in which \f$v\f$ is the flux
 *                      field of a previous time step or non-linear
 *                      iteration.
 * @param[in]     P     Per cell compressible accumulation term.  One
 *                      scalar per cell.
 * @param[in,out] sys   Hybrid system management structure allocated
 *                      using hybsys_allocate_symm() and initialised
 *                      using hybsys_init().
 */
void
hybsys_schur_comp_unsymm(int nc, const int *pconn,
                         const double *Binv, const double *BIV,
                         const double *P, struct hybsys *sys);


/**
 * Compute elemental (per-cell) contributions to unsymmetric Schur
 * system of simultaneous linear equations.
 *
 * This function assumes that the coefficient matrix of the hybrid
 * system of linear equations is that of the introduction with no
 * additional provisions.  In other words, this
 * function assumes that the coefficient matrix is of the form
 * \f[
 * \begin{pmatrix}
 * B              & C_1 & D \\
 * C_2^\mathsf{T} & P   & 0 \\
 * D^\mathsf{T}   & 0   & 0
 * \end{pmatrix}.
 * \f]
 * This function fills the @c F1, @c F2 and @c L fields of
 * the management structure.
 *
 * @param[in]     nc    Total number of grid cells.
 * @param[in]     pconn Cell-to-face start pointers.
 * @param[in]     Binv  Inverse inner product results, usually
 *                      computed using mim_ip_simple_all() and
 *                      mim_ip_mobility_update().
 * @param[in]     C2    Explicit representation of the \f$C_2\f$
 *                      matrix as a linear array.  Assumed to only
 *                      contain the (structurally) non-zero matrix
 *                      elements (that correspond to the non-zero
 *                      structure of \f$C_1\f$).
 * @param[in]     P     Per cell compressible accumulation term.  One
 *                      scalar per cell.
 * @param[in,out] sys   Hybrid system management structure allocated
 *                      using hybsys_allocate_symm() and initialised
 *                      using hybsys_init().
 */
void
hybsys_schur_comp_gen(int nc, const int *pconn,
                      const double *Binv, const double *C2,
                      const double *P, struct hybsys *sys);

/**
 * Compute elemental contributions to global, symmetric system of
 * simultaneous linear equations from cell<->well connections.
 *
 * Specifically, for a well @c w intersecting a cell @c c, this function
 * computes the elemental contributions
 * \f[
 * (F_1)_{wc} = C_{wc}^\mathsf{T} B_{wc}^{-1} D_{wc} = \mathit{WI}_{wc}
 * \f]
 * and
 * \f[
 * L_{wc} = C_{wc}^\mathsf{T} B_{wc}^{-1} C_{wc} = \mathit{WI}_{wc}
 * \f]
 * and incorporates the contributions into the global system quantities
 * as appropriate.
 *
 * This function modifies <CODE>sys->L</CODE> and <CODE>wsys->F1</CODE>.
 *
 * @param[in]     nc    Total number of grid cells.
 * @param[in]     cwpos Indirection array that defines each cell's
 *                      connecting wells.  Values typically computed
 *                      using function derive_cell_wells().
 * @param[in]     WI    Peaceman well connection indices.  Array of
 *                      size <CODE>cwpos[nc]</CODE>.  Must incorporate
 *                      effects of multiple phases (i.e., total mobility)
 *                      if applicable.
 * @param[in,out] sys   Hybrid system management structure allocated
 *                      using hybsys_allocate_symm() and initialised
 *                      using hybsys_init() and/or filled using function
 *                      hybsys_schur_comp_symm().
 * @param[in,out] wsys  Hybrid well-system management structure obtained
 *                      from function hybsys_well_allocate_symm().
 */
void
hybsys_well_schur_comp_symm(int nc, const int *cwpos,
                            double             *WI,
                            struct hybsys      *sys,
                            struct hybsys_well *wsys);

/**
 * Compute final (symmetric) Schur complement contributions to
 * global system of simultaneous linear equations.
 *
 * This function forms the coefficient matrix
 * \f[
 * S_c = D^\mathsf{T}B_c^{-1}D - F_c^\mathsf{T}L_c^{-1}F_c
 * \f]
 * and similar right-hand side \f$r_c\f$ elemental contributions.
 * These values must be subsequently assembled into the global system
 * using function hybsys_global_assemble_cell() after imposing any
 * applicable boundary conditions.
 *
 * This function overwrites the fields @c S and @c r of the hybrid system
 * structure.
 *
 * @param[in]     c      Cell for which to compute local contributions.
 * @param[in]     nconn  Number of connections (faces) of cell @c c.
 * @param[in]     p1     Start address (into @c gpress) of the gravity
 *                       contributions of cell @c c.
 * @param[in]     p2     Start address (into @c Binv) of the inverse
 *                       inner product of cell @c c.
 * @param[in]     gpress Gravity contributions of all cells.  Must
 *                       include effects of multiple phases if applicable.
 * @param[in]     src    Explicit source terms for all cells.
 * @param[in]     Binv   Inverse inner products for all cells.  Must
 *                       include effects of multiple phases if applicable.
 * @param[in,out] sys    Hybrid system management structure allocated
 *                       using hybsys_allocate_symm() and initialised
 *                       using hybsys_init() and/or filled using function
 *                       hybsys_schur_comp_symm() and
 *                       hybsys_well_schur_comp_symm() if applicable.
 */
void
hybsys_cellcontrib_symm(int c, int nconn, int p1, int p2,
                        const double *gpress, const double *src,
                        const double *Binv, struct hybsys *sys);


/**
 * Compute final (non-symmetric) Schur complement contributions to
 * global system of simultaneous linear equations.
 *
 * This function forms the coefficient matrix
 * \f[
 * S_c = D^\mathsf{T}B_c^{-1}D - (F_1)_c^\mathsf{T}L_c^{-1}(F_2)_c
 * \f]
 * and similar right-hand side \f$r_c\f$ elemental contributions.
 * These values must be subsequently assembled into the global system
 * using function hybsys_global_assemble_cell() after imposing any
 * applicable boundary conditions.
 *
 * This function overwrites the fields @c S and @c r of the hybrid system
 * structure.
 *
 * @param[in]     c      Cell for which to compute local contributions.
 * @param[in]     nconn  Number of connections (faces) of cell @c c.
 * @param[in]     p1     Start address (into @c gpress) of the gravity
 *                       contributions of cell @c c.
 * @param[in]     p2     Start address (into @c Binv) of the inverse
 *                       inner product of cell @c c.
 * @param[in]     gpress Gravity contributions of all cells.  Must
 *                       include effects of multiple phases if applicable.
 * @param[in]     src    Explicit source terms for all cells.
 * @param[in]     Binv   Inverse inner products for all cells.  Must
 *                       include effects of multiple phases if applicable.
 * @param[in,out] sys    Hybrid system management structure allocated
 *                       using hybsys_allocate_symm() and initialised
 *                       using hybsys_init() and/or filled using functions
 *                       hybsys_schur_comp_unsymm() or hybsys_schur_comp_gen().
 */
void
hybsys_cellcontrib_unsymm(int c, int nconn, int p1, int p2,
                          const double *gpress, const double *src,
                          const double *Binv, struct hybsys *sys);

void
hybsys_well_cellcontrib_symm(int c, int ngconn, int p1,
                             const int *cwpos,
                             const double *WI, const double *wdp,
                             struct hybsys *sys, struct hybsys_well *wsys);


/**
 * Recover cell pressures and outward fluxes (with respect to cells--i.e., the
 * ``half-face fluxes'') through back substitution after solving a symmetric
 * (i.e., incompressible) Schur complement system of simultaneous linear
 * equations.
 *
 * Specifically, given the solution \f$\pi\f$ to the global system of
 * simultaneous linear equations, \f$A\pi=b\f$, that arises as a result of the
 * Schur complement analysis, this function recovers the cell pressures \f$p\f$
 * and outward fluxes \f$v\f$ defined by
 * \f[
 * \begin{aligned}
 * Lp &= g - C_2^\mathsf{T}B^{-1}G + F_2\pi \\
 * Bv &= G + C_1p - D\pi
 * \end{aligned}.
 * \f]
 *
 * @param[in]     nc     Total number of grid cells.
 * @param[in]     pconn  Cell-to-face start pointers.
 * @param[in]     conn   Cell-to-face mapping.
 * @param[in]     gpress Gravity contributions of all cells.  Must coincide with
 *                       equally named parameter in calls to cell contribution
 *                       functions such as hybsys_cellcontrib_symm().
 * @param[in]     Binv   Inverse inner products for all cells.  Must coincide
 *                       with equally named parameter in calls to contribution
 *                       functions such as hybsys_cellcontrib_symm().
 * @param[in]     sys    Hybrid system management structure coinciding with
 *                       equally named parameter in contribution functions such
 *                       as hybsys_cellcontrib_symm() or
 *                       hybsys_cellcontrib_unsymm().
 * @param[in]     pi     Solution (interface/contact pressure) obtained from
 *                       solving the global system \f$A\pi = b\f$.
 * @param[out]    press  Cell pressures, \f$p\f$.  Array of size @c nc.
 * @param[out]    flux   Outward interface fluxes, \f$v\f$.  Array of size
 *                       <CODE>pconn[nc]</CODE>.
 * @param[in,out] work   Scratch array for temporary results.  Array of size at
 *                       least \f$\max_c \{   \mathit{pconn}_{c + 1}
 *                                          - \mathit{pconn}_c \} \f$.
 */
void
hybsys_compute_press_flux(int nc, const int *pconn, const int *conn,
                          const double *gpress,
                          const double *Binv, const struct hybsys *sys,
                          const double *pi, double *press, double *flux,
                          double *work);


/**
 * Recover well pressures (i.e., bottom-hole pressure values) and well
 * connection (perforation) fluxes.
 *
 * Specifically, this function performs the same role (i.e., back-substitution)
 * for wells as function hybsys_compute_press_flux() does for grid cells and
 * grid contacts (interfaces).
 *
 * @param[in]     nc     Total number of grid cells.
 * @param[in]     pgconn Cell-to-face start pointers.
 * @param[in]     nf     Total number of grid faces.
 * @param[in]     nw     Total number of wells.
 * @param[in]     pwconn Cell-to-well start pointers.  If <CODE>nw > 0</CODE>,
 *                       then this parameter must coincide with the @c cwpos
 *                       array used in call to hybsys_well_schur_comp_symm().
 * @param[in]     wconn  Cell-to-well mapping.
 * @param[in]     Binv   Inverse inner products for all cells.  Must coincide
 *                       with equally named parameter in calls to contribution
 *                       functions such as hybsys_well_cellcontrib_symm().
 * @param[in]     WI     Peaceman well connection indices.  Array of
 *                       size <CODE>pwconn[nc]</CODE>.  Must coincide with
 *                       equally named parameter of contribution function
 *                       hybsys_well_cellcontrib_symm().
 * @param[in]     wdp    Well connection gravity pressure adjustments.
 * @param[in]     sys    Hybrid system management structure coinciding with
 *                       equally named parameter in contribution functions such
 *                       as hybsys_cellcontrib_symm() and
 *                       hybsys_well_cellcontrib_symm().
 * @param[in]     wsys   Hybrid well-system management structure.  Must coincide
 *                       with equally named paramter of contribution function
 *                       hybsys_well_cellcontrib_symm().
 * @param[in]     pi     Solution (interface/contact pressure and well BHPs)
 *                       obtained from solving the global system \f$A\pi = b\f$.
 * @param[in]     cpress Cell pressures, \f$p\f$, obtained from a previous call
 *                       to function hybsys_compute_press_flux().
 * @param[in]     cflux  Outward fluxes, \f$v\f$, obtained from a previous call
 *                       to function hybsys_compute_press_flux().
 * @param[out]    wpress Well (i.e., bottom-hole) pressures.  Array of size
 *                       @c nw.
 * @param[out]    wflux  Well connection (perforation) fluxes.  Array of size
 *                       <CODE>pwconn[nw]</CODE>.
 * @param[in,out] work   Scratch array for storing intermediate results.  Array
 *                       of size at least \f$\max_w \{  \mathit{pwconn}_{w + 1}
 *                                                    - \mathit{pwconn}_w\}\f$.
 */
void
hybsys_compute_press_flux_well(int nc, const int *pgconn, int nf,
                               int nw, const int *pwconn, const int *wconn,
                               const double *Binv,
                               const double *WI,
                               const double *wdp,
                               const struct hybsys      *sys,
                               const struct hybsys_well *wsys,
                               const double             *pi,
                               double *cpress, double *cflux,
                               double *wpress, double *wflux,
                               double *work);


#ifdef __cplusplus
}
#endif

#endif  /* OPM_HYBSYS_HEADER_INCLUDED */
