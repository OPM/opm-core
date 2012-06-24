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

#ifndef OPM_HYBSYS_GLOBAL_HEADER_INCLUDED
#define OPM_HYBSYS_GLOBAL_HEADER_INCLUDED

/**
 * \file
 * Routines to assist in the formation and assembly of a global
 * system of simultaneous linear equations derived from a Schur
 * complement reduction of an original hybrid block system.
 *
 * We assume that the original block system of linear equations
 * is given by
 * \f[
 * \begin{pmatrix}
 *   B            & C & D \\
 *   C^\mathsf{T} & 0 & 0 \\
 *   D^\mathsf{T} & 0 & 0
 * \end{pmatrix}
 * \begin{pmatrix}
 *  v \\ -p \\ \pi
 * \end{pmatrix} = \begin{pmatrix}
 *  f \\ g \\ h
 * \end{pmatrix}
 * \f]
 * in which the block matrices \f$C\f$ and \f$D\f$ are assumed to
 * have a particularly simple structure and the matrix \f$B\f$ is
 * block diagonal.
 *
 * The Schur complement reduction process (a block Gaussian elimination)
 * then produces the following block system of simultaneous linear
 * equations
 * \f[
 * \begin{pmatrix}
 * B &  C &  D \\
 * 0 & -L & -F \\
 * 0 &  0 & A
 * \end{pmatrix}
 * \begin{pmatrix}
 *  v \\ -p \\ \pi
 * \end{pmatrix} = \begin{pmatrix}
 * f \\ g - C^\mathsf{T}B^{-1} f \\ b
 * \end{pmatrix}
 * \f]
 * in which
 * \f[
 * \begin{aligned}
 *   A &= D^\mathsf{T}B^{-1}D - F^\mathsf{T}L^{-1}F \text{ and} \\
 *   b &= D^\mathsf{T}B^{-1}f + F^\mathsf{T}L^{-1} (g - C^\mathsf{T}B^{-1}f) - h.
 * \end{aligned}
 * \f]  The component matrices \f$F\f$
 * and \f$L\f$ are given by
 * \f[
 * \begin{aligned}
 *   L &= C^\mathsf{T} B^{-1} C \\
 *   F &= C^\mathsf{T} B^{-1} D.
 * \end{aligned}
 * \f]
 * The primary degrees of freedom, \f$\pi\f$, may then be recovered
 * by solving the Schur complement system
 * \f[A\pi = b\f]
 * from which the derived quantities \f$p\f$ and \f$v\f$ may be
 * computed through a back-substitution process.
 *
 * The functions in this module assist in the creation of the sparsity
 * pattern of matrix \f$A\f$ and in the global assembling of values
 * into the matrix \f$A\f$ and the right-hand side vector \f$b\f$.
 * Specifically, function hybsys_define_globconn() builds the
 * sparsity pattern of \f$A\f$ while functions hybsys_global_assemble_cell()
 * and hybsys_global_assemble_well_sym() perform the task of
 * adding matrix and right-hand side values from local contributions.
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <opm/core/grid.h>
#include <opm/core/well.h>
#include <opm/core/linalg/sparse_sys.h>


/**
 * Construct sparse matrix capable of managing a (reduced) hybrid
 * system of simultaneous linear equations with one degree of freedom
 * for each grid interface and each well.
 *
 * The return value is suitable for use in pressure solvers based on
 * Schur-complement reductions such as the mimetic finite-difference
 * or multiscale mixed finite-element classes of discretisations.  In
 * typical applications, the matrix will be cleared using function
 * csrmatrix_zero() and then filled using an assembly process similar
 * to the traditional finite-element algorithm.
 *
 * Functions hybsys_global_assemble_cell() and
 * hybsys_global_assemble_well_sys() may be used to assist such an
 * assembly process.
 *
 * @param[in] G Grid.
 * @param[in] W Well topology.  @c NULL in a model without wells.
 * @return Fully formed and structurally consistent sparse matrix
 *         whose number of rows (i.e., degrees-of-freedom) equals
 *         the number of grid faces (<CODE>G->number_of_faces</CODE>)
 *         plus the number of wells (<CODE>W->number_of_wells</CODE>).
 */
struct CSRMatrix *
hybsys_define_globconn(struct UnstructuredGrid *G, well_t *W);


/**
 * Assemble local contributions into global system of simultaneous
 * linear equations.
 *
 * The contributions will typically have been computed using
 * function hybsys_schur_comp_symm() and function
 * hybsys_cellcontrib_symm() or some of the related functions.
 *
 * @param[in]     nconn Number of cell faces.
 * @param[in]     l2g   Local-to-global mapping of cell's primary
 *                      degrees of freedom (i.e., the faces).
 * @param[in]     S     Single cell local contribution to global
 *                      coefficient matrix.  An
 *                      \f$\mathit{nconn}\times\mathit{nconn}\f$
 *                      dense matrix in column major (Fortran) order.
 * @param[in]     r     Single cell local contribution to global
 *                      system right-hand side.  An
 *                      \f$\mathit{nconn}\times 1\f$ dense vector.
 * @param[in,out] A     Global coefficient matrix (of Schur
 *                      complement system).
 * @param[in,out] b     Global system right-hand side (of Schur
 *                      complement system).
 */
void
hybsys_global_assemble_cell(int nconn, int *l2g,
                            const double     *S,
                            const double     *r,
                            struct CSRMatrix *A,
                            double           *b);


/**
 * Assemble local contributions from single cell's well connections
 * (perforations) into global system of simultaneous linear equations.
 *
 * This function assumes that the connection strength from cell to well
 * equals the connection strength from well to cell.  In other words,
 * that the numerical values of the well contributions are symmetric.
 *
 * The contributions are typically computed using functions
 * hybsys_well_schur_comp_symm() and hybsys_well_cellcontrib_symm().
 *
 * @param[in]     ngconn_tot Total number of grid connections.
 *                           Expected to equal
 *                           <CODE>G->number_of_faces</CODE>
 *                           when @c G is the grid used to form the
 *                           original matrix in
 *                           hybsys_define_globconn().
 * @param[in]     ngconn     Number of grid connections referenced by
 *                           given cell.
 * @param[in]     gconn      Actual grid connections (DOFs) referenced
 *                           by given cell.  Pointer to @c ngconn
 *                           consecutive DOF indices.
 * @param[in]     nwconn     Number of well connections intersecting
 *                           given cell.  Typically \f$\mathit{ngconn} = 1\f$.
 * @param[in]     wconn      Actual well connections (DOFs) intersecting
 *                           given cell.  Pointer to @c nwconn consecutive
 *                           DOF indices.
 * @param[in]     r2w        Reservoir-to-well connection strengths.
 * @param[in]     w2w        Well-to-well-connection strenghts.
 * @param[in]     r          Single cell local contribution to global
 *                           system right-hand side.  An
 *                           \f$(\mathit{ngconn} + \mathit{nwconn})
 *                           \times 1\f$ dense vector.
 * @param[in,out] A          Global coefficient matrix (of Schur
 *                           complement system).
 * @param[in,out] b          Global system right-hand side (of Schur
 *                           complement system).
 */
void
hybsys_global_assemble_well_sym(int ngconn_tot,
                                int ngconn, const int *gconn,
                                int nwconn, const int *wconn,
                                const double     *r2w,
                                const double     *w2w,
                                const double     *r,
                                struct CSRMatrix *A,
                                double           *b);



#ifdef __cplusplus
}
#endif

#endif  /* OPM_HYBSYS_GLOBAL_HEADER_INCLUDED */
