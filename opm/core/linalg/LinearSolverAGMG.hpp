/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.

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

#ifndef OPM_LINEARSOLVERAGMG_HEADER_INCLUDED
#define OPM_LINEARSOLVERAGMG_HEADER_INCLUDED

/**
 * \file
 * Gateway to Notay's AGMG package implementing an aggregation-based
 * algebraic multigrid method.
 *
 * References:
 *   * Y. Notay,
 *     An aggregation-based algebraic multigrid method,
 *     Electronic Transactions on Numerical Analysis, vol 37,
 *     pp. 123-146, 2010.
 *
 *   * A. Napov and Y. Notay,
 *     An algebraic multigrid method with guaranteed convergence rate,
 *     Report GANMN 10-03, Université Libre de Bruxelles, Brussels,
 *     Belgium, 2010 (Revised 2011).
 *
 *   * Y. Notay,
 *     Aggregation-based algebraic multigrid for convection-diffusion
 *     equations, Report GANMN 11-01, Université Libre de Bruxelles,
 *     Brussels, Belgium, 2011.
 */

#include <opm/core/linalg/LinearSolverInterface.hpp>

namespace Opm
{
    /// Concrete class encapsulating Notay's AGMG package
    class LinearSolverAGMG : public LinearSolverInterface
    {
    public:
        /**
         * Constructor.
         * \param[in] max_it  Maximum number of solver iterations.
         * \param[in] rtol    Residual reduction tolerance.
         * \param[in] is_spd  Whether or not the matrix is SPD.  SPD
         *                    systems are solved using CG while other
         *                    systems are solved using GCR.
         */
        LinearSolverAGMG(const int    max_it = 100   ,
                         const double rtol   = 1.0e-6,
                         const bool   is_spd = false);

        /// Destructor.
        virtual ~LinearSolverAGMG();

        using LinearSolverInterface::solve;

        /// Solve a linear system, with a matrix given in compressed
        /// sparse row format.
        /// \param[in] size        Number of rows (and colums).
        /// \param[in] nonzeros    Number of (structural) non-zeros.
        /// \param[in] ia          Row pointers.
        /// \param[in] ja          Column indices.
        /// \param[in] sa          (structurally) non-zero elements.
        /// \param[in] rhs         System right-hand side.
        /// \param[inout] solution System solution.
        /// \return Solver meta-data concerning most recent system solve.
        virtual LinearSolverInterface::LinearSolverReport
        solve(const int size, const int nonzeros,
              const int* ia, const int* ja, const double* sa,
              const double* rhs, double* solution) const;

    private:
        int    max_it_;
        double rtol_  ;
        bool   is_spd_;
    };


} // namespace Opm



#endif // OPM_LINEARSOLVERAGMG_HEADER_INCLUDED
