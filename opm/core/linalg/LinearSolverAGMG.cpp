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


#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <opm/core/linalg/LinearSolverAGMG.hpp>
#include <opm/core/utility/ErrorMacros.hpp>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <stdexcept>
#include <vector>

#if HAVE_AGMG
// Manual prototype of main gateway routine to DOUBLE PRECISION,
// serial version of the AGMG software package.
//
// Note that both the matrix entries and column indices are writable.
// The solver may permute the matrix entries within each row during
// the setup phase.
#define DAGMG_ FC_FUNC(dagmg, DAGMG)

extern "C"
void
DAGMG_(const int*    N     ,    // System size
       double*       sa    ,    // Non-zero elements
       int*          ja    ,    // Column indices
       const int*    ia    ,    // Row pointers
       const double* f     ,    // Right-hand side
       double*       x     ,    // Solution
       int*          ijob  ,    // Main control parameter
       int*          iprint,    // Message output unit
       int*          nrest ,    // Number of GCR restarts
       int*          iter  ,    // Maximum (and actual) number of iterations
       const double* tol   );   // Residual reduction tolerance

#endif // HAVE_AGMG

namespace Opm
{
    LinearSolverAGMG::LinearSolverAGMG(const int    max_it,
                                       const double rtol  ,
                                       const bool   is_spd)
        : max_it_(max_it),
          rtol_  (rtol)  ,
          is_spd_(is_spd)
    {
#if !HAVE_AGMG
        THROW("AGMG support is not enabled in this library");
#endif  // HAVE_AGMG
    }

    LinearSolverAGMG::~LinearSolverAGMG() {}
    
    LinearSolverInterface::LinearSolverReport
    LinearSolverAGMG::solve(const int     size    ,
                            const int     nonzeros,
                            const int*    ia      ,
                            const int*    ja      ,
                            const double* sa      ,
                            const double* rhs     ,
                            double*       solution) const
    {
        const std::vector<double>::size_type nnz = ia[size];

        assert (nnz == std::vector<double>::size_type(nonzeros));

#if defined(NDEBUG)
        // Suppress warning about unused parameter.
        static_cast<void>(nonzeros);
#endif

        std::vector<double> a(sa, sa + nnz);

        // Account for 1-based indexing.
        std::vector<int>    i(ia, ia + std::vector<int>::size_type(size + 1));
        std::transform(i.begin(), i.end(), i.begin(),
                       std::bind2nd(std::plus<int>(), 1));

        std::vector<int>    j(ja, ja + nnz);
        std::transform(j.begin(), j.end(), j.begin(),
                       std::bind2nd(std::plus<int>(), 1));

        LinearSolverInterface::LinearSolverReport rpt = {};
        rpt.iterations = max_it_;

        int ijob = 0;           // Setup + solution + cleanup, x0==0.

        int nrest;
        if (is_spd_) {
            nrest = 1;          // Use CG algorithm
        }
        else {
            nrest = 10;         // Suggested default number of GCR restarts.
        }

        int iprint = 0;         // Suppress most output
        DAGMG_(& size, & a[0], & j[0], & i[0], rhs, solution,
               & ijob, & iprint, & nrest, & rpt.iterations, & rtol_);

        rpt.converged = rpt.iterations <= max_it_;
        return rpt;
    }

    void
    LinearSolverAGMG::setTolerance(const double tol)
    {
        rtol_ = tol;
    }

    double
    LinearSolverAGMG::getTolerance() const
    {
        return rtol_;
    }
} // namespace Opm
