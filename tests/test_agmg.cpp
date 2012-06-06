/*
  Copyright 2012 SINTEF ICT, Applied Mathematics.
  Copyright 2012 Statoil ASA.

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

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/linalg/LinearSolverAGMG.hpp>

namespace {
    std::size_t
    compute_nnz(const std::size_t m)
    {
        assert (m > 0);

        std::size_t nnz = m;              // A(i,i)

        nnz += (m > 1) ? 2           : 0; // A(0,1), A(m-1,m-2)
        nnz += (m > 2) ? 2 * (m - 2) : 0; // A(i,i-1), A(i,i+1)

        return nnz;
    }


    boost::shared_ptr<CSRMatrix>
    build_laplace_1d(const std::size_t m)
    {
        assert (m >= 2);

        const std::size_t nnz = compute_nnz(m);

        boost::shared_ptr<CSRMatrix>
            A(csrmatrix_new_known_nnz(m, nnz), csrmatrix_delete);

        A->ia[ 0 ] = 0;

        // First row
        A->ia[       0 + 1     ] = A->ia[ 0 ];
        A->ja[ A->ia[0 + 1]    ] = 0 + 0;
        A->sa[ A->ia[0 + 1] ++ ] =   2.0;
        A->ja[ A->ia[0 + 1]    ] = 0 + 1;
        A->sa[ A->ia[0 + 1] ++ ] = - 1.0;

        // General rows
        for (std::size_t i = 1; i < m - 1; ++i) {
            A->ia[i + 1] = A->ia[i];

            A->ja[ A->ia[i + 1]    ] = int(i) - 1;
            A->sa[ A->ia[i + 1] ++ ] = - 1.0;

            A->ja[ A->ia[i + 1]    ] = int(i)    ;
            A->sa[ A->ia[i + 1] ++ ] =   2.0;

            A->ja[ A->ia[i + 1]    ] = int(i) + 1;
            A->sa[ A->ia[i + 1] ++ ] = - 1.0;
        }

        // Last row
        A->ia[ (m - 1) + 1 ] = A->ia[ m - 1 ];

        A->ja[ A->ia[ (m - 1) + 1 ]    ] = int(m - 1) - 1;
        A->sa[ A->ia[ (m - 1) + 1 ] ++ ] = - 1.0;

        A->ja[ A->ia[ (m - 1) + 1 ]    ] = int(m - 1)    ;
        A->sa[ A->ia[ (m - 1) + 1 ] ++ ] =   2.0;

        return A;
    }
}


int main()
{
    const std::size_t m = 10;

    boost::shared_ptr<CSRMatrix> A = build_laplace_1d(m);

    // Form right-hand side [1, 0, 0, ...., 0, 1]
    std::vector<double> b(m, 0.0);
    b[0] = 1.0; b.back() = 1.0;

    // Allocate solution vector
    std::vector<double> x(m);

    // Create solver for SPD system.
    Opm::LinearSolverAGMG linsolve(100, 1e-9, true);

    Opm::LinearSolverInterface::LinearSolverReport
        rpt = linsolve.solve(A.get(), & b[0], & x[0]);

    double e = 0.0;
    for (std::size_t i = 0; i < m; ++i) {
        const double d = x[i] - 1.0;
        e += d * d;
    }

    std::cerr << "|| e ||_2 = "
              << std::scientific
              << std::setprecision(5)
              << std::sqrt(e) / double(m) << '\n';

    return 0;
}
