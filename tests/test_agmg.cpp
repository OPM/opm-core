#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <opm/core/linalg/sparse_sys.h>
#include <opm/core/linalg/LinearSolverAGMG.hpp>

namespace {
    boost::shared_ptr<CSRMatrix>
    build_laplace_1d(const std::size_t m)
    {
        const std::size_t nnz = 3*(m - 2) + 2*2;

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
