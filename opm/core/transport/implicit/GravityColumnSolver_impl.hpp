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

#include <opm/core/transport/GravityColumnSolver.hpp>
#include <opm/core/linalg/blas_lapack.h>
#include <opm/core/utility/ErrorMacros.hpp>

#include <iostream>

#include <sys/time.h>

namespace Opm
{

    template <class Model>
    GravityColumnSolver<Model>::GravityColumnSolver(Model& model,
						    const UnstructuredGrid& grid,
						    const double tol,
						    const int maxit)
	: model_(model), grid_(grid), tol_(tol), maxit_(maxit)
    {
    }

    namespace {
	struct ZeroVec
	{
	    double operator[](int) const { return 0.0; }
	};
	struct StateWithZeroFlux
	{
	    StateWithZeroFlux(std::vector<double>& s) : sat(s) {}
	    const ZeroVec& faceflux() const { return zv; }
	    const std::vector<double>& saturation() const { return sat; }
	    std::vector<double>& saturation() { return sat; }
	    ZeroVec zv;
	    std::vector<double>& sat;
	};
	struct Vecs
	{
	    Vecs(int sz) : sol(sz, 0.0) {}
	    const std::vector<double>& solution() const { return sol; }
	    std::vector<double>& writableSolution() { return sol; }
	    std::vector<double> sol;
	};
	struct JacSys
	{
	    JacSys(int sz) : v(sz) {}
	    const Vecs& vector() const { return v; }
	    Vecs& vector() { return v; }
	    Vecs v;
	    typedef std::vector<double> vector_type;
	};
    } // anon namespace



    /// \param[in] columns         for each column columns.second,
    ///                            contains the cells on which to solve the segregation
    ///                            problem. For each column, its cells must be in a single
    ///                            vertical column, and ordered
    ///                            (direction doesn't matter).
    template <class Model>
    void GravityColumnSolver<Model>::solve(const std::vector<std::vector<int> >& columns,
					   const double dt,
					   std::vector<double>& s)
    {
	// Initialize model. These things are done for the whole grid!
	StateWithZeroFlux state(s); // This holds s by reference.
	JacSys sys(grid_.number_of_cells);
	std::vector<double> increment(grid_.number_of_cells, 0.0);
	model_.initStep(state, grid_, sys);

	int iter = 0;
	double max_delta = 1e100;

	while (iter < maxit_) {
	    model_.initIteration(state, grid_, sys);
            //   std::map<int, std::vector<int> >::const_iterator it;
            //for (it = columns.begin(); it != columns.end(); ++it) {
            int size = columns.size();

            #pragma omp parallel for schedule(dynamic)
            for(int i = 0; i < size; ++i) {
		solveSingleColumn(columns[i], dt, s, increment);
	    }

	    for (int cell = 0; cell < grid_.number_of_cells; ++cell) {
		sys.vector().writableSolution()[cell] += increment[cell];
	    }

	    const double maxelem = *std::max_element(increment.begin(), increment.end());
	    const double minelem = *std::min_element(increment.begin(), increment.end());
	    max_delta = std::max(maxelem, -minelem);
	    std::cout << "Gravity column solver iteration " << iter << "   max_delta = " << max_delta << std::endl;
	    if (max_delta < tol_) {
		break;
	    }
	    ++iter;
	}
	if (max_delta >= tol_) {
	    OPM_THROW(std::runtime_error, "Failed to converge!");
	}
	// Finalize.
	// model_.finishIteration(); // Doesn't do anything in th 2p model.
	// finishStep() writes to state, which holds s by reference.
	// This will update the entire grid's state...
	model_.finishStep(grid_, sys.vector().solution(), state);
    }



    /// \param[in] column_cells    the cells on which to solve the segregation
    ///                            problem. Must be in a single vertical column,
    ///                            and ordered (direction doesn't matter).
    template <class Model>
    void GravityColumnSolver<Model>::solveSingleColumn(const std::vector<int>& column_cells,
						       const double dt,
						       std::vector<double>& s,
						       std::vector<double>& sol_vec)
    {
	// This is written only to work with SinglePointUpwindTwoPhase,
	// not with arbitrary problem models.
	const int col_size = column_cells.size();
        if (col_size == 1) {
            sol_vec[column_cells[0]] = 0.0;
            return;
        }
	StateWithZeroFlux state(s); // This holds s by reference.

	// Assemble.
	std::vector<double> tridiag_matrix_data(3*col_size - 2, 0.0);
	double* DU = &tridiag_matrix_data[0];
	double* D = DU + col_size - 1;
	double* DL = D + col_size;
	std::vector<double> rhs(col_size, 0.0);
	for (int ci = 0; ci < col_size; ++ci) {
	    double rescontrib, j1contrib, j2contrib;
	    const int cell = column_cells[ci];
	    const int prev_cell = (ci == 0) ? -999 : column_cells[ci - 1];
	    const int next_cell = (ci == col_size - 1) ? -999 : column_cells[ci + 1];
	    // model_.initResidual(cell, F);
	    for (int j = grid_.cell_facepos[cell]; j < grid_.cell_facepos[cell+1]; ++j) {
		const int face = grid_.cell_faces[j];
		const int c1 = grid_.face_cells[2*face + 0];
                const int c2 = grid_.face_cells[2*face + 1];
		if (c1 == prev_cell || c2 == prev_cell || c1 == next_cell || c2 == next_cell) {
		    j1contrib = j2contrib = rescontrib = 0.0;
		    model_.fluxConnection(state, grid_, dt, cell, face, &j1contrib, &j2contrib, &rescontrib);
		    if (c1 == prev_cell || c2 == prev_cell) {
			DL[ci-1] += j2contrib;
		    } else {
			ASSERT(c1 == next_cell || c2 == next_cell);
			DU[ci] += j2contrib;
		    }
		    D[ci] += j1contrib;
		    rhs[ci] += rescontrib;
		}
	    }
	    j1contrib = rescontrib = 0.0;
	    model_.accumulation(grid_, cell, &j1contrib, &rescontrib);
	    D[ci] += j1contrib;
	    rhs[ci] += rescontrib;
	}
	// model_.sourceTerms(); // Not needed
	// Solve.
	const MAT_SIZE_T num_rhs = 1, colSize = col_size;
	MAT_SIZE_T info = 0;
	// Solution will be written to rhs.
	dgtsv_(&colSize, &num_rhs, DL, D, DU, &rhs[0], &colSize, &info);
	if (info != 0) {
	    OPM_THROW(std::runtime_error, "Lapack reported error in dgtsv: " << info);
	}
	for (int ci = 0; ci < col_size; ++ci) {
	    sol_vec[column_cells[ci]] = -rhs[ci];
	}
    }

} // namespace Opm
