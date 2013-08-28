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

#include "config.h"
#include <opm/core/utility/VelocityInterpolation.hpp>
#include <opm/core/grid.h>
#include <opm/core/linalg/blas_lapack.h>

#include <iostream>

namespace Opm
{

    // --------  Methods of class VelocityInterpolationInterface  --------


    VelocityInterpolationInterface::~VelocityInterpolationInterface()
    {
    }



    // --------  Methods of class VelocityInterpolationConstant  --------

    /// Constructor.
    /// \param[in]  grid   A grid.
    VelocityInterpolationConstant::VelocityInterpolationConstant(const UnstructuredGrid& grid)
        : grid_(grid)
    {
    }

    /// Set up fluxes for interpolation.
    /// \param[in]  flux   One signed flux per face in the grid.
    void VelocityInterpolationConstant::setupFluxes(const double* flux)
    {
        flux_ = flux;
    }

    /// Interpolate velocity.
    /// \param[in]  cell   Cell in which to interpolate.
    /// \param[in]  x      Coordinates of point at which to interpolate.
    ///                    Must be array of length grid.dimensions.
    /// \param[out] v      Interpolated velocity.
    ///                    Must be array of length grid.dimensions.
    void VelocityInterpolationConstant::interpolate(const int cell,
                                                    const double* /*x*/,
                                                    double* v) const
    {
        const int dim = grid_.dimensions;
        std::fill(v, v + dim, 0.0);
        const double* cc = grid_.cell_centroids + cell*dim;
        for (int hface = grid_.cell_facepos[cell]; hface < grid_.cell_facepos[cell+1]; ++hface) {
            const int face = grid_.cell_faces[hface];
            const double* fc = grid_.face_centroids + face*dim;
            double face_flux = 0.0;
            if (cell == grid_.face_cells[2*face]) {
                face_flux = flux_[face];
            } else {
                ASSERT(cell == grid_.face_cells[2*face + 1]);
                face_flux = -flux_[face];
            }
            for (int dd = 0; dd < dim; ++dd) {
                v[dd] += face_flux * (fc[dd] - cc[dd]) / grid_.cell_volumes[cell];
            }
        }
    }


    // --------  Methods of class VelocityInterpolationECVI  --------



    /// Constructor.
    /// \param[in]  grid   A grid.
    VelocityInterpolationECVI::VelocityInterpolationECVI(const UnstructuredGrid& grid)
        : bcmethod_(grid), grid_(grid)
    {
    }

    /// Set up fluxes for interpolation.
    /// Computes the corner velocities.
    /// \param[in]  flux   One signed flux per face in the grid.
    void VelocityInterpolationECVI::setupFluxes(const double* flux)
    {
        // We must now update the velocity member of the CornerInfo
        // for each corner.
        const int dim = grid_.dimensions;
        std::vector<double> N(dim*dim); // Normals matrix. Fortran ordering!
        std::vector<double> orig_N(dim*dim); // Normals matrix. Fortran ordering!
        std::vector<double> f(dim);     // Flux vector.
        std::vector<double> orig_f(dim);     // Flux vector.
        std::vector<MAT_SIZE_T> piv(dim); // For LAPACK solve
        const SparseTable<WachspressCoord::CornerInfo>& all_ci = bcmethod_.cornerInfo();
        const std::vector<int>& adj_faces = bcmethod_.adjacentFaces();
        corner_velocity_.resize(dim*all_ci.dataSize());
        const int num_cells = grid_.number_of_cells;
        for (int cell = 0; cell < num_cells; ++cell) {
            const int num_cell_corners = bcmethod_.numCorners(cell);
            for (int cell_corner = 0; cell_corner < num_cell_corners; ++cell_corner) {
                const int cid = all_ci[cell][cell_corner].corner_id;
                for (int adj_ix = 0; adj_ix < dim; ++adj_ix) {
                    const int face = adj_faces[dim*cid + adj_ix];
                    const double* fn = grid_.face_normals + dim*face;
                    for (int dd = 0; dd < dim; ++dd) {
                        N[adj_ix + dd*dim] = fn[dd]; // Row adj_ix, column dd
                    }
                    f[adj_ix] = flux[face];
                }
                // Now we have built N and f. Solve Nv = f.
                // Note that the face orientations do not matter,
                // as changing an orientation would negate both a
                // row in N and the corresponding element of f.
                // Solving linear equation with LAPACK.
                MAT_SIZE_T n = dim;
                MAT_SIZE_T nrhs = 1;
                MAT_SIZE_T lda = n;
                MAT_SIZE_T ldb = n;
                MAT_SIZE_T info = 0;
                orig_N = N;
                orig_f = f;
                dgesv_(&n, &nrhs, &N[0], &lda, &piv[0], &f[0], &ldb, &info);
                if (info != 0) {
                    // Print the local matrix and rhs.
                    std::cerr << "Failed solving single-cell system Nv = f in cell " << cell
                              << " with N = \n";
                    for (int row = 0; row < n; ++row) {
                        for (int col = 0; col < n; ++col) {
                            std::cerr << "    " << orig_N[row + n*col];
                        }
                        std::cerr << '\n';
                    }
                    std::cerr << "and f = \n";
                    for (int row = 0; row < n; ++row) {
                        std::cerr << "    " << orig_f[row] << '\n';
                    }
                    OPM_THROW(std::runtime_error, "Lapack error: " << info << " encountered in cell " << cell);
                }
                // The solution ends up in f, so we must copy it.
                std::copy(f.begin(), f.end(), corner_velocity_.begin() + dim*cid);
            }
        }
    }

    /// Interpolate velocity.
    /// \param[in]  cell   Cell in which to interpolate.
    /// \param[in]  x      Coordinates of point at which to interpolate.
    ///                    Must be array of length grid.dimensions.
    /// \param[out] v      Interpolated velocity.
    ///                    Must be array of length grid.dimensions.
    void VelocityInterpolationECVI::interpolate(const int cell,
                                                const double* x,
                                                double* v) const
    {
        const int n = bcmethod_.numCorners(cell);
        const int dim = grid_.dimensions;
        bary_coord_.resize(n);
        bcmethod_.cartToBary(cell, x, &bary_coord_[0]);
        std::fill(v, v + dim, 0.0);
        const SparseTable<WachspressCoord::CornerInfo>& all_ci = bcmethod_.cornerInfo();
        for (int i = 0; i < n; ++i) {
            const int cid = all_ci[cell][i].corner_id;
            for (int dd = 0; dd < dim; ++dd) {
                v[dd] += corner_velocity_[dim*cid + dd] * bary_coord_[i];
            }
        }
    }


} // namespace Opm
