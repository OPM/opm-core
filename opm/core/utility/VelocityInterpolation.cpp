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

#include <opm/core/utility/VelocityInterpolation.hpp>
#include <opm/core/grid.h>
#include <opm/core/linalg/blas_lapack.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <set>

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
                face_flux = -flux_[face];
            }
            for (int dd = 0; dd < dim; ++dd) {
                v[dd] += face_flux * (fc[dd] - cc[dd]) / grid_.cell_volumes[cell];
            }
        }
    }


    // -------- Helper methods for class VelocityInterpolationECVI --------

    namespace
    {
        /// Calculates the determinant of a 3 x 3 matrix, represented as
        /// three three-dimensional arrays.
        double determinantOf(const double* a0,
                             const double* a1,
                             const double* a2)
        {
            return
                a0[0] * (a1[1] * a2[2] - a2[1] * a1[2]) -
                a0[1] * (a1[0] * a2[2] - a2[0] * a1[2]) +
                a0[2] * (a1[0] * a2[1] - a2[0] * a1[1]);
        }
    } // anonymous namespace



    // --------  Methods of class VelocityInterpolationECVI  --------



    /// Constructor.
    /// \param[in]  grid   A grid.
    VelocityInterpolationECVI::VelocityInterpolationECVI(const UnstructuredGrid& grid)
        : grid_(grid)
    {
        const int dim = grid.dimensions;
        if (dim > Maxdim) {
            THROW("Grid has more than " << Maxdim << " dimensions.");
        }
        // Compute static data for each corner.
        const int num_cells = grid.number_of_cells;
        int corner_id_count = 0;
        for (int cell = 0; cell < num_cells; ++cell) {
            std::set<int> cell_vertices;
            std::vector<int> cell_faces;
            std::multimap<int, int> vertex_adj_faces;
            for (int hface = grid.cell_facepos[cell]; hface < grid.cell_facepos[cell + 1]; ++hface) {
                const int face = grid.cell_faces[hface];
                cell_faces.push_back(face);
                const int fn0 = grid.face_nodepos[face];
                const int fn1 = grid.face_nodepos[face + 1];
                cell_vertices.insert(grid.face_nodes + fn0, grid.face_nodes + fn1);
                for (int fn = fn0; fn < fn1; ++fn) {
                    const int vertex = grid.face_nodes[fn];
                    vertex_adj_faces.insert(std::make_pair(vertex, face));
                }
            }
            std::sort(cell_faces.begin(), cell_faces.end()); // set_difference requires sorted ranges
            std::vector<CornerInfo> cell_corner_info;
            std::set<int>::const_iterator it = cell_vertices.begin();
            for (; it != cell_vertices.end(); ++it) {
                CornerInfo ci;
                ci.corner_id = corner_id_count++;;
                ci.vertex = *it;
                typedef double* DblPtr;
                DblPtr fnorm[3] = { 0, 0, 0 };
                typedef std::multimap<int, int>::const_iterator MMIt;
                std::pair<MMIt, MMIt> frange = vertex_adj_faces.equal_range(ci.vertex);
                int fi = 0;
                std::vector<int> vert_adj_faces(dim);
                for (MMIt face_it = frange.first; face_it != frange.second; ++face_it, ++fi) {
                    if (fi >= dim) {
                        THROW("In cell " << cell << ", vertex " << ci.vertex << " has "
                              << " more than " << dim << " adjacent faces.");
                    }
                    fnorm[fi] = grid_.face_normals + dim*(face_it->second);
                    vert_adj_faces[fi] = face_it->second;
                }
                ASSERT(fi == dim);
                adj_faces_.insert(adj_faces_.end(), vert_adj_faces.begin(), vert_adj_faces.end());
                const double corner_vol = std::fabs(determinantOf(fnorm[0], fnorm[1], fnorm[2]));
                ci.volume = corner_vol;
                cell_corner_info.push_back(ci);
                std::sort(vert_adj_faces.begin(), vert_adj_faces.end());
                std::vector<int> vert_nonadj_faces(cell_faces.size() - vert_adj_faces.size());
                std::set_difference(cell_faces.begin(), cell_faces.end(),
                                    vert_adj_faces.begin(), vert_adj_faces.end(),
                                    vert_nonadj_faces.begin());
                nonadj_faces_.appendRow(vert_nonadj_faces.begin(), vert_nonadj_faces.end());
            }
            corner_info_.appendRow(cell_corner_info.begin(), cell_corner_info.end());
        }
        ASSERT(corner_id_count == corner_info_.dataSize());
    }

    /// Set up fluxes for interpolation.
    /// \param[in]  flux   One signed flux per face in the grid.
    void VelocityInterpolationECVI::setupFluxes(const double* flux)
    {
        // We must now update the velocity member of the CornerInfo
        // for each corner.
        const int dim = grid_.dimensions;
        std::vector<double> N(dim*dim); // Normals matrix. Fortran ordering!
        std::vector<double> f(dim);     // Flux vector.
        std::vector<MAT_SIZE_T> piv(dim); // For LAPACK solve
        const int num_cells = grid_.number_of_cells;
        for (int cell = 0; cell < num_cells; ++cell) {
            const int num_cell_corners = corner_info_[cell].size();
            for (int cell_corner = 0; cell_corner < num_cell_corners; ++cell_corner) {
                CornerInfo& ci = corner_info_[cell][cell_corner];
                for (int adj_ix = 0; adj_ix < dim; ++adj_ix) {
                    const int face = adj_faces_[dim*ci.corner_id + adj_ix];
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
                dgesv_(&n, &nrhs, &N[0], &lda, &piv[0], &f[0], &ldb, &info);
                if (info != 0) {
                    // Print the local matrix and rhs.
                    std::cerr << "Failed solving single-cell system Nv = f in cell " << cell
                              << " with N = \n";
                    for (int row = 0; row < n; ++row) {
                        for (int col = 0; col < n; ++col) {
                            std::cerr << "    " << N[row + n*col];
                        }
                        std::cerr << '\n';
                    }
                    std::cerr << "and f = \n";
                    for (int row = 0; row < n; ++row) {
                        std::cerr << "    " << f[row] << '\n';
                    }
                    THROW("Lapack error: " << info << " encountered in cell " << cell);
                }
                // The solution ends up in f, so we must copy it.
                std::copy(f.begin(), f.end(), ci.velocity);
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
        const int n = corner_info_[cell].size();
        const int dim = grid_.dimensions;
        bary_coord_.resize(n);
        cartToBaryWachspress(cell, x, &bary_coord_[0]);
        std::fill(v, v + dim, 0.0);
        for (int i = 0; i < n; ++i) {
            const CornerInfo& ci = corner_info_[cell][i];
            for (int dd = 0; dd < dim; ++dd) {
                v[dd] += ci.velocity[dd] * bary_coord_[i];
            }
        }
    }

    // Compute generalized barycentric coordinates for some point x
    // with respect to the vertices of a cell.
    void VelocityInterpolationECVI::cartToBaryWachspress(const int cell,
                                                         const double* x,
                                                         double* xb) const
    {
        const int n = corner_info_[cell].size();
        const int dim = grid_.dimensions;
        double totw = 0.0;
        for (int i = 0; i < n; ++i) {
            const CornerInfo& ci = corner_info_[cell][i];
            // Weight (unnormalized) is equal to:
            // V_i * (prod_{j \in nonadjacent faces} n_j * (c_j - x) )
            // ^^^                                   ^^^    ^^^
            // corner "volume"                    normal    centroid
            xb[i] = ci.volume;
            const int num_nonadj_faces = nonadj_faces_[ci.corner_id].size();
            for (int j = 0; j < num_nonadj_faces; ++j) {
                const int face = nonadj_faces_[ci.corner_id][j];
                double factor = 0.0;
                for (int dd = 0; dd < dim; ++dd) {
                    factor += grid_.face_normals[dim*face + dd]*(grid_.face_centroids[dim*face + dd] - x[dd]);
                }
                // Assumes outward-pointing normals, so negate factor if necessary.
                if (grid_.face_cells[2*face] != cell) {
                    ASSERT(grid_.face_cells[2*face + 1] == cell);
                    factor = -factor;
                }
                xb[i] *= factor;
            }
            totw += xb[i];
        }
        for (int i = 0; i < n; ++i) {
            xb[i] /= totw;
        }
    }


} // namespace Opm
