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
#include <algorithm>
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




    // --------  Methods of class VelocityInterpolationECVI  --------


    /// Constructor.
    /// \param[in]  grid   A grid.
    VelocityInterpolationECVI::VelocityInterpolationECVI(const UnstructuredGrid& grid)
        : grid_(grid)
    {
        if (grid.dimensions > Maxdim) {
            THROW("Grid has more than " << Maxdim << " dimensions.");
        }
        // Compute static data for each corner.
        const int num_cells = grid.number_of_cells;
        int corner_id_count = 0;
        for (int cell = 0; cell < num_cells; ++cell) {
            std::set<int> cell_vertices;
            std::multimap<int, int> vertex_adj_faces;
            for (int hface = grid.cell_facepos[cell]; hface < grid.cell_facepos[cell + 1]; ++hface) {
                const int face = grid.cell_faces[hface];
                const int fn0 = grid.face_nodepos[face];
                const int fn1 = grid.face_nodepos[face + 1];
                cell_vertices.insert(grid.face_nodes + fn0, grid.face_nodes + fn1);
                for (int fn = fn0; fn < fn1; ++fn) {
                    const int vertex = grid.face_nodes[fn];
                    vertex_adj_faces.insert(std::make_pair(vertex, face));
                }
            }
            std::vector<CornerInfo> cell_corner_info;
            std::set<int>::const_iterator it = cell_vertices.begin();
            for (; it != cell_vertices.end(); ++it) {
                CornerInfo ci;
                ci.corner_id = corner_id_count++;;
                ci.vertex = *it;
                cell_corner_info.push_back(ci);
                // nonadj_faces_.appendRow();
            }
            corner_info_.appendRow(cell_corner_info.begin(), cell_corner_info.end());
        }
        ASSERT(corner_id_count == corner_info_.dataSize());
    }

    /// Set up fluxes for interpolation.
    /// \param[in]  flux   One signed flux per face in the grid.
    void VelocityInterpolationECVI::setupFluxes(const double* flux)
    {
        THROW("Not implemented");
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
            const int num_nonadj_faces = nonadj_faces[ci.corner_id].size();
            for (int j = 0; j < num_nonadj_faces; ++j) {
                const int face = nonadj_faces[ci.corner_id][j];
                double factor = 0.0;
                for (int dd = 0; dd < dim; ++dd) {
                    factor += grid_.face_normals[dim*face + dd]*(grid_.face_centroids[dim*face + dd] - x[dd]);
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
