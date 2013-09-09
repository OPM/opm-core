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
#include <opm/core/utility/WachspressCoord.hpp>
#include <opm/core/grid.h>
#include <cmath>
#include <map>
#include <set>

namespace Opm
{


    // -------- Helper methods for class WachspressCoord --------

    namespace
    {
        /// Calculates the determinant of a 2 x 2 matrix, represented as
        /// two two-dimensional arrays.
        double determinantOf(const double* a0,
                             const double* a1)
        {
            return
                a0[0] * a1[1] - a0[1] * a1[0];
        }

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

        /// Calculates the volume of the parallelepiped given by
        /// the vectors n[i] for i = 0..(dim-1), each n[i] is of size dim.
        double cornerVolume(double** n, const int dim)
        {
            assert(dim == 2 || dim == 3);
            double det = (dim == 2) ? determinantOf(n[0], n[1]) : determinantOf(n[0], n[1], n[2]);
            return std::fabs(det);
        }

    } // anonymous namespace



    // --------  Methods of class WachspressCoord  --------

    // The formula used is a modification of the formula given in:
    // M. Meyer, A. Barr, H. Lee, and M. Desbrun.
    // Generalized barycentric coordinates on irregular poly-
    // gons. Journal of Graphics Tools, 7(1):13–22, 2002.
    //
    // The formula given there is, for a corner i,
    //     b_i = w_i / sum_{k} w_k
    //     w_i = V_i / (prod_{j \in adjacent faces} n_j * (x_i - x) )
    //           ^^^                                ^^^    ^^^
    //           corner "volume"                 normal    corner coordinates
    //     V_i = |Det({n_j}_{j \in adjacent faces})|
    // The corner coordinate x_i above can be replaced with any point on face j
    // without changing the value of w_i, and we replace it with c_j, the face
    // centroid.
    // However, this formula has the problem that the denominator of w_i becomes zero
    // close to the boundary. Our solution is to multiply all w_i by
    /// prod_{all j} n_j * (c_j - x), resulting in the formula:
    //     w_i = V_i * (prod_{j \in nonadjacent faces} n_j * (c_j - x) ).
    // Another implementation note is that the above formulas assumes that
    // the normals have length 1 (are unit normals). It is easy to see that this
    // can be relaxed, since each normal occurs once in the formula for w_i, and
    // all w_i will be scaled by the same number. In our implementation we therefore
    // use the area-scaled normals directly as provided by the UnstructuredGrid.

    /// Constructor.
    /// \param[in]  grid   A grid.
    WachspressCoord::WachspressCoord(const UnstructuredGrid& grid)
        : grid_(grid)
    {
        enum { Maxdim = 3 };
        const int dim = grid.dimensions;
        if (dim > Maxdim) {
            OPM_THROW(std::runtime_error, "Grid has more than " << Maxdim << " dimensions.");
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
                double* fnorm[Maxdim] = { 0 };
                typedef std::multimap<int, int>::const_iterator MMIt;
                std::pair<MMIt, MMIt> frange = vertex_adj_faces.equal_range(ci.vertex);
                int fi = 0;
                std::vector<int> vert_adj_faces(dim);
                for (MMIt face_it = frange.first; face_it != frange.second; ++face_it, ++fi) {
                    if (fi >= dim) {
                        OPM_THROW(std::runtime_error, "In cell " << cell << ", vertex " << ci.vertex << " has "
                              << " more than " << dim << " adjacent faces.");
                    }
                    fnorm[fi] = grid_.face_normals + dim*(face_it->second);
                    vert_adj_faces[fi] = face_it->second;
                }
                assert(fi == dim);
                adj_faces_.insert(adj_faces_.end(), vert_adj_faces.begin(), vert_adj_faces.end());
                const double corner_vol = cornerVolume(fnorm, dim);
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
        assert(corner_id_count == corner_info_.dataSize());
    }




    /// Count of vertices adjacent to a call.
    /// \param[in]  cell   A cell index.
    /// \return            Number of corners of cell.
    int WachspressCoord::numCorners(const int cell) const
    {
        return corner_info_[cell].size();
    }


    /// The class stores some info for each corner.
    /// \return            The corner info container.
    const SparseTable<WachspressCoord::CornerInfo>& WachspressCoord::cornerInfo() const
    {
        return corner_info_;
    }



    /// The class stores some info for each corner.
    /// \return            The corner info container.
    const std::vector<int>& WachspressCoord::adjacentFaces() const
    {
        return adj_faces_;
    }



    /// Compute generalized barycentric coordinates for some point x
    /// with respect to the vertices of a grid cell.
    /// \param[in]  cell   Cell in which to compute coordinates.
    /// \param[in]  x      Coordinates of point in cartesian coordinates.
    ///                    Must be array of length grid.dimensions.
    /// \param[out] xb     Coordinates of point in barycentric coordinates.
    ///                    Must be array of length numCorners(cell).
    void WachspressCoord::cartToBary(const int cell,
                                     const double* x,
                                     double* xb) const
    {
        // Note:
        // A possible optimization is: compute all n_j * (c_j - x) factors
        // once, instead of repeating computation for all corners (for
        // which j is a nonadjacent face).
        const int n = numCorners(cell);
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
                    assert(grid_.face_cells[2*face + 1] == cell);
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
