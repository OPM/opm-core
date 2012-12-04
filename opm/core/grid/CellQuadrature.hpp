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

#ifndef OPM_CELLQUADRATURE_HEADER_INCLUDED
#define OPM_CELLQUADRATURE_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/utility/ErrorMacros.hpp>
#include <algorithm>
#include <cmath>

namespace Opm
{


    namespace {
        /// Calculates the determinant of a 3 x 3 matrix, represented as
        /// three three-dimensional arrays.
        inline double determinantOf(const double* a0,
                                    const double* a1,
                                    const double* a2)
        {
            return
                a0[0] * (a1[1] * a2[2] - a2[1] * a1[2]) -
                a0[1] * (a1[0] * a2[2] - a2[0] * a1[2]) +
                a0[2] * (a1[0] * a2[1] - a2[0] * a1[1]);
        }

        /// Computes the volume of a tetrahedron consisting of 4 vertices
        /// with 3-dimensional coordinates
        inline double tetVolume(const double* p0,
                                const double* p1,
                                const double* p2,
                                const double* p3)
        {
            double a[3] = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
            double b[3] = { p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2] };
            double c[3] = { p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2] };
            return std::fabs(determinantOf(a, b, c) / 6.0);
        }
    } // anonymous namespace


    /// A class providing numerical quadrature for cells.
    /// In general: \int_{cell} g(x) dx = \sum_{i=0}^{n-1} w_i g(x_i).
    /// Note that this class does multiply weights by cell volume,
    /// so weights always sum to cell volume.
    /// Degree 1 method:
    ///     Midpoint (centroid) method.
    ///         n = 1, w_0 = cell volume, x_0 = cell centroid
    /// Degree 2 method:
    ///    Based on subdivision of each cell face into triangles
    ///    with the face centroid as a common vertex, and then
    ///    subdividing the cell into tetrahedra with the cell
    ///    centroid as a common vertex. Then apply the tetrahedron
    ///    rule with the following 4 nodes (uniform weights):
    ///        a = 0.138196601125010515179541316563436
    ///        x_i has all barycentric coordinates = a, except for
    ///            the i'th coordinate which is = 1 - 3a.
    ///    This rule is from http://nines.cs.kuleuven.be/ecf,
    ///    it is the second degree 2 4-point rule for tets,
    ///    referenced to Stroud(1971).
    ///    The tetrahedra are numbered T_{i,j}, and are given by the
    ///    cell centroid C, the face centroid FC_i, and two nodes
    ///    of face i: FN_{i,j}, FN_{i,j+1}.
    class CellQuadrature
    {
    public:
        CellQuadrature(const UnstructuredGrid& grid,
                       const int cell,
                       const int degree)
            : grid_(grid), cell_(cell), degree_(degree)
        {
            if (grid.dimensions != 3) {
                THROW("CellQuadrature only implemented for 3D case.");
            }
            if (degree > 2) {
                THROW("CellQuadrature exact for polynomial degrees > 1 not implemented.");
            }
        }

        int numQuadPts() const
        {
            if (degree_ < 2) {
                return 1;
            }
            // Degree 2 case.
            int sumnodes = 0;
            for (int hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
                const int face = grid_.cell_faces[hf];
                sumnodes += grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
            }
            return 4*sumnodes;
        }

        void quadPtCoord(const int index, double* coord) const
        {
            const int dim = grid_.dimensions;
            const double* cc = grid_.cell_centroids + dim*cell_;
            if (degree_ < 2) {
                std::copy(cc, cc + dim, coord);
                return;
            }
            // Degree 2 case.
            int tetindex = index / 4;
            const int subindex = index % 4;
            const double* nc = grid_.node_coordinates;
            for (int hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
                const int face = grid_.cell_faces[hf];
                const int nfn = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
                if (nfn <= tetindex) {
                    // Our tet is not associated with this face.
                    tetindex -= nfn;
                    continue;
                }
                const double* fc = grid_.face_centroids + dim*face;
                const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face];
                const int node0 = fnodes[tetindex];
                const int node1 = fnodes[(tetindex + 1) % nfn];
                const double* n0c = nc + dim*node0;
                const double* n1c = nc + dim*node1;
                const double a = 0.138196601125010515179541316563436;
                // Barycentric coordinates of our point in the tet.
                double baryc[4] = { a, a, a, a };
                baryc[subindex] = 1.0 - 3.0*a;
                for (int dd = 0; dd < dim; ++dd) {
                    coord[dd] = baryc[0]*cc[dd] + baryc[1]*fc[dd] + baryc[2]*n0c[dd] + baryc[3]*n1c[dd];
                }
                return;
            }
            THROW("Should never reach this point.");
        }

        double quadPtWeight(const int index) const
        {
            if (degree_ < 2) {
                return grid_.cell_volumes[cell_];
            }
            // Degree 2 case.
            const int dim = grid_.dimensions;
            const double* cc = grid_.cell_centroids + dim*cell_;
            int tetindex = index / 4;
            const double* nc = grid_.node_coordinates;
            for (int hf = grid_.cell_facepos[cell_]; hf < grid_.cell_facepos[cell_ + 1]; ++hf) {
                const int face = grid_.cell_faces[hf];
                const int nfn = grid_.face_nodepos[face + 1] - grid_.face_nodepos[face];
                if (nfn <= tetindex) {
                    // Our tet is not associated with this face.
                    tetindex -= nfn;
                    continue;
                }
                const double* fc = grid_.face_centroids + dim*face;
                const int* fnodes = grid_.face_nodes + grid_.face_nodepos[face];
                const int node0 = fnodes[tetindex];
                const int node1 = fnodes[(tetindex + 1) % nfn];
                const double* n0c = nc + dim*node0;
                const double* n1c = nc + dim*node1;
                return 0.25*tetVolume(cc, fc, n0c, n1c);
            }
            THROW("Should never reach this point.");
        }

    private:
        const UnstructuredGrid& grid_;
        const int cell_;
        const int degree_;
    };

} // namespace Opm

#endif // OPM_CELLQUADRATURE_HEADER_INCLUDED
