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

#ifndef OPM_WACHSPRESSCOORD_HEADER_INCLUDED
#define OPM_WACHSPRESSCOORD_HEADER_INCLUDED

#include <opm/core/utility/SparseTable.hpp>
#include <vector>

struct UnstructuredGrid;

namespace Opm
{

    /// Class capable of computing Wachspress coordinates in 2d and 3d.
    /// The formula used is a modification of the formula given in:
    /// M. Meyer, A. Barr, H. Lee, and M. Desbrun.
    /// Generalized barycentric coordinates on irregular poly-
    /// gons. Journal of Graphics Tools, 7(1):13–22, 2002.
    class WachspressCoord
    {
    public:
        /// Constructor.
        /// \param[in]  grid   A grid.
        explicit WachspressCoord(const UnstructuredGrid& grid);

        /// Count of vertices adjacent to a call.
        /// \param[in]  cell   A cell index.
        /// \return            Number of corners of cell.
        int numCorners(const int cell) const;

        /// Compute generalized barycentric coordinates for some point x
        /// with respect to the vertices of a grid cell.
        /// \param[in]  cell   Cell in which to compute coordinates.
        /// \param[in]  x      Coordinates of point in cartesian coordinates.
        ///                    Must be array of length grid.dimensions.
        /// \param[out] xb     Coordinates of point in barycentric coordinates.
        ///                    Must be array of length numCorners(cell).
        void cartToBary(const int cell,
                        const double* x,
                        double* xb) const;

        // A corner is here defined as a {cell, vertex} pair where the
        // vertex is adjacent to the cell.
        struct CornerInfo
        {
            int corner_id;         // Unique for each corner.
            int vertex;            // Shared between corners belonging to different cells.
            double volume;         // Defined as det(N) where N is the matrix of adjacent face normals.
        };

        /// The class stores some info for each corner.
        /// \return            The corner info container.
        const SparseTable<CornerInfo>& cornerInfo() const;

    private:
        const UnstructuredGrid& grid_;
        SparseTable<CornerInfo> corner_info_;   // Corner info by cell.
        std::vector<int> adj_faces_;    // Set of adjacent faces, by corner id. Contains dim face indices per corner.
        SparseTable<int> nonadj_faces_; // Set of nonadjacent faces, by corner id.
    };

} // namespace Opm

#endif // OPM_WACHSPRESSCOORD_HEADER_INCLUDED
