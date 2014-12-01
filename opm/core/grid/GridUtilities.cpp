/*
  Copyright 2014 SINTEF ICT, Applied Mathematics.

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

#include <opm/core/grid/GridUtilities.hpp>
#include <opm/core/grid/GridHelpers.hpp>
#include <boost/math/constants/constants.hpp>
#include <set>
#include <vector>
#include <cmath>
#include <algorithm>

namespace Opm
{
    /// For each cell, find indices of all other cells sharing a vertex with it.
    /// \param[in] grid    A grid object.
    /// \return            A table of neighbour cell-indices by cell.
    SparseTable<int> cellNeighboursAcrossVertices(const UnstructuredGrid& grid)
    {
        // 1. Create vertex->cell mapping. We do this by iterating
        //    over all faces, and adding both its cell neighbours
        //    to each of its vertices' data.
        using namespace UgGridHelpers;
        const int num_vertices = grid.number_of_nodes;
        std::vector<std::set<int>> v2c(num_vertices);
        const int num_faces = numFaces(grid);
        const auto fc = faceCells(grid);
        for (int face = 0; face < num_faces; ++face) {
            for (int nodepos = grid.face_nodepos[face]; nodepos < grid.face_nodepos[face + 1]; ++nodepos) {
                const int vertex = grid.face_nodes[nodepos];
                for (int face_nb = 0; face_nb < 2; ++face_nb) {
                    const int face_nb_cell = fc(face, face_nb);
                    if (face_nb_cell >= 0) {
                        v2c[vertex].insert(face_nb_cell);
                    }
                }
            }
        }

        // 2. For each cell, iterate over its faces, iterate over
        //    their vertices, and collect all those vertices' cell
        //    neighbours. Add as row to sparse table.
        SparseTable<int> cell_nb;
        const int num_cells = numCells(grid);
        const auto c2f = cell2Faces(grid);
        // Reserve sufficient room for cartesian grids in 2 and 3
        // dimensions. Note that this is not a limit, just an
        // optimization similar to std::vector.
        cell_nb.reserve(num_cells, (dimensions(grid) == 2 ? 8 : 26) * num_cells);
        std::set<int> nb;
        for (int cell = 0; cell < num_cells; ++cell) {
            nb.clear();
            const auto cell_faces = c2f[cell];
            const int num_cell_faces = cell_faces.size();
            for (int local_face = 0; local_face < num_cell_faces; ++local_face) {
                const int face = cell_faces[local_face];
                for (int nodepos = grid.face_nodepos[face]; nodepos < grid.face_nodepos[face + 1]; ++nodepos) {
                    const int vertex = grid.face_nodes[nodepos];
                    nb.insert(v2c[vertex].begin(), v2c[vertex].end());
                }
            }
            nb.erase(cell);
            cell_nb.appendRow(nb.begin(), nb.end());
        }

        // 3. Done. Return.
        return cell_nb;
    }






    /// For each cell, order the (cell) neighbours counterclockwise.
    /// \param[in] grid    A 2d grid object.
    /// \param[in, out] nb A cell-cell neighbourhood table, such as from cellNeighboursAcrossVertices().
    void orderCounterClockwise(const UnstructuredGrid& grid,
                               SparseTable<int>& nb)
    {
        if (grid.dimensions != 2) {
            OPM_THROW(std::logic_error, "Cannot use orderCounterClockwise in " << grid.dimensions << " dimensions.");
        }
        const int num_cells = grid.number_of_cells;
        if (nb.size() != num_cells) {
            OPM_THROW(std::logic_error, "Inconsistent arguments for orderCounterClockwise().");
        }

        // For each cell, compute each neighbour's angle with the x axis,
        // sort that to find the correct permutation of the neighbours.
        typedef std::pair<double, int> AngleAndPos;
        std::vector<AngleAndPos> angle_and_pos;
        std::vector<int> original;
        for (int cell = 0; cell < num_cells; ++cell) {
            const int num_nb = nb[cell].size();
            angle_and_pos.clear();
            angle_and_pos.resize(num_nb);
            for (int ii = 0; ii < num_nb; ++ii) {
                const int cell2 = nb[cell][ii];
                const double v[2] = { grid.cell_centroids[2*cell2] - grid.cell_centroids[2*cell],
                                      grid.cell_centroids[2*cell2 + 1] - grid.cell_centroids[2*cell + 1] };
                // The formula below gives an angle in [0, 2*pi] with the positive x axis.
                const double angle = boost::math::double_constants::pi - std::atan2(v[1], -v[0]);
                angle_and_pos[ii] = std::make_pair(angle, ii);
            }
            original.assign(nb[cell].begin(), nb[cell].end());
            std::sort(angle_and_pos.begin(), angle_and_pos.end());
            for (int ii = 0; ii < num_nb; ++ii) {
                nb[cell][ii] = original[angle_and_pos[ii].second];
            }
        }
    }

} // namespace Opm
