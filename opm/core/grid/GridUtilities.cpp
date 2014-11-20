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
#include <set>
#include <vector>

namespace Opm
{
    /// For each cell, find indices of all other cells sharing a vertex with it.
    /// \param[in] grid    A grid object.
    /// \return            A table of neighbour cell-indices by cell.
    SparseTable<int> vertexNeighbours(const UnstructuredGrid& grid)
    {
	// 1. Create vertex->cell mapping. We do this by iterating
	//    over all faces, and adding both its cell neighbours
	//    to each of its vertices' data.
	const int num_vertices = grid.number_of_nodes;
	std::vector<std::set<int>> v2c(num_vertices);
	const int num_faces = grid.number_of_faces;
	for (int face = 0; face < num_faces; ++face) {
	    for (int nodepos = grid.face_nodepos[face]; nodepos < grid.face_nodepos[face + 1]; ++nodepos) {
		const int vertex = grid.face_nodes[nodepos];
		for (int face_nb = 0; face_nb < 2; ++face_nb) {
		    const int face_nb_cell = grid.face_cells[2*face + face_nb];
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
	const int num_cells = grid.number_of_cells;
	// Reserve sufficient room for cartesian grids in 2 and 3
	// dimensions. Note that this is not a limit, just an
	// optimization similar to std::vector.
	cell_nb.reserve(num_cells, (grid.dimensions == 2 ? 8 : 26) * num_cells);
	std::set<int> nb;
	for (int cell = 0; cell < num_cells; ++cell) {
	    nb.clear();
	    for (int hf = grid.cell_facepos[cell]; hf < grid.cell_facepos[cell + 1]; ++hf) {
		const int face = grid.cell_faces[hf];
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
} // namespace Opm
