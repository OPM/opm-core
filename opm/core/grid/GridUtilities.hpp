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

#ifndef OPM_GRIDUTILITIES_HEADER_INCLUDED
#define OPM_GRIDUTILITIES_HEADER_INCLUDED

#include <opm/core/grid.h>
#include <opm/core/utility/SparseTable.hpp>

namespace Opm
{

    /// For each cell, find indices of all cells sharing a vertex with it.
    /// \param[in] grid    A grid object.
    /// \return            A table of neighbour cell-indices by cell.
    SparseTable<int> vertexNeighbours(const UnstructuredGrid& grid);

} // namespace Opm

#endif // OPM_GRIDUTILITIES_HEADER_INCLUDED
