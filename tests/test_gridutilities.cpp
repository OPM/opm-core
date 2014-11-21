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


#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE GridUtilitiesTest
#include <boost/test/unit_test.hpp>

#include <opm/core/grid/GridUtilities.hpp>
#include <opm/core/grid/GridManager.hpp>

using namespace Opm;

BOOST_AUTO_TEST_CASE(cartesian_2d_vertexNeighbours)
{
    const GridManager gm(2, 2);
    const UnstructuredGrid& grid = *gm.c_grid();
    const SparseTable<int> vnb = vertexNeighbours(grid);

    const int num_elem = 12;
    const int elem[num_elem] = { 1, 2, 3, 0, 2, 3, 0, 1, 3, 0, 1, 2 };
    const int num_rows = 4;
    const int rowsizes[num_rows] = { 3, 3, 3, 3 };
    const SparseTable<int> truth(elem, elem + num_elem, rowsizes, rowsizes + num_rows);
    BOOST_CHECK(vnb == truth);
}

BOOST_AUTO_TEST_CASE(cartesian_3d_vertexNeighbours)
{
    const GridManager gm(3, 2, 2);
    const UnstructuredGrid& grid = *gm.c_grid();
    const SparseTable<int> vnb = vertexNeighbours(grid);

    BOOST_CHECK_EQUAL(int(vnb.size()), grid.number_of_cells);
    BOOST_REQUIRE(!vnb.empty());
    const int n = 7;
    BOOST_CHECK_EQUAL(int(vnb[0].size()), n);
    const int nb[n] = { 1, 3, 4, 6, 7, 9, 10 };
    BOOST_CHECK_EQUAL_COLLECTIONS(vnb[0].begin(), vnb[0].end(), nb, nb + n);
}

BOOST_AUTO_TEST_CASE(cartesian_2d_orderCounterClockwise)
{
    const GridManager gm(2, 2);
    const UnstructuredGrid& grid = *gm.c_grid();
    SparseTable<int> vnb = vertexNeighbours(grid);
    orderCounterClockwise(grid, vnb);

    BOOST_REQUIRE(!vnb.empty());
    const int num_elem = 12;
    const int elem[num_elem] = { 1, 3, 2, 3, 2, 0, 3, 0, 1, 2, 0, 1 };
    const int num_rows = 4;
    const int rowsizes[num_rows] = { 3, 3, 3, 3 };
    const SparseTable<int> truth(elem, elem + num_elem, rowsizes, rowsizes + num_rows);
    BOOST_CHECK(vnb == truth);
    for (int c = 0; c < num_rows; ++c) {
	BOOST_CHECK_EQUAL_COLLECTIONS(vnb[c].begin(), vnb[c].end(), truth[c].begin(), truth[c].end());
    }
}
