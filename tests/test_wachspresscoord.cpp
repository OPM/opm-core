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

#include <config.h>

#if HAVE_DYNAMIC_BOOST_TEST
#define BOOST_TEST_DYN_LINK
#endif
#define NVERBOSE // to suppress our messages when throwing

#define BOOST_TEST_MODULE WachspressCoordTest
#include <boost/test/unit_test.hpp>

#include <opm/core/utility/WachspressCoord.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid.h>
#include <cmath>

using namespace Opm;


namespace
{
    class Interpolator
    {
    public:
        explicit Interpolator(const UnstructuredGrid& grid)
        : bcmethod_(grid), grid_(grid)
        {
        }
        template <class Func>
        double interpolate(const Func& f,
                           const int cell,
                           const std::vector<double>& x) const
        {
            const int ncor = bcmethod_.numCorners(cell);
            bary_coord_.resize(ncor);
            bcmethod_.cartToBary(cell, &x[0], &bary_coord_[0]);
            double val = 0.0;
            for (int cor = 0; cor < ncor; ++cor) {
                const int vertex = bcmethod_.cornerInfo()[cell][cor].vertex;
                const double vval = f(grid_.node_coordinates + grid_.dimensions*vertex);
                val += vval*bary_coord_[cor];
            }
            return val;
        }
    private:
        WachspressCoord bcmethod_;
        const UnstructuredGrid grid_;
        mutable std::vector<double> bary_coord_;
    };

} // anonymous namespace

struct LinearFunc
{
    double operator()(const double* x) const
    {
        return 1.0*x[0] + 2.0*x[1] + 3.0;
    }
};

static void test2dCart()
{
    // Set up 2d 1-cell cartesian case.
    GridManager g(1, 1);
    const UnstructuredGrid& grid = *g.c_grid();
    Interpolator interp(grid);
    LinearFunc f;

    // Test a few points
    std::vector<double> x(2);
    x[0] = 0.23456;
    x[1] = 0.87654;
    double val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.5;
    x[1] = 0.5;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 1.0;
    x[1] = 0.5;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 1.0;
    x[1] = 1.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
}


namespace
{

    // Data for a pyramid. Node numbering goes
    // lexicographic on bottom, then top.
    // Face numbering goes xmin, xmax, ymin, ymax, bottom.
    namespace Pyramid
    {
        static int face_nodes[]   = { 0, 4, 2,    3, 4, 1,    0, 1, 4,    4, 3, 2,    0, 2, 3, 1,       };
        static int face_nodepos[] = { 0,          3,          6,          9,          12,            16 };
        static int face_cells[]   = { 0, -1,      0, -1,      0, -1,      0, -1,      0, -1             };
        static int cell_faces[]   = { 0, 1, 2, 3, 4 };
        static int cell_facepos[] = { 0, 5 };
        static double node_coordinates[] = { 0.0, 0.0, 0.0,   1.0, 0.0, 0.0,   0.0, 1.0, 0.0,   1.0, 1.0, 0.0,   0.0, 0.0, 1.0 };
        static double face_centroids[]   = { 0,       1.0/3.0, 1.0/3.0,
                                             2.0/3.0, 1.0/3.0, 1.0/3.0,
                                             1.0/3.0, 0,       1.0/3.0,
                                             1.0/3.0, 2.0/3.0, 1.0/3.0,
                                             0.5,     0.5,     0        };
        static double face_areas[] = { 0.5, std::sqrt(2.0), 0.5, std::sqrt(2.0), 1.0 };
        static double face_normals[] = {   -0.5000,         0,         0,
                                            0.5000,         0,    0.5000,
                                            0,        -0.5000,         0,
                                            0,         0.5000,    0.5000,
                                            0,              0,   -1.0000    };
        static double cell_centroids[] = { 0.375, 0.375, 0.25 };
        static double cell_volumes[] = { 1.0/3.0 };

    } // namespace Pyramid

    UnstructuredGrid makePyramid()
    {
        // Make a 3d 1-cell grid, where the single cell is a pyramid.
        UnstructuredGrid grid;
        grid.dimensions = 3;
        grid.number_of_cells = 1;
        grid.number_of_faces = 5;
        grid.number_of_nodes = 5;
        grid.face_nodes = Pyramid::face_nodes;
        grid.face_nodepos = Pyramid::face_nodepos;
        grid.face_cells = Pyramid::face_cells;
        grid.cell_faces = Pyramid::cell_faces;
        grid.cell_facepos = Pyramid::cell_facepos;
        grid.node_coordinates = Pyramid::node_coordinates;
        grid.face_centroids = Pyramid::face_centroids;
        grid.face_areas = Pyramid::face_areas;
        grid.face_normals = Pyramid::face_normals;
        grid.cell_centroids = Pyramid::cell_centroids;
        grid.cell_volumes = Pyramid::cell_volumes;
        return grid;
    }

} // anonymous namespace


static void testPyramid()
{
    // Set up a 3d 1-cell non-cartesian case (a pyramid).
    UnstructuredGrid grid = makePyramid();
    Interpolator interp(grid);
    LinearFunc f;

    // Test a few points
    std::vector<double> x(3);
    x[0] = 0.123;
    x[1] = 0.0123;
    x[2] = 0.213;
    double val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 1.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 0.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 0.1;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
}


namespace
{
    // Data for an irregular 2d polygon.
    namespace Irreg2d
    {
        static int face_nodes[]   = { 0, 1,    1, 2,    2, 3,    3, 4,    4, 0        };
        static int face_nodepos[] = { 0,       2,       4,       6,       8,       10 };
        static int face_cells[]   = { 0, -1,   0, -1,   0, -1,   0, -1,   0, -1       };
        static int cell_faces[]   = { 0, 1, 2, 3, 4 };
        static int cell_facepos[] = { 0, 5 };
        static double node_coordinates[] = { 0, 0,    3, 0,    3, 2,    1, 3,    0, 2 };
        static double face_centroids[]   = { 1.5, 0,    3, 1,    2, 2.5,    0.5, 2.5,    0, 1 };
        static double face_areas[] = { 3, 2, std::sqrt(5.0), std::sqrt(2.0), 2 };
        static double face_normals[] = { 0, -3,    2, 0,    1, 2,    -1, 1,    -2, 0 };
        static double cell_centroids[] = {  22.0/15.0, 19.0/15.0 };
        static double cell_volumes[] = { 7.5 };
    } // namespace Irreg2d

    UnstructuredGrid makeIrreg2d()
    {
        // Make a 2d 1-cell grid, where the single cell is a polygon.
        UnstructuredGrid grid;
        grid.dimensions = 2;
        grid.number_of_cells = 1;
        grid.number_of_faces = 5;
        grid.number_of_nodes = 5;
        grid.face_nodes = Irreg2d::face_nodes;
        grid.face_nodepos = Irreg2d::face_nodepos;
        grid.face_cells = Irreg2d::face_cells;
        grid.cell_faces = Irreg2d::cell_faces;
        grid.cell_facepos = Irreg2d::cell_facepos;
        grid.node_coordinates = Irreg2d::node_coordinates;
        grid.face_centroids = Irreg2d::face_centroids;
        grid.face_areas = Irreg2d::face_areas;
        grid.face_normals = Irreg2d::face_normals;
        grid.cell_centroids = Irreg2d::cell_centroids;
        grid.cell_volumes = Irreg2d::cell_volumes;
        return grid;
    }

} // anonymous namespace


static void testIrreg2d()
{
    // Set up a 2d 1-cell where the single cell is a polygon.
    UnstructuredGrid grid = makeIrreg2d();
    Interpolator interp(grid);
    LinearFunc f;

    // Test a few points
    std::vector<double> x(2);
    x[0] = 1.2345;
    x[1] = 2.0123;
    double val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.0;
    x[1] = 0.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 1.0;
    x[1] = 3.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 3.0;
    x[1] = 1.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
}


namespace
{
    // Data for an irregular 3d prism.
    namespace IrregPrism
    {
        static int face_nodes[]   = { 0, 4, 2, 1, 3, 5, 0, 1, 5, 4, 2, 4, 5, 3, 2, 3, 0, 1};
        static int face_nodepos[] = { 0, 3, 6, 10, 14, 18 };
        static int face_cells[]   = { 0, -1,   0, -1,   0, -1,   0, -1,   0, -1 };
        static int cell_faces[]   = { 0, 1, 2, 3, 4 };
        static int cell_facepos[] = { 0, 5 };
        static double node_coordinates[] = { 0, 0, 0,
                                             2, 0, 0,
                                             0, 1, 0,
                                             2, 1, 0,
                                             0, 0, 1,
                                             1, 0, 1 };
        static double face_centroids[]   = { 0,       1.0/3.0, 1.0/3.0,
                                             5.0/3.0, 1.0/3.0, 1.0/3.0,
                                             7.0/9.0, 0,       4.0/9.0,
                                             7.0/9.0, 5.0/9.0, 4.0/9.0,
                                             1,       0.5,     0        };
        static double face_areas[] = {    0.500000000000000,
                                          0.707106781186548,
                                          1.500000000000000,
                                          2.121320343559642,
                                          2.000000000000000 };
        static double face_normals[] = {   -0.500000000000000,                   0,                   0,
                                            0.500000000000000,   0.000000000000000,   0.500000000000000,
                                                            0,  -1.500000000000000,                   0,
                                                            0,   1.500000000000000,   1.500000000000000,
                                                            0,                   0,   -2.000000000000000        };
        static double cell_centroids[] = { 0.85, 0.35, 0.3 };
        static double cell_volumes[] = { 5.0/6.0 };
    } // namespace IrregPrism

    UnstructuredGrid makeIrregPrism()
    {
        // Make a 3d 1-cell grid, where the single cell is a prism.
        UnstructuredGrid grid;
        grid.dimensions = 3;
        grid.number_of_cells = 1;
        grid.number_of_faces = 5;
        grid.number_of_nodes = 6;
        grid.face_nodes = IrregPrism::face_nodes;
        grid.face_nodepos = IrregPrism::face_nodepos;
        grid.face_cells = IrregPrism::face_cells;
        grid.cell_faces = IrregPrism::cell_faces;
        grid.cell_facepos = IrregPrism::cell_facepos;
        grid.node_coordinates = IrregPrism::node_coordinates;
        grid.face_centroids = IrregPrism::face_centroids;
        grid.face_areas = IrregPrism::face_areas;
        grid.face_normals = IrregPrism::face_normals;
        grid.cell_centroids = IrregPrism::cell_centroids;
        grid.cell_volumes = IrregPrism::cell_volumes;
        return grid;
    }

} // anonymous namespace

static void testIrregPrism()
{
    // Set up a 3d 1-cell non-cartesian case (a prism).
    UnstructuredGrid grid = makeIrregPrism();
    Interpolator interp(grid);
    LinearFunc f;

    // Test a few points
    std::vector<double> x(3);
    x[0] = 0.123;
    x[1] = 0.0123;
    x[2] = 0.213;
    double val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 1.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 1.0;
    x[1] = 0.0;
    x[2] = 1.0;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 0.5;
    val = interp.interpolate(f, 0, x);
    BOOST_CHECK(std::fabs(val - f(&x[0])) < 1e-12);
}


BOOST_AUTO_TEST_CASE(test_WachspressCoord)
{
    test2dCart();
    BOOST_CHECK_THROW(testPyramid(), std::exception);
    testIrreg2d();
    testIrregPrism();
}
