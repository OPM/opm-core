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

#define BOOST_TEST_MODULE QuadratureTest
#include <boost/test/unit_test.hpp>

#include <opm/core/grid/CellQuadrature.hpp>
#include <opm/core/grid/FaceQuadrature.hpp>
#include <opm/core/grid/GridManager.hpp>
#include <opm/core/grid.h>
#include <cmath>

using namespace Opm;


namespace
{

    template <class Quadrature>
    class Integrator
    {
    public:
        explicit Integrator(const UnstructuredGrid& grid, const int entity, const int degree)
            : quad_(grid, entity, degree),
              pt(grid.dimensions)
        {
        }
        template <class Func>
        double integrate(const Func& f) const
        {
            double res = 0;
            for (int quad_pt = 0; quad_pt < quad_.numQuadPts(); ++quad_pt) {
                quad_.quadPtCoord(quad_pt, &pt[0]);
                const double w = quad_.quadPtWeight(quad_pt);
                const double fval = f(&pt[0]);
                res += w*fval;
            }
            return res;
        }
    private:
        const Quadrature quad_;
        mutable std::vector<double> pt;
    };


    template <class Quadrature, class Func>
    void testSingleCase(const UnstructuredGrid& grid,
                        const int entity,
                        const int degree,
                        const double expected_answer)
    {
        Integrator<Quadrature> integrator(grid, entity, degree);
        Func f;
        const double val = integrator.integrate(f);
        BOOST_CHECK(std::fabs(val - expected_answer) < 1e-12);
    }

} // anonymous namespace


namespace cart2d
{
    struct ConstantFunc
    {
        double operator()(const double*) const
        {
            return 1.234;
        }
    };

    struct LinearFunc
    {
        double operator()(const double* x) const
        {
            return 1.0*x[0] + 2.0*x[1] + 3.0;
        }
    };

    struct QuadraticFunc
    {
        double operator()(const double* x) const
        {
            return 3.0*x[0]*x[0] + 1.0*x[0]*x[1] + 2.0*x[1] + 3.0;
        }
    };

    static void test()
    {
        // Set up 2d 1-cell cartesian case.
        GridManager g(1, 1);
        const UnstructuredGrid& grid = *g.c_grid();

        // CellQuadrature tests.
        testSingleCase<CellQuadrature, ConstantFunc>(grid, 0, 1, 1.234);
        testSingleCase<CellQuadrature, LinearFunc>(grid, 0, 1, 4.5);
        testSingleCase<CellQuadrature, ConstantFunc>(grid, 0, 2, 1.234);
        testSingleCase<CellQuadrature, LinearFunc>(grid, 0, 2, 4.5);
        testSingleCase<CellQuadrature, QuadraticFunc>(grid, 0, 2, 5.25);

        // FaceQuadrature tests, degree 1 precision.
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 0, 1, 4);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 1, 1, 5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 2, 1, 3.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 3, 1, 5.5);

        // FaceQuadrature tests, degree 2 precision.
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 0, 2, 4);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 1, 2, 5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 2, 2, 3.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 3, 2, 5.5);

        // FaceQuadrature tests, quadratic function, degree 2 precision.
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 0, 2, 4.0);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 1, 2, 7.5);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 2, 2, 4.0);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 3, 2, 6.5);
    }
} // namespace cart2d


namespace cart3d
{
    struct LinearFunc
    {
        double operator()(const double* x) const
        {
            return 1.0*x[0] + 2.0*x[1] + 1.0*x[2] + 3.0;
        }
    };

    struct QuadraticFunc
    {
        double operator()(const double* x) const
        {
            return 1.0*x[0]*x[1] + 2.0*x[1] + 4.0*x[2] + 3.0;
        }
    };

    static void test()
    {
        // Set up 3d 1-cell cartesian case.
        GridManager g(1, 1, 1);
        const UnstructuredGrid& grid = *g.c_grid();

        // CellQuadrature tests.
        testSingleCase<CellQuadrature, LinearFunc>(grid, 0, 1, 5.0);
        testSingleCase<CellQuadrature, LinearFunc>(grid, 0, 2, 5.0);
        testSingleCase<CellQuadrature, QuadraticFunc>(grid, 0, 2, 6.25);

        // FaceQuadrature tests, degree 1 precision.
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 0, 1, 4.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 1, 1, 5.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 2, 1, 4.0);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 3, 1, 6.0);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 4, 1, 4.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 5, 1, 5.5);

        // FaceQuadrature tests, degree 2 precision.
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 0, 2, 4.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 1, 2, 5.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 2, 2, 4.0);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 3, 2, 6.0);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 4, 2, 4.5);
        testSingleCase<FaceQuadrature, LinearFunc>(grid, 5, 2, 5.5);

        // FaceQuadrature tests, quadratic function, degree 2 precision.
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 0, 2, 6.0);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 1, 2, 6.5);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 2, 2, 5.0);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 3, 2, 7.5);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 4, 2, 4.25);
        testSingleCase<FaceQuadrature, QuadraticFunc>(grid, 5, 2, 8.25);
    }

} // namespace cart3d

BOOST_AUTO_TEST_CASE(test_quadratures)
{
    cart2d::test();
    cart3d::test();
}
