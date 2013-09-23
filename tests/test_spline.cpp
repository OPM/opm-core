// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010-2012 by Andreas Lauser                               *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief This is a program to test the polynomial spline interpolation.
 *
 * It just prints some function to stdout. You can look at the result
 * using the following commands:
 *
----------- snip -----------
./test_spline > spline.csv
gnuplot

gnuplot> plot "spline.csv" using 1:2 w l ti "Curve", \
              "spline.csv" using 1:3 w l ti "Derivative", \
              "spline.csv" using 1:4 w p ti "Monotonical"
----------- snap -----------


*/
#include "config.h"

#include <array>
#include <opm/core/utility/Spline.hpp>

#define GCC_VERSION (__GNUC__ * 10000 \
               + __GNUC_MINOR__ * 100 \
               + __GNUC_PATCHLEVEL__)

template <class Spline>
void testCommon(const Spline &sp,
                const double *x,
                const double *y)
{
    static double eps = 1e-10;

    int n = sp.numSamples();
    for (int i = 0; i < n; ++i) {
        // sure that we hit all sampling points
        double y0 = (i>0)?sp.eval(x[i]-eps):y[0];
        double y1 = sp.eval(x[i]);
        double y2 = (i<n-1)?sp.eval(x[i]+eps):y[n-1];

        if (std::abs(y0 - y[i]) > 100*eps || std::abs(y2 - y[i]) > 100*eps)
            OPM_THROW(std::runtime_error,
                       "Spline seems to be discontinuous at sampling point " << i << "!");
        if (std::abs(y1 - y[i]) > eps)
            OPM_THROW(std::runtime_error,
                       "Spline does not capture sampling point " << i << "!");

        // make sure the derivative is continuous (assuming that the
        // second derivative is smaller than 1000)
        double d1 = sp.evalDerivative(x[i]);
        double d0 = (i>0)?sp.evalDerivative(x[i]-eps):d1;
        double d2 = (i<n-1)?sp.evalDerivative(x[i]+eps):d1;

        if (std::abs(d1 - d0) > 1000*eps || std::abs(d2 - d0) > 1000*eps)
            OPM_THROW(std::runtime_error,
                       "Spline seems to exhibit a discontinuous derivative at sampling point " << i << "!");
    }
}

template <class Spline>
void testFull(const Spline &sp,
              const double *x,
              const double *y,
              double m0,
              double m1)
{
    // test the common properties of splines
    testCommon(sp, x, y);

    static double eps = 1e-5;
    int n = sp.numSamples();

    // make sure the derivative at both end points is correct
    double d0 = sp.evalDerivative(x[0]);
    double d1 = sp.evalDerivative(x[n-1]);
    if (std::abs(d0 - m0) > eps)
        OPM_THROW(std::runtime_error,
                   "Invalid derivative at beginning of interval: is "
                   << d0 << " ought to be " << m0);
    if (std::abs(d1 - m1) > eps)
        OPM_THROW(std::runtime_error,
                   "Invalid derivative at end of interval: is "
                   << d1 << " ought to be " << m1);
}

template <class Spline>
void testNatural(const Spline &sp,
                 const double *x,
                 const double *y)
{
    // test the common properties of splines
    testCommon(sp, x, y);

    static double eps = 1e-5;
    int n = sp.numSamples();

    // make sure the second derivatives at both end points are 0
    double d0 = sp.evalDerivative(x[0]);
    double d1 = sp.evalDerivative(x[0] + eps);

    double d2 = sp.evalDerivative(x[n-1] - eps);
    double d3 = sp.evalDerivative(x[n-1]);

    if (std::abs(d1 - d0)/eps > 1000*eps)
        OPM_THROW(std::runtime_error,
                   "Invalid second derivative at beginning of interval: is "
                   << (d1 - d0)/eps << " ought to be 0");

    if (std::abs(d3 - d2)/eps > 1000*eps)
        OPM_THROW(std::runtime_error,
                   "Invalid second derivative at end of interval: is "
                   << (d3 - d2)/eps << " ought to be 0");
}

void testAll()
{
    double x[] = { 0, 5, 7.5, 8.75, 9.375 };
    double y[] = { 10, 0, 10, 0, 10 };
    double m0 = 10;
    double m1 = -10;
    double points[][2] =
        {
            {x[0], y[0]},
            {x[1], y[1]},
            {x[2], y[2]},
            {x[3], y[3]},
            {x[4], y[4]},
        };


#if GCC_VERSION >= 40500
    std::initializer_list<const std::pair<double, double> > pointsInitList =
        {
            {x[0], y[0]},
            {x[1], y[1]},
            {x[2], y[2]},
            {x[3], y[3]},
            {x[4], y[4]},
        };
#endif

    std::vector<double> xVec;
    std::vector<double> yVec;
    std::vector<double*> pointVec;
    for (int i = 0; i < 5; ++i) {
        xVec.push_back(x[i]);
        yVec.push_back(y[i]);
        pointVec.push_back(points[i]);
    }

    /////////
    // test spline with two sampling points
    /////////

    // full spline
    { Opm::Spline<double> sp(x[0], x[1], y[0], y[1], m0, m1); sp.set(x[0],x[1],y[0],y[1],m0, m1); testFull(sp, x, y, m0, m1); };
    { Opm::Spline<double> sp(2, x, y, m0, m1); sp.setXYArrays(2, x, y, m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp(2, points, m0, m1); sp.setArrayOfPoints(2, points, m0, m1); testFull(sp, x, y, m0, m1); };

    /////////
    // test variable length splines
    /////////

    // full spline
    { Opm::Spline<double> sp(5, x, y, m0, m1); sp.setXYArrays(5,x,y,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp(xVec, yVec, m0, m1); sp.setXYContainers(xVec,yVec,m0, m1); testFull(sp, x, y, m0, m1);  };
    { Opm::Spline<double> sp; sp.setArrayOfPoints(5,points,m0, m1); testFull(sp, x, y, m0, m1); };
    { Opm::Spline<double> sp; sp.setContainerOfPoints(pointVec,m0, m1); testFull(sp, x, y, m0, m1);  };
#if GCC_VERSION >= 40500
    { Opm::Spline<double> sp; sp.setContainerOfTuples(pointsInitList,m0, m1); testFull(sp, x, y, m0, m1); };
#endif

    // natural spline
    { Opm::Spline<double> sp(5, x, y); sp.setXYArrays(5,x,y); testNatural(sp, x, y);  };
    { Opm::Spline<double> sp(xVec, yVec); sp.setXYContainers(xVec,yVec); testNatural(sp, x, y); };
    { Opm::Spline<double> sp; sp.setArrayOfPoints(5,points); testNatural(sp, x, y); };
    { Opm::Spline<double> sp; sp.setContainerOfPoints(pointVec); testNatural(sp, x, y); };
#if GCC_VERSION >= 40500
    { Opm::Spline<double> sp; sp.setContainerOfTuples(pointsInitList); testNatural(sp, x, y); };
#endif
}

void plot()
{
    const int numSamples = 5;
    const int n = numSamples - 1;
    typedef std::array<double, numSamples> FV;

    double x_[] = { 0, 5, 7.5, 8.75, 10 };
    double y_[] = { 10, 0, 10, 0, 10 };
    double m1 = 10;
    double m2 = -10;
    FV &xs = *reinterpret_cast<FV*>(x_);
    FV &ys = *reinterpret_cast<FV*>(y_);

    Opm::Spline<double> spFull(xs, ys, m1, m2);
    Opm::Spline<double> spNatural(xs, ys);
    Opm::Spline<double> spPeriodic(xs, ys, /*type=*/Opm::Spline<double>::Periodic);
    Opm::Spline<double> spMonotonic(xs, ys, /*type=*/Opm::Spline<double>::Monotonic);

    spFull.printCSV(x_[0] - 1.00001,
                    x_[n] + 1.00001,
                    1000);
    std::cout << "\n";
    spNatural.printCSV(x_[0] - 1.00001,
                       x_[n] + 1.00001,
                       1000);
    std::cout << "\n";

    spPeriodic.printCSV(x_[0] - 1.00001,
                        x_[n] + 1.00001,
                       1000);
    std::cout << "\n";

    spMonotonic.printCSV(x_[0] - 1.00001,
                         x_[n] + 1.00001,
                         1000);
    std::cout << "\n";

    std::cerr << "Spline is monotonic: " << spFull.monotonic(x_[0], x_[n]) << "\n";
}

int main(int argc, char** argv)
{
    try {
        testAll();

        plot();
    }
    catch (const std::exception &e) {
        std::cout << "Caught OPM exception: " << e.what() << "\n";
    }

    return 0;
}
