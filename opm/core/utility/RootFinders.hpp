//===========================================================================
//
// File: RootFinders.hpp
//
// Created: Thu May  6 19:59:42 2010
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Jostein R Natvig    <jostein.r.natvig@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2010 Statoil ASA.

  This file is part of The Open Reservoir Simulator Project (OpenRS).

  OpenRS is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OpenRS is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OpenRS.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPENRS_ROOTFINDERS_HEADER
#define OPENRS_ROOTFINDERS_HEADER



#include <algorithm>
#include <limits>
#include <cmath>
#include <opm/core/utility/ErrorMacros.hpp>

namespace Opm
{

    struct ThrowOnError
    {
        static double handleBracketingFailure(const double x0, const double x1, const double f0, const double f1)
        {
            THROW("Error in parameters, zero not bracketed: [a, b] = ["
                  << x0 << ", " << x1 << "]    f(a) = " << f0 << "   f(b) = " << f1);
            return -1e100; // Never reached.
        }
        static double handleTooManyIterations(const double x0, const double x1, const int maxiter)
        {
            THROW("Maximum number of iterations exceeded: " << maxiter << "\n"
                  << "Current interval is [" << std::min(x0, x1) << ", "
                  << std::max(x0, x1) << "]");
        }
    };


    struct ContinueOnError
    {
    };



    template <class ErrorPolicy = ThrowOnError>
    class RegulaFalsi
    {
    public:


	/// Implements a modified regula falsi method as described in
	/// "Improved algorithms of Illinois-type for the numerical
	/// solution of nonlinear equations"
	/// by J. A. Ford.
	/// Current variant is the 'Pegasus' method.
        template <class Functor>
	inline static double solve(const Functor& f,
                            const double a,
                            const double b,
                            const int max_iter,
                            const double tolerance,
                            int& iterations_used)
	{
	    using namespace std;
	    const double macheps = numeric_limits<double>::epsilon();
	    const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);

	    double x0 = a;
	    double x1 = b;
	    double f0 = f(x0);
	    const double epsF = tolerance + macheps*max(fabs(f0), 1.0);
	    if (fabs(f0) < epsF) {
		return x0;
	    }
	    double f1 = f(x1);
	    if (fabs(f1) < epsF) {
		return x1;
	    }
	    if (f0*f1 > 0.0) {
                return ErrorPolicy::handleBracketingFailure(a, b, f0, f1);
	    }
	    iterations_used = 0;
	    // In every iteraton, x1 is the last point computed,
	    // and x0 is the last point computed that makes it a bracket.
	    while (fabs(x1 - x0) >= 1e-9*eps) {
		double xnew = regulaFalsiStep(x0, x1, f0, f1);
		double fnew = f(xnew);
// 		cout << "xnew = " << xnew << "    fnew = " << fnew << endl;
		++iterations_used;
		if (iterations_used > max_iter) {
                    return ErrorPolicy::handleTooManyIterations(x0, x1, max_iter);
		}
		if (fabs(fnew) < epsF) {
		    return xnew;
		}
		// Now we must check which point we must replace.
		if ((fnew > 0.0) == (f0 > 0.0)) {
		    // We must replace x0.
		    x0 = x1;
		    f0 = f1;
		} else {
		    // We must replace x1, this is the case where
		    // the modification to regula falsi kicks in,
		    // by modifying f0.
		    // 1. The classic Illinois method
//  		    const double gamma = 0.5;
		    // @afr: The next two methods do not work??!!?
		    // 2. The method called 'Method 3' in the paper.
// 		    const double phi0 = f1/f0;
// 		    const double phi1 = fnew/f1;
// 		    const double gamma = 1.0 - phi1/(1.0 - phi0);
		    // 3. The method called 'Method 4' in the paper.
// 		    const double phi0 = f1/f0;
// 		    const double phi1 = fnew/f1;
// 		    const double gamma = 1.0 - phi0 - phi1;
// 		    cout << "phi0 = " << phi0 <<" phi1 = " << phi1 <<
// 		    " gamma = " << gamma << endl;
		    // 4. The 'Pegasus' method
 		    const double gamma = f1/(f1 + fnew);
		    f0 *= gamma;
		}
		x1 = xnew;
		f1 = fnew;
	    }
	    return 0.5*(x0 + x1);
	}


	/// Implements a modified regula falsi method as described in
	/// "Improved algorithms of Illinois-type for the numerical
	/// solution of nonlinear equations"
	/// by J. A. Ford.
	/// Current variant is the 'Pegasus' method.
        /// This version takes an extra parameter for the initial guess.
	template <class Functor>
	inline static double solve(const Functor& f,
                            const double initial_guess,
                            const double a,
                            const double b,
                            const int max_iter,
                            const double tolerance,
                            int& iterations_used)
	{
	    using namespace std;
	    const double macheps = numeric_limits<double>::epsilon();
	    const double eps = tolerance + macheps*max(max(fabs(a), fabs(b)), 1.0);

            double f_initial = f(initial_guess);
	    const double epsF = tolerance + macheps*max(fabs(f_initial), 1.0);
            if (fabs(f_initial) < epsF) {
                return initial_guess;
            }
	    double x0 = a;
	    double x1 = b;
            double f0 = f_initial;
            double f1 = f_initial;
            if (x0 != initial_guess) {
                f0 = f(x0);
                if (fabs(f0) < epsF) {
                    return x0;
                }
	    }
            if (x1 != initial_guess) {
                f1 = f(x1);
                if (fabs(f1) < epsF) {
                    return x1;
                }
            }
            if (f0*f_initial < 0.0) {
                x1 = initial_guess;
                f1 = f_initial;
            } else {
                x0 = initial_guess;
                f0 = f_initial;
            }
	    if (f0*f1 > 0.0) {
		THROW("Error in parameters, zero not bracketed: [a, b] = ["
		      << a << ", " << b << "]    fa = " << f0 << "   fb = " << f1);
	    }
	    iterations_used = 0;
	    // In every iteraton, x1 is the last point computed,
	    // and x0 is the last point computed that makes it a bracket.
	    while (fabs(x1 - x0) >= 1e-9*eps) {
		double xnew = regulaFalsiStep(x0, x1, f0, f1);
		double fnew = f(xnew);
// 		cout << "xnew = " << xnew << "    fnew = " << fnew << endl;
		++iterations_used;
		if (iterations_used > max_iter) {
		    THROW("Maximum number of iterations exceeded.\n"
			  << "Current interval is [" << min(x0, x1) << ", "
			  << max(x0, x1) << "]");
		}
		if (fabs(fnew) < epsF) {
		    return xnew;
		}
		// Now we must check which point we must replace.
		if ((fnew > 0.0) == (f0 > 0.0)) {
		    // We must replace x0.
		    x0 = x1;
		    f0 = f1;
		} else {
		    // We must replace x1, this is the case where
		    // the modification to regula falsi kicks in,
		    // by modifying f0.
		    // 1. The classic Illinois method
//  		    const double gamma = 0.5;
		    // @afr: The next two methods do not work??!!?
		    // 2. The method called 'Method 3' in the paper.
// 		    const double phi0 = f1/f0;
// 		    const double phi1 = fnew/f1;
// 		    const double gamma = 1.0 - phi1/(1.0 - phi0);
		    // 3. The method called 'Method 4' in the paper.
// 		    const double phi0 = f1/f0;
// 		    const double phi1 = fnew/f1;
// 		    const double gamma = 1.0 - phi0 - phi1;
// 		    cout << "phi0 = " << phi0 <<" phi1 = " << phi1 <<
// 		    " gamma = " << gamma << endl;
		    // 4. The 'Pegasus' method
 		    const double gamma = f1/(f1 + fnew);
		    f0 *= gamma;
		}
		x1 = xnew;
		f1 = fnew;
	    }
	    return 0.5*(x0 + x1);
	}


    private:
	inline static double regulaFalsiStep(const double a,
                                             const double b,
                                             const double fa,
                                             const double fb)
	{
	    ASSERT(fa*fb < 0.0);
	    return (b*fa - a*fb)/(fa - fb);
	}


    };




	template <class Functor>
	inline void bracketZero(const Functor& f,
                                const double x0,
                                const double dx,
                                double& a,
                                double& b)
        {
            const int max_iters = 100;
            double f0 = f(x0);
            double cur_dx = dx;
            int i = 0;
            for (; i < max_iters; ++i) {
                double x = x0 + cur_dx;
                double f_new = f(x);
                if (f0*f_new <= 0.0) {
                    break;
                }
                cur_dx = -2.0*cur_dx;
            }
            if (i == max_iters) {
                THROW("Could not bracket zero in " << max_iters << "iterations.");
            }
            if (cur_dx < 0.0) {
                a = x0 + cur_dx;
                b = i < 2 ? x0 : x0 + 0.25*cur_dx;
            } else {
                a = i < 2 ? x0 : x0 + 0.25*cur_dx;
                b = x0 + cur_dx;
            }
        }


} // namespace Opm




#endif // OPENRS_ROOTFINDERS_HEADER
