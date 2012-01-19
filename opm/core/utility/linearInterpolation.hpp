//===========================================================================
//
// File: linearInterpolation.hpp
//
// Created: Tue Sep  9 12:49:39 2008
//
// Author(s): Atgeirr F Rasmussen <atgeirr@sintef.no>
//
// $Date$
//
// $Revision$
//
//===========================================================================

/*
  Copyright 2009, 2010 SINTEF ICT, Applied Mathematics.
  Copyright 2009, 2010 Statoil ASA.

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

#ifndef OPENRS_LINEARINTERPOLATION_HEADER
#define OPENRS_LINEARINTERPOLATION_HEADER


#include <vector>
#include <algorithm>

namespace Opm
{

    /** Linear interpolation. 
     *  Given an increasing vector xv of parameter values and
     *  a vector yv of point values of the same size,
     *  the function returns ...
     */
    template <typename T>
    T linearInterpolation(const std::vector<double>& xv,
			  const std::vector<T>& yv,
			  double x)
    {
	std::vector<double>::const_iterator lb = std::lower_bound(xv.begin(), xv.end(), x);
	int lb_ix = lb - xv.begin();
	if (lb_ix == 0) {
	    return yv[0];
	} else if (lb_ix == int(xv.size())) {
	    return yv.back();
	} else {
	    double w = (x - xv[lb_ix - 1])/(xv[lb_ix] - xv[lb_ix - 1]);
	    return (1.0 - w)*yv[lb_ix - 1] + w*yv[lb_ix];
	}
    }

    /// @brief
    /// @todo Doc me!
    /// @tparam
    /// @param
    /// @return
    template <typename T>
    T linearInterpolationDerivative(const std::vector<double>& xv,
				    const std::vector<T>& yv,
				    double x)
    {
	double epsilon = 1e-4; // @@ Ad hoc, should choose based on xv.
	double x_low = std::max(xv[0], x - epsilon);
	double x_high = std::min(xv.back(), x + epsilon);
	T low = linearInterpolation(xv, yv, x_low);
	T high = linearInterpolation(xv, yv, x_high);
	return (high - low)/(x_high - x_low);
    }


} // namespace Opm



#endif // OPENRS_LINEARINTERPOLATION_HEADER
