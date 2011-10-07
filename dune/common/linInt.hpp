//===========================================================================
//
// File: linInt.hpp
//
// Created: Tue Feb 16 14:44:10 2010
//
// Author(s): Bjørn Spjelkavik    <bsp@sintef.no>
//            Atgeirr F Rasmussen <atgeirr@sintef.no>
//            Bård Skaflestad     <bard.skaflestad@sintef.no>
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

#ifndef OPENRS_LININT_HEADER
#define OPENRS_LININT_HEADER


namespace Dune
{

    inline int tableIndex(const std::vector<double>& table, double x)
    {
	// Returns an index in an ordered table such that x is between
	// table[j] and table[j+1]. If x is out of range, first or last
	// interval is returned; Binary search.
	int n = table.size() - 1;
	if (n < 2) {
	    return 0;
	}
	int jl = 0;
	int ju = n;
	bool ascend = (table[n] > table[0]);
	while (ju - jl > 1) {
	    int jm = (ju + jl)/2;   // Compute a midpoint
	    if ( (x >= table[jm]) == ascend) {
		jl = jm;     // Replace lower limit
	    } else {
		ju = jm;     // Replace upper limit
	    }
	}
	return jl;
    }

    inline double linearInterpolation(const std::vector<double>& xv,
                                      const std::vector<double>& yv, double x)
    {
	// Returns end point if x is outside xv
	std::vector<double>::const_iterator lb = lower_bound(xv.begin(), xv.end(), x);
	int ix2 = lb - xv.begin();
	if (ix2 == 0) {
	    return yv[0];
	} else if (ix2 == int(xv.size())) {
	    return yv[ix2-1];
	}
	int ix1 = ix2 - 1;
	return  (yv[ix2] - yv[ix1])/(xv[ix2] - xv[ix1])*(x - xv[ix1]) + yv[ix1];
    }

    inline double linearInterpolDerivative(const std::vector<double>& xv,
                                           const std::vector<double>& yv, double x)
    {
	int ix1 = tableIndex(xv, x);
	int ix2 = ix1 + 1;
	return  (yv[ix2] - yv[ix1])/(xv[ix2] - xv[ix1]);
    }

    inline double linearInterpolationExtrap(const std::vector<double>& xv,
                                            const std::vector<double>& yv, double x)
    {
	// Extrapolates if x is outside xv
	int ix1 = tableIndex(xv, x);
	int ix2 = ix1 + 1;
	return  (yv[ix2] - yv[ix1])/(xv[ix2] - xv[ix1])*(x - xv[ix1]) + yv[ix1];
    }
    
    inline double linearInterpolationExtrap(const std::vector<double>& xv,
                                            const std::vector<double>& yv,
                                            double x, int& ix1)
    {
	// Extrapolates if x is outside xv
	ix1 = tableIndex(xv, x);
	int ix2 = ix1 + 1;
	return (yv[ix2] - yv[ix1])/(xv[ix2] - xv[ix1])*(x - xv[ix1]) + yv[ix1];
    }

} // namespace Dune





#endif // OPENRS_LININT_HEADER
