//===========================================================================
//                                                                           
// File: MiscibilityLiveGas.cpp                                               
//                                                                           
// Created: Wed Feb 10 09:21:53 2010                                         
//                                                                           
// Author: Bjørn Spjelkavik <bsp@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
//===========================================================================
/*
  Copyright 2010 SINTEF ICT, Applied Mathematics.

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


#include <algorithm>
#include "MiscibilityLiveGas.hpp"
#include <dune/common/ErrorMacros.hpp>
#include <dune/common/linInt.hpp>
#include <dune/common/Units.hpp>

using namespace std;
using namespace Dune;

namespace Opm
{

    //------------------------------------------------------------------------
    // Member functions
    //-------------------------------------------------------------------------

    /// Constructor
    MiscibilityLiveGas::MiscibilityLiveGas(const table_t& pvtg, const EclipseUnits& units)
    {
	// GAS, PVTG
	const int region_number = 0;
	if (pvtg.size() != 1) {
	    THROW("More than one PVD-region");
	}
	saturated_gas_table_.resize(4);
	const int sz = pvtg[region_number].size();
	for (int k=0; k<4; ++k) {
	    saturated_gas_table_[k].resize(sz);
	}
        using namespace Dune::unit;
        const double bunit = units.gasvol_r/units.gasvol_s;
        const double runit = units.liqvol_s/units.gasvol_s;
	for (int i=0; i<sz; ++i) {
	    saturated_gas_table_[0][i] = convert::from(pvtg[region_number][i][0], units.pressure); // p
	    saturated_gas_table_[1][i] = convert::from(pvtg[region_number][i][2], bunit); // Bg
	    saturated_gas_table_[2][i] = convert::from(pvtg[region_number][i][3], units.viscosity);  // mu_g
	    saturated_gas_table_[3][i] = convert::from(pvtg[region_number][i][1], runit); // Rv
	}

	undersat_gas_tables_.resize(sz);
	for (int i=0; i<sz; ++i) {
	    undersat_gas_tables_[i].resize(3);
	    int tsize = (pvtg[region_number][i].size() - 1)/3;
	    undersat_gas_tables_[i][0].resize(tsize);
	    undersat_gas_tables_[i][1].resize(tsize);
	    undersat_gas_tables_[i][2].resize(tsize);
	    for (int j=0, k=0; j<tsize; ++j) {
		undersat_gas_tables_[i][0][j] = convert::from(pvtg[region_number][i][++k], runit); // Rv
		undersat_gas_tables_[i][1][j] = convert::from(pvtg[region_number][i][++k], bunit); // Bg
		undersat_gas_tables_[i][2][j] = convert::from(pvtg[region_number][i][++k], units.viscosity); // mu_g
	    }
	}
    }

    // Destructor
     MiscibilityLiveGas::~MiscibilityLiveGas()
    {
    }

    double MiscibilityLiveGas::getViscosity(int /*region*/, double press, const surfvol_t& surfvol) const
    {
	return miscible_gas(press, surfvol, 2, false);
    }

    // Vaporised oil-gas ratio.
    double MiscibilityLiveGas::R(int /*region*/, double press, const surfvol_t& surfvol) const
    {
        if (surfvol[Liquid] == 0.0) {
            return 0.0;
        }
	double R = linearInterpolationExtrap(saturated_gas_table_[0],
					     saturated_gas_table_[3], press);
	double maxR = surfvol[Liquid]/surfvol[Vapour];
	if (R < maxR ) {  // Saturated case
	    return R;
	} else {
	    return maxR;  // Undersaturated case
	}
    }

    // Vaporised oil-gas ratio derivative
    double MiscibilityLiveGas::dRdp(int /*region*/, double press, const surfvol_t& surfvol) const
    {
	double R = linearInterpolationExtrap(saturated_gas_table_[0],
					     saturated_gas_table_[3], press);
	double maxR = surfvol[Liquid]/surfvol[Vapour];
	if (R < maxR ) {  // Saturated case
	    return linearInterpolDerivative(saturated_gas_table_[0],
					    saturated_gas_table_[3],
					    press);
	} else {
	    return 0.0;  // Undersaturated case
	}	
    }

    double MiscibilityLiveGas::B(int /*region*/, double press, const surfvol_t& surfvol) const
    {
        if (surfvol[Vapour] == 0.0) return 1.0; // To handle no-gas case.
        return  miscible_gas(press, surfvol, 1, false);
    }

    double MiscibilityLiveGas::dBdp(int /*region*/, double press, const surfvol_t& surfvol) const
    {	
        if (surfvol[Vapour] == 0.0) return 0.0; // To handle no-gas case.
        return miscible_gas(press, surfvol, 1, true);
    }

    double MiscibilityLiveGas::miscible_gas(double press, const surfvol_t& surfvol, int item,
					    bool deriv) const
    {
	int section;
	double R = linearInterpolationExtrap(saturated_gas_table_[0],
					     saturated_gas_table_[3], press,
					     section);
	double maxR = surfvol[Liquid]/surfvol[Vapour];
	if (deriv) {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolDerivative(saturated_gas_table_[0],
						saturated_gas_table_[item],
						press);
	    } else {  // Undersaturated case
		int is = section;
		if (undersat_gas_tables_[is][0].size() < 2) {
		    double val = (saturated_gas_table_[item][is+1]
				  - saturated_gas_table_[item][is]) /
			(saturated_gas_table_[0][is+1] -
			 saturated_gas_table_[0][is]);
		    return val;
		}
		double val1 =
		    linearInterpolationExtrap(undersat_gas_tables_[is][0],
					      undersat_gas_tables_[is][item],
					      maxR);
		double val2 = 
		    linearInterpolationExtrap(undersat_gas_tables_[is+1][0],
					      undersat_gas_tables_[is+1][item],
					      maxR);
		double val = (val2 - val1)/
		    (saturated_gas_table_[0][is+1] - saturated_gas_table_[0][is]);
		return val;
	    }
	} else {
	    if (R < maxR ) {  // Saturated case
		return linearInterpolationExtrap(saturated_gas_table_[0],
						 saturated_gas_table_[item],
						 press);
	    } else {  // Undersaturated case
		int is = section;
		// Extrapolate from first table section
		if (is == 0 && press < saturated_gas_table_[0][0]) {
		    return linearInterpolationExtrap(undersat_gas_tables_[0][0],
						     undersat_gas_tables_[0][item],
						     maxR);
		}

		// Extrapolate from last table section
		int ltp = saturated_gas_table_[0].size() - 1;
		if (is+1 == ltp && press > saturated_gas_table_[0][ltp]) {
		    return linearInterpolationExtrap(undersat_gas_tables_[ltp][0],
						     undersat_gas_tables_[ltp][item],
						     maxR);
		}

		// Interpolate between table sections
		double w = (press - saturated_gas_table_[0][is]) /
		    (saturated_gas_table_[0][is+1] - 
		     saturated_gas_table_[0][is]);
		if (undersat_gas_tables_[is][0].size() < 2) {
		    double val = saturated_gas_table_[item][is] +
			w*(saturated_gas_table_[item][is+1] -
			   saturated_gas_table_[item][is]);
		    return val;
		}
		double val1 =
		    linearInterpolationExtrap(undersat_gas_tables_[is][0],
					      undersat_gas_tables_[is][item],
					      maxR);
		double val2 = 
		    linearInterpolationExtrap(undersat_gas_tables_[is+1][0],
					      undersat_gas_tables_[is+1][item],
					      maxR);
		double val = val1 + w*(val2 - val1);
		return val;
	    }
	}
    }


} // namespace Opm
